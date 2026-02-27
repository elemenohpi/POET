"""
CGP (Cartesian Genetic Programming) post-processing layer for POET.

Takes a POET model (patterns + weights) and learns a non-linear combination
of rule-match signals using CGP, replacing the default additive prediction
(ŷ = Σ w_j) with a learned function ŷ = f(w_1, w_2, …, w_n).

Based on Rachel Sarabia's integration work (Rachel2026 branch).

Usage as standalone:
    python cgp_runner.py --config config.ini

Usage from within POET (called automatically when use_cgp = True):
    from cgp_runner import run_cgp_on_model
    run_cgp_on_model(config, model_path)
"""

import os
import copy
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

import fitness as F
import individual as I

try:
    from cgp import CartesianGP
    from cgp.cgp_operators import add, sub, mul, div

    CGP_AVAILABLE = True
except ImportError:
    CGP_AVAILABLE = False


def rules_to_constants(
    config: dict,
    model_path: str,
    peptides: pd.DataFrame,
):
    """Convert a POET model + peptide dataset into CGP input/output arrays.

    Uses the POET Fitness evaluation engine to match rules properly,
    supporting all features (gaps, char classes, variable gaps, etc.).

    For each peptide, build a feature vector where each element is the
    rule's weight if the rule matched the sequence, else 0.

    Parameters
    ----------
    config : dict
        POET config dict (needed for matching mode / experimental flags).
    model_path : str
        Path to the POET model CSV.
    peptides : DataFrame
        Must have a sequence column (col 0) and a fitness/score column (col 1).

    Returns
    -------
    input_data : np.ndarray, shape (n_peptides, n_rules)
    target_output : np.ndarray, shape (n_peptides,)
    poet_predictions : np.ndarray, shape (n_peptides,)
        The raw additive POET predictions for each peptide.
    """
    cols = peptides.columns.tolist()
    seq_col = cols[0]
    fit_col = cols[1]

    # Build a Fitness evaluator and an Individual from the model file
    fitness_obj = F.Fitness(config)
    individual = I.Individual(config)
    individual.makeFromFile(model_path)
    active_indices = [
        i for i, r in enumerate(individual.rules) if isinstance(r.pattern, str)
    ]
    n_rules = len(individual.rules)

    input_data = []
    target_output = []
    poet_predictions = []

    for _, row in peptides.iterrows():
        target_val = float(row[fit_col])
        sequence = str(row[seq_col])
        target_output.append(target_val)

        # Reset the individual for a fresh evaluation
        fitness_obj.resetIndividual(individual)
        # Get the POET prediction (this uses the proper matching engine)
        _, prediction = fitness_obj.eval(sequence, target_val, individual, True)
        poet_predictions.append(prediction)

        # Build the per-rule feature vector: weight if matched, else 0
        datum = []
        for rule in individual.rules:
            if rule.status == 1:
                datum.append(float(rule.weight))
            else:
                datum.append(0.0)
        input_data.append(datum)

    return (
        np.array(input_data),
        np.array(target_output),
        np.array(poet_predictions),
    )


def setup_cgp(n_inputs: int, seed: int, config: dict | None = None):
    """Construct and return a CartesianGP evolver.

    Parameters
    ----------
    n_inputs : int
        Number of input features (= number of active rules in the model).
    seed : int
        Random seed for reproducibility.
    config : dict, optional
        POET config dict.  CGP hyperparameters are read from keys prefixed
        with ``cgp_`` when present, otherwise sensible defaults are used.
    """
    if not CGP_AVAILABLE:
        raise ImportError(
            "CGP package not found. Clone it from "
            "https://github.com/MarkKocherovsky/cgp_crossover (branch With-SMAC) "
            "and place the src/cgp/ directory inside the POET project root."
        )

    config = config or {}

    max_generations = int(config.get("cgp_generations", "1000"))
    model_size = int(config.get("cgp_model_size", "32"))
    max_parents = int(config.get("cgp_parents", "1"))
    max_children = int(config.get("cgp_children", "4"))
    mutation_type = config.get("cgp_mutation_type", "full")
    mutation_rate = float(config.get("cgp_mutation_rate", "1.0"))
    selection_type = config.get("cgp_selection_type", "paretoelite")
    fitness_function = config.get("cgp_fitness_function", "correlation")
    n_elites = int(config.get("cgp_n_elites", "1"))
    tournament_size = int(config.get("cgp_tournament_size", "4"))

    model_parameters = {
        "max_size": model_size,
        "inputs": n_inputs,
        "outputs": 1,
        "arity": 2,
        "constants": np.array([]),
    }
    function_bank = {"add": add, "sub": sub, "mul": mul, "div": div}

    asex = max_parents < max_children
    mutation_breeding = asex

    checkpoint_path = os.path.join(os.environ.get("SCRATCH", "output"), "ckpt")
    Path(checkpoint_path).mkdir(parents=True, exist_ok=True)
    checkpoint_file = os.path.join(checkpoint_path, f"poet_cgp_ckpt_{seed}.pkl")

    evolver = CartesianGP(
        parents=max_parents,
        children=max_children,
        max_generations=max_generations,
        mutation=mutation_type,
        selection=selection_type,
        xover=None,
        fixed_length=True,
        fitness_function=fitness_function,
        model_parameters=model_parameters,
        n_points=1,
        n_elites=n_elites,
        tournament_size=tournament_size,
        function_bank=function_bank,
        mutation_breeding=mutation_breeding,
        checkpoint_filename=checkpoint_file,
        one_dimensional_xover=False,
        seed=seed,
        tuning=False,
    )
    return evolver


def run_cgp_on_model(config: dict, model_path: str | None = None):
    """Run CGP post-processing on a POET model.

    Parameters
    ----------
    config : dict
        POET config dictionary.
    model_path : str, optional
        Path to the POET model CSV.  If None, uses ``config["output_model"]``.

    Returns
    -------
    best_train_model, best_test_model
        The best CGP models found on training and test data respectively.
    """
    if not CGP_AVAILABLE:
        raise ImportError(
            "CGP package not found. See cgp_runner.py docstring for installation."
        )

    model_path = model_path or config["output_model"]
    data_path = config["learn_data"]
    unseen_path = config.get("unseen_data", data_path)
    seed = int(config["seed"])
    step_size = int(config.get("cgp_step_size", "10"))

    print(f"\n--- CGP Post-Processing ---")
    print(f"Model: {model_path}")
    print(f"Training data: {data_path}")
    print(f"Test data: {unseen_path}")

    data = pd.read_csv(data_path)
    unseen_data = pd.read_csv(unseen_path)

    # Convert rules ↔ peptides to feature vectors using POET's matching engine
    input_data, target_fitness, poet_train_preds = rules_to_constants(
        config, model_path, data
    )
    input_test, test_fitness, poet_test_preds = rules_to_constants(
        config, model_path, unseen_data
    )

    if input_data.shape[1] == 0:
        print("WARNING: No active rules in model. Skipping CGP.")
        return None, None

    evolver = setup_cgp(input_data.shape[1], seed, config)

    # Output directory for CGP results
    output_dir = os.path.dirname(model_path) or "output"
    cgp_dir = os.path.join(output_dir, "cgp")
    Path(cgp_dir).mkdir(parents=True, exist_ok=True)

    print(
        f"Running CGP evolution ({config.get('cgp_generations', '1000')} generations)..."
    )
    best_train, best_test = evolver.fit(
        input_data, input_test, target_fitness, test_fitness, step_size=step_size
    )

    evolver.save_metrics(cgp_dir)
    print("---")
    print("Best test model:")
    best_test.print_model()
    print(f"Complexity: {best_test.count_active_nodes()}")

    # Save CGP model
    df = pd.DataFrame(best_test.model)
    df.to_csv(os.path.join(cgp_dir, "best_cgp_model.csv"), index=True)
    print(f"CGP model saved to {cgp_dir}/best_cgp_model.csv")

    # ── Compare POET (additive) vs CGP (non-linear) ──
    comparison = evaluate_poet_vs_cgp(
        poet_train_preds,
        poet_test_preds,
        best_train,
        best_test,
        input_data,
        target_fitness,
        input_test,
        test_fitness,
    )
    # Save comparison to CSV
    comp_df = pd.DataFrame([comparison])
    comp_path = os.path.join(cgp_dir, "poet_vs_cgp.csv")
    comp_df.to_csv(comp_path, index=False)
    print(f"Comparison saved to {comp_path}")

    return best_train, best_test


def evaluate_poet_vs_cgp(
    poet_train_preds: np.ndarray,
    poet_test_preds: np.ndarray,
    cgp_train_model,
    cgp_test_model,
    cgp_train_inputs: np.ndarray,
    train_targets: np.ndarray,
    cgp_test_inputs: np.ndarray,
    test_targets: np.ndarray,
) -> dict:
    """Compare POET additive predictions vs CGP non-linear predictions.

    Parameters
    ----------
    poet_train_preds : np.ndarray
        POET predictions on training data (from ``rules_to_constants``).
    poet_test_preds : np.ndarray
        POET predictions on test data (from ``rules_to_constants``).
    cgp_train_model, cgp_test_model
        Best CGP models from ``evolver.fit()``.
    cgp_train_inputs, cgp_test_inputs : np.ndarray
        Rule-match feature matrices for CGP.
    train_targets, test_targets : np.ndarray
        Actual fitness values.

    Returns
    -------
    dict
        Metrics for both models and deltas.
    """

    # Align POET predictions (linear fit, same as POET uses internally)
    def align_predictions(preds, targets):
        if np.std(preds) < 1e-8:
            return preds  # can't align constant predictions
        slope, intercept = np.polyfit(preds, targets, 1)
        return preds * slope + intercept

    poet_train_aligned = align_predictions(poet_train_preds, train_targets)
    poet_test_aligned = align_predictions(poet_test_preds, test_targets)

    # ── CGP predictions: ŷ = f(w_1, …, w_n) ──
    cgp_train_preds = cgp_train_model(cgp_train_inputs, mutable=False).flatten()
    cgp_test_preds = cgp_test_model(cgp_test_inputs, mutable=False).flatten()

    # ── Metrics ──
    def compute_metrics(preds, targets):
        mask = np.isfinite(preds) & np.isfinite(targets)
        p, t = preds[mask], targets[mask]
        if len(p) < 2 or np.std(p) < 1e-8 or np.std(t) < 1e-8:
            return 0.0, float("inf")
        r, _ = stats.pearsonr(p, t)
        r2 = r**2
        rmse = np.sqrt(np.mean((p - t) ** 2))
        return round(r2, 6), round(rmse, 6)

    poet_train_r2, poet_train_rmse = compute_metrics(poet_train_aligned, train_targets)
    poet_test_r2, poet_test_rmse = compute_metrics(poet_test_aligned, test_targets)
    cgp_train_r2, cgp_train_rmse = compute_metrics(cgp_train_preds, train_targets)
    cgp_test_r2, cgp_test_rmse = compute_metrics(cgp_test_preds, test_targets)

    delta_train_r2 = round(cgp_train_r2 - poet_train_r2, 6)
    delta_test_r2 = round(cgp_test_r2 - poet_test_r2, 6)
    delta_train_rmse = round(cgp_train_rmse - poet_train_rmse, 6)
    delta_test_rmse = round(cgp_test_rmse - poet_test_rmse, 6)

    # ── Print comparison table ──
    print("\n+------------------+---------------+---------------+---------+")
    print("|           POET (additive) vs CGP (non-linear)              |")
    print("+------------------+---------------+---------------+---------+")
    print("|     Metric       |     POET      |      CGP      |  Delta  |")
    print("+------------------+---------------+---------------+---------+")
    print(
        f"|  Train r2        |  {poet_train_r2:<12} |  {cgp_train_r2:<12} | {delta_train_r2:+.4f} |"
    )
    print(
        f"|  Test  r2        |  {poet_test_r2:<12} |  {cgp_test_r2:<12} | {delta_test_r2:+.4f} |"
    )
    print(
        f"|  Train RMSE      |  {poet_train_rmse:<12} |  {cgp_train_rmse:<12} | {delta_train_rmse:+.4f} |"
    )
    print(
        f"|  Test  RMSE      |  {poet_test_rmse:<12} |  {cgp_test_rmse:<12} | {delta_test_rmse:+.4f} |"
    )
    print("+------------------+---------------+---------------+---------+")

    # Interpretation
    if delta_test_r2 > 0.01:
        print(f"  -> CGP IMPROVED test r2 by {delta_test_r2:+.4f}")
    elif delta_test_r2 < -0.01:
        print(f"  -> CGP WORSENED test r2 by {delta_test_r2:+.4f}")
    else:
        print(f"  -> CGP had NEGLIGIBLE effect on test r2 ({delta_test_r2:+.4f})")

    return {
        "poet_train_r2": poet_train_r2,
        "poet_test_r2": poet_test_r2,
        "poet_train_rmse": poet_train_rmse,
        "poet_test_rmse": poet_test_rmse,
        "cgp_train_r2": cgp_train_r2,
        "cgp_test_r2": cgp_test_r2,
        "cgp_train_rmse": cgp_train_rmse,
        "cgp_test_rmse": cgp_test_rmse,
        "delta_train_r2": delta_train_r2,
        "delta_test_r2": delta_test_r2,
        "delta_train_rmse": delta_train_rmse,
        "delta_test_rmse": delta_test_rmse,
    }


# ── Standalone CLI ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    import eletility

    parser = argparse.ArgumentParser(
        description="Run CGP post-processing on a POET model"
    )
    parser.add_argument(
        "--config",
        default="config.ini",
        help="Path to config.ini file (default: config.ini)",
    )
    parser.add_argument(
        "--model",
        default=None,
        help="Path to POET model CSV (overrides config output_model)",
    )
    args = parser.parse_args()

    configparser = eletility.ConfigParser()
    cfg = configparser.read(args.config)
    run_cgp_on_model(cfg, args.model)

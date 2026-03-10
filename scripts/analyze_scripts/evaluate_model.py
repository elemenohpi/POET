"""Evaluate a POET model against a dataset.

Reports:
  - Pearson r, r², and fitness (1 - r²)
  - Per-sequence predicted fitness values (aligned to real scale via linear
    regression, i.e. the actual predicted values, not the raw scores)

Usage:
  python scripts/analyze_scripts/evaluate_model.py <model_csv> [options]

Examples:
  # Use default config (config.ini) and its learn_data
  python scripts/analyze_scripts/evaluate_model.py output/Mar-02-2026-upar_epoch_0/models/cleaned/model_40.csv

  # Specify a config
  python scripts/analyze_scripts/evaluate_model.py output/.../model_40.csv -c configs/exp_charclass_vargaps.ini

  # Override the dataset
  python scripts/analyze_scripts/evaluate_model.py output/.../model_40.csv -d data/unseen.csv

  # Save predictions to CSV
  python scripts/analyze_scripts/evaluate_model.py output/.../model_40.csv -o predictions.csv
"""

import argparse
import os
import sys
import math

import numpy as np
import pandas as pd
from scipy import stats

# Allow imports from the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import eletility
import fitness as F
import individual as I


def evaluate(model_path, config, data_path=None):
    """Return correlation stats and aligned predictions for a single model.

    Parameters
    ----------
    model_path : str
        Path to the model CSV file.
    config : dict
        Parsed POET config dictionary.
    data_path : str or None
        Optional override for the dataset CSV (must have sequence, fitness columns).

    Returns
    -------
    dict with keys:
        pearson_r, r_squared, fitness, sequences, actual, raw_predictions,
        aligned_predictions, align_slope, align_intercept
    """
    # Optionally override dataset
    if data_path is not None:
        config = dict(config)
        config["learn_data"] = data_path

    # Load model
    ind = I.Individual(config)
    ind.makeFromFile(model_path)

    # Build fitness evaluator
    objF = F.Fitness(config)

    # Get raw predictions
    objF.resetIndividual(ind)
    raw_predictions = []
    for seq, actual in zip(objF.sequences, objF.fitness_values):
        _, prediction = objF.eval(seq, actual, ind, True)
        raw_predictions.append(prediction)

    actuals = objF.fitness_values

    # Pearson correlation
    r, p_value = stats.pearsonr(raw_predictions, actuals)
    if math.isnan(r):
        r = 0.0
    r_squared = r ** 2
    fitness = 1 - r_squared

    # Linear alignment: map raw predictions to the real fitness scale
    align = np.polyfit(raw_predictions, actuals, 1)
    slope, intercept = align[0], align[1]
    aligned_predictions = [p * slope + intercept for p in raw_predictions]

    return {
        "pearson_r": r,
        "r_squared": r_squared,
        "fitness": fitness,
        "p_value": p_value,
        "sequences": objF.sequences,
        "actual": actuals,
        "raw_predictions": raw_predictions,
        "aligned_predictions": aligned_predictions,
        "align_slope": slope,
        "align_intercept": intercept,
        "num_rules": len(ind.rules),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate a POET model against a dataset."
    )
    parser.add_argument("model", help="Path to the model CSV file.")
    parser.add_argument(
        "-c", "--config", default="config.ini",
        help="Path to POET config file (default: config.ini)."
    )
    parser.add_argument(
        "-d", "--data", default=None,
        help="Override dataset CSV (must have sequence, fitness columns)."
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Save per-sequence predictions to this CSV file."
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="Only print correlation stats, skip per-sequence table."
    )
    args = parser.parse_args()

    configparser = eletility.ConfigParser()
    config = configparser.read(args.config)

    result = evaluate(args.model, config, args.data)

    # Print correlation summary
    print(f"Model:       {args.model}")
    print(f"Rules:       {result['num_rules']}")
    print(f"Dataset:     {args.data or config.get('learn_data', '?')} "
          f"({len(result['sequences'])} sequences)")
    print(f"Pearson r:   {result['pearson_r']:.6f}")
    print(f"r²:          {result['r_squared']:.6f}")
    print(f"Fitness:     {result['fitness']:.6f}  (1 - r²)")
    print(f"p-value:     {result['p_value']:.2e}")
    print(f"Alignment:   predicted_fitness = raw * {result['align_slope']:.6f} "
          f"+ {result['align_intercept']:.6f}")

    # Build predictions dataframe
    df = pd.DataFrame({
        "sequence": result["sequences"],
        "actual_fitness": result["actual"],
        "predicted_fitness": [round(v, 6) for v in result["aligned_predictions"]],
        "raw_score": [round(v, 6) for v in result["raw_predictions"]],
    })

    if not args.quiet:
        print()
        print(df.to_string(index=False))

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nPredictions saved to {args.output}")


if __name__ == "__main__":
    main()

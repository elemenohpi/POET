"""Experiment mode: compare two configurations across multiple seeds.

Usage (via CLI):
    python poet.py --experiment configA.ini configB.ini --replicates 10 --gens 100
"""

import copy
import os
import random as rand
import time

import archivist
import optimizer
import pop as population


def run_experiment(
    config_a: dict, config_b: dict, replicates: int, gens: int, workers: int
):
    """Run both configs for *replicates* seeds × *gens* generations each.

    Collects best fitness per generation for every replicate, then prints
    a summary table comparing the two configurations.
    """
    label_a = config_a.get("matching_mode", "substring")
    label_b = config_b.get("matching_mode", "substring")
    # Disambiguate labels if they happen to be the same
    if label_a == label_b:
        label_a += "_A"
        label_b += "_B"

    print("=" * 70)
    print("EXPERIMENT MODE")
    print(
        "  Config A : {} ({})".format(
            label_a, config_a.get("matching_mode", "substring")
        )
    )
    print(
        "  Config B : {} ({})".format(
            label_b, config_b.get("matching_mode", "substring")
        )
    )
    print("  Replicates : {}".format(replicates))
    print("  Generations: {}".format(gens))
    print("  Workers    : {}".format(workers if workers else "auto"))
    print("=" * 70)

    results_a: list[list[float]] = []  # [replicate][generation] = best_fitness
    results_b: list[list[float]] = []

    base_seed = int(config_a.get("seed", "1"))

    for rep in range(replicates):
        seed = base_seed + rep
        print("\n--- Replicate {}/{} (seed={}) ---".format(rep + 1, replicates, seed))

        best_a = _run_single(config_a, label_a, seed, gens, workers)
        results_a.append(best_a)

        best_b = _run_single(config_b, label_b, seed, gens, workers)
        results_b.append(best_b)

    # ── Summary ─────────────────────────────────────────────────────────
    _print_summary(label_a, label_b, results_a, results_b, gens, replicates)
    _save_csv(label_a, label_b, results_a, results_b, gens, replicates)
    _plot(label_a, label_b, results_a, results_b, gens)


def _run_single(
    config: dict, label: str, seed: int, gens: int, workers: int
) -> list[float]:
    """Run one configuration for one seed. Returns list of best fitness per gen."""
    cfg = copy.deepcopy(config)
    cfg["seed"] = str(seed)
    cfg["runs"] = str(gens)
    cfg["workers"] = str(workers)

    # Unique output paths so runs don't clobber each other
    os.makedirs("output/experiment", exist_ok=True)
    cfg["output_evo"] = "output/experiment/evo_{}_{}.csv".format(label, seed)
    cfg["output_model"] = "output/experiment/model_{}_{}.csv".format(label, seed)

    rand.seed(int(cfg["seed"]))

    arch = archivist.Archivist(cfg)
    arch.setup()

    t0 = time.time()
    pop = population.Population(cfg)
    opt = optimizer.Optimizer(cfg, pop)
    opt.optimize()
    elapsed = time.time() - t0

    # Extract best-fitness-per-generation from the evo log
    best_per_gen = _parse_evo(cfg["output_evo"])
    print(
        "  [{}] seed={} done in {:.1f}s  final_best={:.4f}".format(
            label, seed, elapsed, best_per_gen[-1] if best_per_gen else float("nan")
        )
    )
    return best_per_gen


def _parse_evo(path: str) -> list[float]:
    """Read the evo CSV and return the best_fitness column (col index 1)."""
    values: list[float] = []
    with open(path, "r") as f:
        for line_no, line in enumerate(f):
            if line_no == 0:
                continue  # skip header
            parts = line.strip().split(",")
            if len(parts) >= 2:
                try:
                    values.append(float(parts[1]))
                except ValueError:
                    pass
    return values


def _print_summary(la: str, lb: str, ra, rb, gens: int, reps: int):
    """Print a table of mean best fitness at key generation milestones."""
    import statistics

    print("\n" + "=" * 70)
    print("SUMMARY  ({} replicates × {} generations)".format(reps, gens))
    print("=" * 70)
    print("{:>6s}  {:>14s}  {:>14s}  {:>10s}".format("Gen", la, lb, "Winner"))
    print("-" * 50)

    # Show milestones: every 10% of gens, plus first and last
    milestones = sorted(set([0, gens - 1] + [int(gens * p / 10) for p in range(1, 10)]))

    for g in milestones:
        vals_a = [r[g] for r in ra if g < len(r)]
        vals_b = [r[g] for r in rb if g < len(r)]
        if not vals_a or not vals_b:
            continue
        mean_a = statistics.mean(vals_a)
        mean_b = statistics.mean(vals_b)
        winner = la if mean_a < mean_b else lb if mean_b < mean_a else "tie"
        print("{:>6d}  {:>14.4f}  {:>14.4f}  {:>10s}".format(g, mean_a, mean_b, winner))

    # Final stats
    final_a = [r[-1] for r in ra if r]
    final_b = [r[-1] for r in rb if r]
    mean_a = statistics.mean(final_a)
    mean_b = statistics.mean(final_b)
    std_a = statistics.stdev(final_a) if len(final_a) > 1 else 0
    std_b = statistics.stdev(final_b) if len(final_b) > 1 else 0

    print("-" * 50)
    print("Final mean ± std:")
    print("  {}: {:.4f} ± {:.4f}".format(la, mean_a, std_a))
    print("  {}: {:.4f} ± {:.4f}".format(lb, mean_b, std_b))
    overall = la if mean_a < mean_b else lb
    print("  Overall winner: {}".format(overall))
    print("=" * 70)


def _save_csv(la: str, lb: str, ra, rb, gens: int, reps: int):
    """Save raw results to a CSV for further analysis."""
    path = "output/experiment/experiment_results.csv"
    with open(path, "w") as f:
        f.write("replicate,generation,{}_best,{}_best\n".format(la, lb))
        for rep in range(reps):
            max_g = min(len(ra[rep]), len(rb[rep]))
            for g in range(max_g):
                f.write("{},{},{},{}\n".format(rep, g, ra[rep][g], rb[rep][g]))
    print("\nRaw results saved to {}".format(path))


def _plot(la: str, lb: str, ra, rb, gens: int):
    """Plot mean best fitness ± std error for both configs and save as PNG."""
    import statistics

    import matplotlib.pyplot as plt

    generations = list(range(gens))

    def _stats(results):
        means, errs = [], []
        for g in generations:
            vals = [r[g] for r in results if g < len(r)]
            if vals:
                m = statistics.mean(vals)
                se = (statistics.stdev(vals) / len(vals) ** 0.5) if len(vals) > 1 else 0
                means.append(m)
                errs.append(se)
            else:
                means.append(float("nan"))
                errs.append(0)
        return means, errs

    means_a, errs_a = _stats(ra)
    means_b, errs_b = _stats(rb)

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(
        generations,
        means_a,
        yerr=errs_a,
        label=la,
        capsize=3,
        elinewidth=1,
        markeredgewidth=1,
        errorevery=max(1, gens // 20),
    )
    ax.errorbar(
        generations,
        means_b,
        yerr=errs_b,
        label=lb,
        capsize=3,
        elinewidth=1,
        markeredgewidth=1,
        errorevery=max(1, gens // 20),
    )

    ax.set_xlabel("Generation")
    ax.set_ylabel("Best Fitness (lower is better)")
    ax.set_title("{} vs {} ({} replicates)".format(la, lb, len(ra)))
    ax.legend()
    ax.grid(True, alpha=0.3)

    os.makedirs("output/experiment", exist_ok=True)
    out_path = "output/experiment/experiment_plot.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Plot saved to {}".format(out_path))

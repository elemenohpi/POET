"""Experiment mode: compare two or more configurations across multiple seeds.

Usage (via CLI):
    python poet.py --experiment configA.ini configB.ini ... --replicates 10 --gens 50
"""

import copy
import os
import random as rand
import time

import archivist
import optimizer
import pop as population


# ── Colour palette for up to 10 lines ──────────────────────────────────────
_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]


def run_experiment(configs: list[dict], replicates: int, gens: int, workers: int):
    """Run all configs for *replicates* seeds × *gens* generations each.

    Collects best fitness per generation for every replicate, then prints
    a summary table and saves a multi-line comparison plot.
    """
    # ── Derive labels ───────────────────────────────────────────────────
    labels: list[str] = []
    for i, cfg in enumerate(configs):
        lbl = cfg.get("experiment_label", "")
        if not lbl:
            lbl = cfg.get("matching_mode", "substring")
            lbl += "_{}".format(i)
        labels.append(lbl)

    # Ensure uniqueness
    seen: dict[str, int] = {}
    for i, lbl in enumerate(labels):
        if lbl in seen:
            seen[lbl] += 1
            labels[i] = "{}_{}".format(lbl, seen[lbl])
        else:
            seen[lbl] = 0

    print("=" * 70)
    print("EXPERIMENT MODE")
    for i, (lbl, cfg) in enumerate(zip(labels, configs)):
        print("  Config {:>2d} : {}".format(i + 1, lbl))
    print("  Replicates : {}".format(replicates))
    print("  Generations: {}".format(gens))
    print("  Workers    : {}".format(workers if workers else "auto"))
    print("=" * 70)

    # results[config_idx][replicate] = list[float] per generation
    all_results: list[list[list[float]]] = [[] for _ in configs]

    base_seed = int(configs[0].get("seed", "1"))

    for rep in range(replicates):
        seed = base_seed + rep
        print("\n--- Replicate {}/{} (seed={}) ---".format(rep + 1, replicates, seed))

        for ci, (cfg, lbl) in enumerate(zip(configs, labels)):
            best = _run_single(cfg, lbl, seed, gens, workers)
            all_results[ci].append(best)

    # ── Summary ─────────────────────────────────────────────────────────
    _print_summary(labels, all_results, gens, replicates)
    _save_csv(labels, all_results, gens, replicates)
    _plot(labels, all_results, gens)


def _run_single(
    config: dict, label: str, seed: int, gens: int, workers: int
) -> list[float]:
    """Run one configuration for one seed. Returns list of best fitness per gen."""
    cfg = copy.deepcopy(config)
    cfg["seed"] = str(seed)
    cfg["runs"] = str(gens)
    cfg["workers"] = str(workers)

    # Sanitise label for filenames
    safe_label = label.replace(" ", "_").replace("/", "-")

    # Unique output paths so runs don't clobber each other
    os.makedirs("output/experiment", exist_ok=True)
    cfg["output_evo"] = "output/experiment/evo_{}_{}.csv".format(safe_label, seed)
    cfg["output_model"] = "output/experiment/model_{}_{}.csv".format(safe_label, seed)

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


def _print_summary(labels: list[str], all_results, gens: int, reps: int):
    """Print a table of mean best fitness at key generation milestones."""
    import statistics

    print("\n" + "=" * 70)
    print("SUMMARY  ({} replicates x {} generations)".format(reps, gens))
    print("=" * 70)

    # Column widths
    cw = max(14, max(len(l) + 2 for l in labels))
    header = "{:>6s}".format("Gen")
    for lbl in labels:
        header += "  {:>{w}s}".format(lbl, w=cw)
    header += "  {:>12s}".format("Best")
    print(header)
    print("-" * len(header))

    milestones = sorted(set([0, gens - 1] + [int(gens * p / 10) for p in range(1, 10)]))

    for g in milestones:
        row = "{:>6d}".format(g)
        means = []
        for ci, lbl in enumerate(labels):
            vals = [r[g] for r in all_results[ci] if g < len(r)]
            if vals:
                m = statistics.mean(vals)
                means.append((m, lbl))
                row += "  {:>{w}.4f}".format(m, w=cw)
            else:
                means.append((float("inf"), lbl))
                row += "  {:>{w}s}".format("N/A", w=cw)
        # "Best" = lowest RMSE
        best_lbl = min(means, key=lambda x: x[0])[1]
        row += "  {:>12s}".format(best_lbl)
        print(row)

    # Final stats per config
    print("-" * len(header))
    print("Final generation mean +/- std:")
    final_stats = []
    for ci, lbl in enumerate(labels):
        finals = [r[-1] for r in all_results[ci] if r]
        m = statistics.mean(finals) if finals else float("nan")
        s = statistics.stdev(finals) if len(finals) > 1 else 0
        final_stats.append((m, s, lbl))
        print("  {}: {:.4f} +/- {:.4f}".format(lbl, m, s))
    overall = min(final_stats, key=lambda x: x[0])
    print("  Overall best: {}".format(overall[2]))
    print("=" * 70)


def _save_csv(labels: list[str], all_results, gens: int, reps: int):
    """Save raw results to a CSV for further analysis."""
    path = "output/experiment/experiment_results.csv"
    with open(path, "w") as f:
        header = "replicate,generation," + ",".join("{}_best".format(l) for l in labels)
        f.write(header + "\n")
        for rep in range(reps):
            max_g = min(len(all_results[ci][rep]) for ci in range(len(labels)))
            for g in range(max_g):
                vals = ",".join(
                    str(all_results[ci][rep][g]) for ci in range(len(labels))
                )
                f.write("{},{},{}\n".format(rep, g, vals))
    print("\nRaw results saved to {}".format(path))


def _plot(labels: list[str], all_results, gens: int):
    """Plot mean best fitness +/- std error for all configs and save as PNG."""
    import statistics

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    generations = list(range(gens))

    def _stats(results):
        means, errs = [], []
        for g in generations:
            vals = [r[g] for r in results if g < len(r)]
            if vals:
                m = statistics.mean(vals)
                se = statistics.stdev(vals) / len(vals) ** 0.5 if len(vals) > 1 else 0
                means.append(m)
                errs.append(se)
            else:
                means.append(float("nan"))
                errs.append(0)
        return means, errs

    fig, ax = plt.subplots(figsize=(12, 7))

    for ci, lbl in enumerate(labels):
        means, errs = _stats(all_results[ci])
        color = _COLORS[ci % len(_COLORS)]
        ax.errorbar(
            generations,
            means,
            yerr=errs,
            label=lbl,
            color=color,
            capsize=2,
            elinewidth=0.8,
            markeredgewidth=0.8,
            errorevery=max(1, gens // 15),
            linewidth=1.5,
        )

    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("Best Fitness (lower is better)", fontsize=12)
    ax.set_title(
        "Experimental Feature Comparison ({} replicates x {} gens)".format(
            len(all_results[0]), gens
        ),
        fontsize=14,
    )
    ax.legend(fontsize=9, loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    os.makedirs("output/experiment", exist_ok=True)
    out_path = "output/experiment/experiment_plot.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Plot saved to {}".format(out_path))

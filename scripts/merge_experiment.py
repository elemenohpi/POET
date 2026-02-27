"""Merge existing experiment evo files and regenerate the combined CSV + plot.

Reads all evo_{label}_{seed}.csv files from output/experiment/, groups by label,
rebuilds the combined experiment_results.csv and experiment_plot.png.
"""

import os
import re
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EVO_DIR = "output/experiment"

# Ordered labels to include (matching the safe-label used in filenames)
LABELS_ORDER = [
    "Gaps_Only",
    "+_Char_Classes",
    "+_Variable_Gaps",
    "+_Weighted_Pos",
    "+_Composition",
    "+_Match_Count",
    "+_Circular",
    "All_Features",
    "+_CC_&_VarGaps",
    "+_CC_&_VG_&_OLS",
]

# Display labels (nicer for plot legend / CSV header)
DISPLAY = {
    "Gaps_Only": "Gaps Only",
    "+_Char_Classes": "+ Char Classes",
    "+_Variable_Gaps": "+ Variable Gaps",
    "+_Weighted_Pos": "+ Weighted Pos",
    "+_Composition": "+ Composition",
    "+_Match_Count": "+ Match Count",
    "+_Circular": "+ Circular",
    "All_Features": "All Features",
    "+_CC_&_VarGaps": "+ CC & VarGaps",
    "+_CC_&_VG_&_OLS": "+ CC & VG & OLS",
}

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


def _parse_evo(path: str) -> list[float]:
    """Read evo CSV, return best_fitness column (index 1)."""
    vals: list[float] = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.strip().split(",")
            if len(parts) >= 2:
                try:
                    vals.append(float(parts[1]))
                except ValueError:
                    pass
    return vals


def main():
    # ── Discover evo files ──────────────────────────────────────────────
    # filename pattern:  evo_{label}_{seed}.csv
    pattern = re.compile(r"^evo_(.+)_(\d+)\.csv$")

    groups: dict[str, list[tuple[int, list[float]]]] = {}

    for fname in sorted(os.listdir(EVO_DIR)):
        m = pattern.match(fname)
        if not m:
            continue
        label, seed = m.group(1), int(m.group(2))
        if label not in LABELS_ORDER:
            continue
        vals = _parse_evo(os.path.join(EVO_DIR, fname))
        groups.setdefault(label, []).append((seed, vals))

    # Order
    ordered_labels = [l for l in LABELS_ORDER if l in groups]
    display_labels = [DISPLAY.get(l, l) for l in ordered_labels]

    gens = min(len(v) for l in ordered_labels for _, v in groups[l])
    reps = min(len(groups[l]) for l in ordered_labels)

    print(f"Labels : {display_labels}")
    print(f"Gens   : {gens}")
    print(f"Reps   : {reps}")

    # ── Build all_results[config_idx][rep] = list[float] ────────────────
    all_results: list[list[list[float]]] = []
    for lbl in ordered_labels:
        sorted_reps = sorted(groups[lbl], key=lambda x: x[0])[:reps]
        all_results.append([v[:gens] for _, v in sorted_reps])

    # ── Save combined CSV ───────────────────────────────────────────────
    csv_path = os.path.join(EVO_DIR, "experiment_results.csv")
    with open(csv_path, "w") as f:
        header = "replicate,generation," + ",".join(
            f"{dl}_best" for dl in display_labels
        )
        f.write(header + "\n")
        for rep in range(reps):
            for g in range(gens):
                vals = ",".join(
                    str(all_results[ci][rep][g]) for ci in range(len(ordered_labels))
                )
                f.write(f"{rep},{g},{vals}\n")
    print(f"Saved  : {csv_path}")

    # ── Plot ────────────────────────────────────────────────────────────
    generations = list(range(gens))

    def _stats(results):
        means, errs = [], []
        for g in generations:
            vals = [r[g] for r in results if g < len(r)]
            m = statistics.mean(vals)
            se = statistics.stdev(vals) / len(vals) ** 0.5 if len(vals) > 1 else 0
            means.append(m)
            errs.append(se)
        return means, errs

    fig, ax = plt.subplots(figsize=(12, 7))

    for ci, dl in enumerate(display_labels):
        means, errs = _stats(all_results[ci])
        color = _COLORS[ci % len(_COLORS)]
        ax.errorbar(
            generations,
            means,
            yerr=errs,
            label=dl,
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
        f"Experimental Feature Comparison ({reps} replicates × {gens} gens)",
        fontsize=14,
    )
    ax.legend(fontsize=9, loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    plot_path = os.path.join(EVO_DIR, "experiment_plot.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot   : {plot_path}")

    # ── Summary table ───────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"SUMMARY  ({reps} replicates × {gens} generations)")
    print(f"{'='*70}")
    cw = max(14, max(len(dl) + 2 for dl in display_labels))
    header = f"{'Gen':>6s}" + "".join(f"  {dl:>{cw}s}" for dl in display_labels)
    header += f"  {'Best':>12s}"
    print(header)
    print("-" * len(header))
    milestones = sorted(set([0, gens - 1] + [int(gens * p / 10) for p in range(1, 10)]))
    for g in milestones:
        row = f"{g:>6d}"
        means = []
        for ci, dl in enumerate(display_labels):
            vals = [r[g] for r in all_results[ci] if g < len(r)]
            m = statistics.mean(vals)
            means.append((m, dl))
            row += f"  {m:>{cw}.4f}"
        best_lbl = min(means, key=lambda x: x[0])[1]
        row += f"  {best_lbl:>12s}"
        print(row)
    print("-" * len(header))
    print("Final generation mean ± std:")
    final_stats = []
    for ci, dl in enumerate(display_labels):
        finals = [r[-1] for r in all_results[ci]]
        m = statistics.mean(finals)
        s = statistics.stdev(finals) if len(finals) > 1 else 0
        final_stats.append((m, s, dl))
        print(f"  {dl}: {m:.4f} ± {s:.4f}")
    overall = min(final_stats, key=lambda x: x[0])
    print(f"  Overall best: {overall[2]}")
    print("=" * 70)


if __name__ == "__main__":
    main()

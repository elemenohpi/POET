"""Figure + paired statistics for the mutation-selection ablation.

Run after:
    python poet.py --experiment configs/ablation_a_independent.ini \
        configs/ablation_b_feasible.ini configs/ablation_c_feasible_matched.ini \
        --replicates 10 --gens 100

Reads output/experiment/ and writes output/experiment/ablation.png.
"""

import csv
import glob
import os
import re
import statistics as st

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXP = "output/experiment"

# Chart chrome (light surface).
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
# Categorical slots 1-3: validated all-pairs for CVD separation, light mode.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]

SHORT = {
    "mut_A_independent": "A  independent",
    "mut_B_feasible": "B  feasible",
    "mut_C_feasible_matched": "C  feasible, rate-matched",
}


def spread(values, min_gap):
    """Nudge near-identical label positions apart, preserving their order."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = list(values)
    for k in range(1, len(order)):
        prev, cur = order[k - 1], order[k]
        if out[cur] - out[prev] < min_gap:
            out[cur] = out[prev] + min_gap
    return out


def paired(x, y):
    """Return (mean diff, SE, t) for paired samples x - y."""
    d = [xi - yi for xi, yi in zip(x, y)]
    if len(d) < 2:
        return st.mean(d), float("nan"), float("nan")
    se = st.stdev(d) / len(d) ** 0.5
    return st.mean(d), se, (st.mean(d) / se if se else float("nan"))


# ── Load ────────────────────────────────────────────────────────────────────
rows = list(csv.DictReader(open(os.path.join(EXP, "experiment_results.csv"))))
arms = [k[:-5] for k in rows[0] if k.endswith("_best")]
gens = sorted({int(r["generation"]) for r in rows})
reps = sorted({int(r["replicate"]) for r in rows})

series = {a: {g: [] for g in gens} for a in arms}
for r in rows:
    for a in arms:
        series[a][int(r["generation"])].append(float(r[a + "_best"]))
finals = {a: series[a][gens[-1]] for a in arms}

# Rule counts come from the evo log, not the model CSV: the model CSV only
# holds rules with status != 0, so counting its rows gives USED rules and
# misses the total that parsimony_pressure actually penalises.
# evo columns: gen, best fitness, best test, best rule count, best used rules, ...
_seed_re = re.compile(r"_(\d+)\.csv$")  # skips the sibling *_mutstats.csv files
seeds = sorted(
    int(m.group(1))
    for m in (
        _seed_re.search(p)
        for p in glob.glob(os.path.join(EXP, "evo_{}_*.csv".format(arms[0])))
    )
    if m
)


def last_evo_row(arm, seed):
    path = os.path.join(EXP, "evo_{}_{}.csv".format(arm, seed))
    data = [l.strip().split(",") for l in open(path) if l.strip()][1:]
    return data[-1]


total_rules = {a: [int(last_evo_row(a, s)[3]) for s in seeds] for a in arms}
used_rules = {a: [int(last_evo_row(a, s)[4]) for s in seeds] for a in arms}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    }
)
fig, (ax1, ax2, ax3) = plt.subplots(
    1, 3, figsize=(17.5, 5.8), gridspec_kw={"width_ratios": [1.75, 1, 1]}
)


def style(ax):
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=10)


# ── Panel 1: mean trajectory ± 1 SE ─────────────────────────────────────────
ends = [st.mean(series[a][gens[-1]]) for a in arms]
allm = [st.mean(series[a][g]) for a in arms for g in gens]
label_y = spread(ends, (max(allm) - min(allm)) * 0.075)

for i, a in enumerate(arms):
    means = [st.mean(series[a][g]) for g in gens]
    ses = [st.stdev(series[a][g]) / len(series[a][g]) ** 0.5 for g in gens]
    ax1.fill_between(
        gens,
        [m - s for m, s in zip(means, ses)],
        [m + s for m, s in zip(means, ses)],
        color=SERIES[i],
        alpha=0.13,
        linewidth=0,
    )
    ax1.plot(gens, means, color=SERIES[i], linewidth=2, solid_capstyle="round")
    ax1.plot([gens[-1]], [means[-1]], "o", color=SERIES[i], markersize=8, zorder=5)
    # Direct label (relief for the sub-3:1 slot): ink text, colored mark.
    ax1.annotate(
        SHORT[a].split("  ")[0],
        (gens[-1], means[-1]),
        xytext=(gens[-1] * 1.03, label_y[i]),
        textcoords="data",
        va="center",
        fontsize=11,
        color=INK,
        fontweight="bold",
        arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.8, shrinkA=2),
    )

ax1.set_xlabel("Generation", fontsize=11, color=INK2)
ax1.set_ylabel("Best fitness  (lower is better)", fontsize=11, color=INK2)
ax1.set_title(
    "Mean best fitness +/- 1 SE   ·   {} replicates".format(len(reps)),
    fontsize=12.5,
    color=INK,
    pad=12,
    loc="left",
)
ax1.set_xlim(0, gens[-1] * 1.11)
ax1.grid(True, color=GRID, linewidth=0.8)
style(ax1)
ax1.legend(
    handles=[
        plt.Line2D([], [], color=SERIES[i], linewidth=2, label=SHORT[a])
        for i, a in enumerate(arms)
    ],
    fontsize=10,
    loc="upper right",
    frameon=True,
    facecolor=SURFACE,
    edgecolor=GRID,
    labelcolor=INK2,
)


# ── Panels 2 & 3: paired per-seed outcomes ──────────────────────────────────
def paired_panel(ax, data, title, ylabel, fmt):
    xs = list(range(len(arms)))
    for j in range(len(seeds)):
        ax.plot(
            xs,
            [data[a][j] for a in arms],
            color=AXIS,
            linewidth=1.1,
            zorder=1,
            alpha=0.85,
        )
    for i, a in enumerate(arms):
        ax.plot(
            [xs[i]] * len(data[a]),
            data[a],
            "o",
            color=SERIES[i],
            markersize=8,
            zorder=3,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
        )
        m = st.mean(data[a])
        ax.plot([xs[i] - 0.2, xs[i] + 0.2], [m, m], color=INK, linewidth=2.4, zorder=4)
        ax.annotate(
            fmt.format(m),
            (xs[i] + 0.22, m),
            fontsize=10,
            va="center",
            color=INK,
            fontweight="bold",
        )
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT[a].split("  ")[0] for a in arms], fontsize=11, color=INK)
    ax.set_xlim(-0.45, len(arms) - 0.25)
    ax.set_ylabel(ylabel, fontsize=11, color=INK2)
    ax.set_title(title, fontsize=12.5, color=INK, pad=12, loc="left")
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    style(ax)


paired_panel(
    ax2,
    finals,
    "Final best fitness   ·   bar = mean",
    "Final best fitness",
    "{:.3f}",
)
paired_panel(
    ax3,
    used_rules,
    "Used rules in best model   ·   bar = mean",
    "Rules with status != 0",
    "{:.1f}",
)

fig.suptitle(
    "Mutation-selection ablation - {} replicates x {} generations".format(
        len(reps), gens[-1] + 1
    ),
    fontsize=14,
    color=INK,
    x=0.006,
    ha="left",
    y=0.985,
)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(EXP, "ablation.png")
fig.savefig(out, dpi=160, facecolor=SURFACE)
print("saved {}".format(out))

# ── Paired statistics ───────────────────────────────────────────────────────
crit = {5: 2.78, 10: 2.26, 20: 2.09}.get(len(seeds), 2.26)
for name, data, fmt in [
    ("FINAL BEST FITNESS", finals, "{:+.4f}"),
    ("TOTAL RULES  (what parsimony_pressure penalises)", total_rules, "{:+.2f}"),
    ("USED RULES  (status != 0)", used_rules, "{:+.2f}"),
]:
    print("\n{}  ({} paired seeds)".format(name, len(seeds)))
    for a in arms:
        print(
            "  {:<28s} mean {:.4f}   sd {:.4f}".format(
                SHORT[a], st.mean(data[a]), st.stdev(data[a])
            )
        )
    for x, y in [(1, 0), (2, 0), (2, 1)]:
        m, se, t = paired(data[arms[x]], data[arms[y]])
        wins = sum(
            1 for i in range(len(seeds)) if data[arms[x]][i] < data[arms[y]][i]
        )
        print(
            "  {} vs {}:  diff {}  SE {:.4f}  t({}) = {:+.2f}  {}  "
            "[first lower in {}/{}]".format(
                arms[x][4],
                arms[y][4],
                fmt.format(m),
                se,
                len(seeds) - 1,
                t,
                "SIGNIFICANT" if abs(t) > crit else "n.s.",
                wins,
                len(seeds),
            )
        )
print("\n  (|t| > {:.2f} needed for p < .05, two-tailed)".format(crit))

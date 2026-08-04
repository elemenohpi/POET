# Mutation-selection ablation

Tests whether POET should determine which mutations are *possible* before
rolling, instead of letting an impossible mutation waste its roll.

## The question

Every generation, each rule rolls an independent Bernoulli per operator. Some
of those rolls hit an operator that structurally cannot change the rule:

- `mut_add_to_pattern` on a rule already at `maximum_rule_size`
- `mut_insert_gap` on a rule shorter than 3 tokens
- `mut_remove_gap` on a rule that contains no gaps
- `mut_add_rule` on an individual already at `maximum_rule_count`

Those rolls return without doing anything. The proposal is to compute the
feasible operator set first and spend the roll on something that can land.

**This is not a free optimisation.** Removing wasted rolls means more
mutations actually apply, so the effective mutation rate goes up. The rates in
`config.ini` were tuned with that attrition already baked in. A three-arm
comparison is needed to separate "smarter operator choice" from "simply more
mutation".

## Config parameters

| Parameter | Values | Meaning |
|---|---|---|
| `mutation_selection` | `independent` (default), `feasible` | `independent` is the historical behaviour, bit-for-bit. `feasible` computes the feasible operator set per rule first, then redistributes the infeasible operators' probability mass proportionally over the feasible ones. |
| `mutation_rate_scale` | float, default `1.0` | Multiplies every mutation rate. Used to build the rate-matched control arm. |
| `mutation_stats` | `True` / `False` (default) | Counts attempted vs. changed per operator. Prints a table at the end of the run and writes `<output_evo>_mutstats.csv`. |

Both `feasible` and `mutation_stats` require `matching_mode = substring`; they
raise a `ValueError` in regex mode rather than silently doing nothing.

`mutation_selection = independent`, `mutation_rate_scale = 1.0`,
`mutation_stats = False` reproduces pre-existing runs exactly — verified by
diffing evo and model CSVs against the previous commit.

## The three arms

| Arm | Config | Selection | Scale |
|---|---|---|---|
| A — baseline | `configs/ablation_a_independent.ini` | `independent` | 1.0 |
| B — treatment | `configs/ablation_b_feasible.ini` | `feasible` | 1.0 |
| C — rate-matched control | `configs/ablation_c_feasible_matched.ini` | `feasible` | 0.7906 |

Arm C is arm B with the rates scaled down so the number of mutations that
actually land matches arm A.

**How to read the result:**

- B beats A **and** B beats C → the feasibility logic itself helps. Adopt it.
- B beats A **but** B ≈ C → all the change did was raise the mutation rate.
  Skip the refactor and edit the rates in `config.ini` instead.
- B ≈ A → no effect. Keep the simpler code.

## Running it

### Step 1 — measure the waste (already done for the chimera config)

```bash
python poet.py -config configs/ablation_a_independent.ini
```

100 generations, seed 333, on `data/chimera_multihance_relative.csv`:

```
operator                 attempted   changed    wasted   wasted%
------------------------------------------------------------------
change_weight               157318    156589       729      0.5%
add_to_pattern               93807     46867     46940     50.0%
change_character             91352     91352         0      0.0%
remove_from_pattern          62564     62564         0      0.0%
insert_gap                   38062      8631     29431     77.3%
remove_gap                   37940     15873     22067     58.2%
add_rule                      2554       116      2438     95.5%
remove_rule                   1510      1510         0      0.0%
------------------------------------------------------------------
TOTAL                       485107    383502    101605     20.9%
```

So ~1 roll in 5 does nothing, and the waste is concentrated in four
operators. `add_rule` is the worst: the population saturates
`maximum_rule_count = 80` within a few generations and 95.5% of add-rule rolls
are dead thereafter. The `change_weight` 0.5% is unrelated to feasibility —
`round(uniform(0, 1), 2)` rounds to 0.00 about that often.

The run prints the value to use for arm C:

```
Rate-matched control arm: set mutation_rate_scale = 0.7906
```

Re-measure and update `configs/ablation_c_feasible_matched.ini` if you change
the mutation rates, `maximum_rule_size`, `maximum_rule_count`, or the dataset.

### Step 2 — run the comparison

```bash
python poet.py --experiment \
  configs/ablation_a_independent.ini \
  configs/ablation_b_feasible.ini \
  configs/ablation_c_feasible_matched.ini \
  --replicates 5 --gens 100
```

Seeds are `333..337` (base seed + replicate index), identical across arms, so
the comparison is paired. Runtime is about 75 s per run, so 3 arms × 5
replicates ≈ 20 minutes single-threaded. Add `-w 0` to use all cores.

Outputs land in `output/experiment/`:

- `experiment_results.csv` — best fitness per generation per replicate per arm
- `experiment_plot.png` — mean ± standard error curves
- `evo_<label>_<seed>_mutstats.csv` — per-operator waste for every run
- a summary table printed at the end

### Step 3 — check the manipulation worked

Before reading the fitness result, confirm the arms did what they claim.
Compare the `TOTAL changed` row across arms' `_mutstats.csv` files:

Measured at 100 generations, seed 333:

- Arm A — 383,502 applied, 20.9% wasted
- Arm B — 480,015 applied, 0.7% wasted (**25% more mutation than A**)
- Arm C — 382,217 applied, 0.7% wasted (matched to A within 0.3%)

If arm C's applied count is not close to arm A's, the scale is stale —
re-measure and re-run.

## Result (10 replicates × 100 generations, seeds 333–342)

Figure: `output/experiment/ablation.png` (regenerate with
`python scripts/plot_ablation.py`).

| Arm | Final best fitness | 1 − r² | parsimony | mean r | used rules |
|---|---|---|---|---|---|
| A independent | **0.2592** ± 0.048 | 0.1814 | 0.0778 | 0.905 | 39.2 |
| B feasible | 0.2932 ± 0.052 | 0.2174 | 0.0758 | 0.885 | 42.1 |
| C feasible, rate-matched | 0.2655 ± 0.053 | 0.1957 | 0.0698 | 0.897 | 37.0 |

Paired t-tests (|t| > 2.26 for p < .05):

```
fitness      B vs A  +0.0340  t(9) = +1.56  n.s.   [B better in 4/10 seeds]
fitness      C vs A  +0.0063  t(9) = +0.25  n.s.   [C better in 5/10 seeds]
used rules   B vs A   +2.90   t(9) = +2.37  SIG    [B lower in 1/10 seeds]
```

**Verdict: do not adopt the change.**

- **C ≈ A.** Once the mutation-rate increase is controlled for, feasibility-aware
  selection is indistinguishable from the current code — a 5/5 seed split and
  t = 0.25. The selection logic itself contributes nothing.
- **B is directionally worse**, and the whole deficit is in fit quality, not
  parsimony: mean r drops 0.905 → 0.885 while its parsimony penalty is
  *lower* than A's. Not significant on its own (p ≈ 0.15), but it points the
  same way as the prior experience that removing "wasted" variation hurts POET.
- **B's models are both more complex and worse**: +2.9 used rules (the only
  significant result in the study, 9/10 seeds) at a worse correlation.

Read together: B ≈ C plus extra mutation, and C ≈ A, so the only thing the
refactor actually does is raise the effective mutation rate — which is
achievable by editing `config.ini` and does not help here.

Confirming B's harm at 80% power would need ≈ 32 seeds. That is not worth
running: C ≈ A already establishes there is no upside to find.

Two incidental observations from the runs:

- Total rule count sits at 70–78 against `maximum_rule_count = 80` in every
  arm, which is why `mut_add_rule` wastes 95.5% of its rolls. Raising the cap
  is a much cheaper lever than changing the mutation scheme.
- Arm C on seed 340 collapsed to a 5-rule model (3 used) scoring 0.233 —
  identical to arm A's 78-rule model on the same seed, with r = 0.879 vs
  0.919. One seed, not a reproducible effect, but a reminder that
  `parsimony_pressure = 0.001` makes 73 rules worth only 0.073 of fitness.

## Caveats

**Residual waste in `feasible` mode is expected.** The feasible set is
computed once per rule, before any operator is applied — the literal form of
the proposal. An operator applied early in the pass can still invalidate a
later one, which is where the residual ~0.6% comes from. The counters measure
this rather than hiding it.

**Feasibility predicates are structural, not exact.** They mirror the
early-return branches in the `mut_*` methods. An operator reported as feasible
can still pick a sub-action that happens to be a no-op (this is why
`mut_char_class` and `mut_variable_gap` cannot be made fully waste-free
without also making their internal action choice feasibility-aware). Those two
operators are disabled (rate 0.0) in the chimera config.

**The RNG stream differs between arms.** Any change to how many `R.random()`
calls fire desynchronizes everything downstream, so arms are only comparable
statistically across replicates, never run-for-run.

**Redistribution changes the operator mix, not just the total.** At
`maximum_rule_count`, `add_rule`'s probability mass flows to `remove_rule`; at
`maximum_rule_size`, `add_to_pattern`'s mass flows to `remove_from_pattern`
and `change_character`. Arm B therefore applies more shrink pressure at the
boundaries than arm A does. Watch the `-arc` (average rule count) column in
the evo log, not only fitness — this is a real distributional side effect and
`parsimony_pressure` interacts with it.

## Related: a separate bug this measurement exposed

In token mode, `mut_add_to_pattern` returns early when the rule is at
`maximum_rule_size` ([optimizer.py](optimizer.py) — the `len(tokens) >=
self.ruleSize` guard). That guard sits *before* the branch selection, so it
also kills the 70% "resample to a real observed motif" path, not just the
add-a-token path. A rule sitting at 4/4 tokens can never be swapped for a
different 4-token observed motif; it can only drift one token at a time via
`mut_change_character`.

That is a genuine loss of a distinct move type and is worth fixing on its own,
independently of how this ablation turns out. It is deliberately **not** fixed
here, so it does not confound the comparison.

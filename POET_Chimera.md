# POET-Chimera

POET-Chimera is a variant of POET that predicts and ranks **chimeric protein/DNA
constructs** built by swapping regions between two parent variants (labelled
`B2` and `B3` — e.g. the two species backbones, 1B2 and 1B3). It learns, from a
small set of measured constructs, a rule-based model that correlates
region-swap patterns with a measured fitness value (MultiHance ICP signal), then
uses that model to score **any-to-any** candidate swaps and rank the most
promising ones.

It lives on the `POET-Chimera` branch. Everything here is layered on top of the
standard POET evolutionary engine that lives on the `POET2.0` branch — the
chimera work adds a new *sequence representation* and a *data/ranking pipeline*
around that engine, but does **not** change the core genetic algorithm.

---

## 1. The problem it solves

A chimera construct is one of the two parent backbones (B2 or B3) with one or
more of its 12 regions replaced by the corresponding region from the other
parent. Each construct has a measured lab value (metal-binding / contrast
signal). Only a small number of constructs have been measured:

- the two natives (all-B2, all-B3), and
- the single-region swaps (each of the 12 slots swapped one at a time, on each
  backbone).

The goal is to predict the fitness of **combinations of swaps that were never
measured** (e.g. "swap slots 6 and 3 together") and rank them, so the lab can
prioritise which multi-swap constructs to synthesise next.

---

## 2. Sequence representation: position-aware tokens

This is the central idea that distinguishes POET-Chimera. Standard POET treats a
sequence as a **string of amino-acid characters**. POET-Chimera instead treats
each construct as a fixed list of **12 tokens**, one per target slot. Each token
records both *where* it sits (the target slot) and *what* fills it (the donor
species and the donor region):

```
S06:B2_09
└┬┘ └┬┘ └┬┘
 │   │   └── donor region 9
 │   └────── donor species B2
 └────────── target slot 6
```

A full construct is a space-separated string of 12 such tokens, always in slot
order 1..12. A native B2 backbone is:

```
S01:B2_01 S02:B2_02 ... S12:B2_12
```

and a single swap (B2 backbone, slot 6 receives B3's region 6) is:

```
S01:B2_01 ... S06:B3_06 ... S12:B2_12
```

Because the slot index is baked into every token, POET's motif matcher becomes
**position-aware for free**: the token `S06:B3_06` can only ever match at slot 6,
so a learned motif like "`S06:B3_06`" means "region 6 came from B3", which is
exactly the biological event we want to reason about. A single alphabet of
`12 slots × 2 species × 12 regions = 288` tokens covers every possible swap; the
24 tokens where `donor_region == target_slot` are the "native-position" swaps.

The codecs live in:
- `chimera_codec.py` — build/parse tokens, encode a construct from a backbone +
  list of swaps, decode tokens back to structured swaps, enumerate the alphabet,
  and generate any-to-any candidate constructs.
- `sequence_codec.py` — the generic token-mode plumbing (splitting/joining
  tokens, token-level pattern matching with `_` wildcards, drawing rule seeds
  from observed training motifs, slot-aware replacement pools).

---

## 3. The three-stage pipeline

### Stage 1 — Build the dataset (`scripts/chimera/build_chimera_dataset.py`)

Reads the raw workbook `chimera_init_data/ChimeraICPdata.xlsx` directly (no
Excel dependency — it parses the `.xlsx` XML), pulls the MultiHance replicate
measurements for both the B2 and B3 sheets, and:

1. maps each construct name (`B2`, `B3`, `TM6`, `T11`, …) to a 12-token sequence
   via `chimera_codec`;
2. averages replicates per construct;
3. normalises fitness as a **ratio to the B2 native baseline within each sheet**
   (`fitness = mean_raw / B2_mean_raw`), so native B2 = 1.0;
4. writes the POET training CSV (`data/chimera_multihance_relative.csv`), the
   token alphabet (`chimera_init_data/chimera_token_alphabet.csv`), a manifest
   describing every token, and a list of skipped/ambiguous constructs.

The training CSV's first two columns are `sequence` (the 12 tokens) and
`fitness` (the normalised value) — the format POET already expects.

### Stage 2 — Evolve a model (`python poet.py`)

Runs the standard POET evolutionary loop in **token mode** (configured by
`config.ini` / `configs/chimera_multihance.ini`). It evolves a population of
rule-sets; each rule is a token motif with a weight, and a model's prediction for
a construct is the **sum of the weights of the rules that match it**. Fitness of
a model is how well its predictions correlate with the measured values across the
training set (see §4). The best model is written to `output/chimera/model.csv`
and the per-generation log to `output/chimera/evo.csv`.

### Stage 3 — Rank candidates (`scripts/chimera/rank_chimera_candidates.py`)

Loads a trained model and:

1. **Calibrates** it: runs the model over the training constructs and fits a
   single linear map `predicted = raw_score * slope + intercept` (least squares)
   so raw additive scores land in real fitness units.
2. **Enumerates** any-to-any candidates for a chosen backbone and swap count
   (`chimera_codec.generate_candidates`) — every combination of target slots ×
   donor regions, up to `--max-swaps`.
3. **Scores** each candidate with the model, applies the linear calibration, and
   sorts descending. Emits the top-N as CSV.

Example:

```powershell
python scripts/chimera/rank_chimera_candidates.py output/chimera/model.csv \
  -c config.ini --backbone B3 --max-swaps 2 --top 50 -o output/chimera/ranked_B3.csv
```

---

## 4. The evolutionary algorithm (shared POET core)

POET-Chimera reuses POET's genetic-programming engine unchanged; only the
matching/mutation code branches into token-aware paths. The core loop
(`optimizer.py`) is:

- **Individuals & rules.** An individual is a set of `Rule`s. A rule is a
  `(pattern, weight, status)` triple. `pattern` is a motif (a token motif here);
  `weight ∈ [rule_weight_min, rule_weight_max]`; `status` flags whether the rule
  matched anything on the last evaluation.
- **Prediction.** For a construct, scan every position; each rule that matches
  contributes its weight once (first match wins per position). The model's raw
  prediction is the sum of matched weights. Token matching (`_eval_token`) is
  positional and treats `_` as a single-token wildcard.
- **Fitness (model quality).** With `fitness_alg = correlation` (the chimera
  default), a model's fitness is `1 - r²` between its predictions and the
  measured values (lower is better; a perfect linear fit → 0). `RMSE` is the
  alternative. A small `parsimony_pressure` term adds `pressure × rule_count` to
  discourage bloat.
- **Selection.** Tournament selection (`tournament_size`), with optional
  **diversity selection** that penalises a second parent for sharing too many
  patterns with the first (`diversity_weight`), to keep the population varied.
- **Crossover.** Cluster crossover: the offspring inherits the *used* (expressed)
  rules of both parents, plus each *unused* rule with probability
  `crossover_unused_selection_chance`; the result is trimmed back to
  `maximum_rule_count`, preferring to drop unused rules.
- **Mutation.** Per-generation operators add/remove rules, change weights, and
  grow/shrink/edit patterns. In token mode these are slot-aware (see §5).
- **Elitism.** The best individual is copied unchanged into the next generation,
  and reverted if a mutation made the elite worse.

---

## 5. What POET-Chimera changes vs. the POET2.0 algorithm

POET-Chimera is a strict **superset** of the POET2.0 branch: it contains all of
POET2.0's code and adds a token-mode representation plus the chimera pipeline.
When `sequence_mode = char` (the POET2.0 default) the two behave identically; the
chimera behaviour is switched on by `sequence_mode = token`.

| Aspect | POET2.0 | POET-Chimera |
|---|---|---|
| Unit of a sequence | Amino-acid **characters** in a string | **Tokens**, one per target slot (`S06:B2_09`), space-separated |
| Alphabet | ~20 amino-acid letters (`data/translation/...`) | 288 position/species/region tokens (`chimera_token_alphabet.csv`) |
| Position awareness | None — a motif can match anywhere | Built-in — slot index is part of each token, so a motif matches only at its slot |
| Matching | Substring scan over characters, gaps via `_` | Token-list scan (`_eval_token`), gaps via `_` token |
| Reverse matching | On by default | **Off** by default in token mode (`allow_reverse_match = False`) — slot order is meaningful, so a reversed motif is meaningless |
| Rule seeding | Random characters from the alphabet | Drawn from **observed training motifs** (`random_observed_token_pattern`); random only as fallback |
| Point mutation | Any character → any other amino acid | **Slot-constrained** — a token is replaced only by another token for the *same slot* (`token_replacement_pool`) |
| Gap mutations | Character ↔ `_` | Token ↔ `_`, edge wildcards stripped, native slot pools respected |
| Regex / char-class / variable-gap / position-weight / composition / CGP experimental features | Available | Force-disabled in token mode (they assume character semantics) |
| Input data | A peptide fitness CSV | Built from the chimera workbook by `build_chimera_dataset.py`, B2-baseline normalised |
| Output use | A fitness model | A fitness model **plus** an any-to-any candidate ranker with linear calibration |
| Core GA (selection, crossover, elitism, correlation fitness) | — | **Unchanged** |

New files added on the chimera branch:

```
chimera_codec.py                              token/construct/candidate codec
sequence_codec.py                             generic token-mode plumbing
scripts/chimera/build_chimera_dataset.py      workbook → training CSV + alphabet
scripts/chimera/rank_chimera_candidates.py    model → ranked any-to-any candidates
tests/verify_chimera_token_mode.py            end-to-end smoke test
configs/chimera_multihance.ini                chimera config (also the root config.ini default)
configs/default_peptide.ini                   preserved POET2.0-style peptide default
chimera_init_data/                            source workbook, alphabet, manifest
data/chimera_multihance_relative.csv          generated training data
```

Modified files (`fitness.py`, `individual.py`, `optimizer.py`, `archivist.py`)
only add token-mode branches guarded by `SC.is_token_mode(config)`; the
character-mode paths are untouched.

---

## 6. Key configuration (`config.ini` / `configs/chimera_multihance.ini`)

```ini
sequence_mode = token                 # switch that turns on the chimera behaviour
token_delimiter = space
alphabet_data = chimera_init_data/chimera_token_alphabet.csv
alphabet_column = code
allow_reverse_match = False           # slot order matters — no reversed motifs
learn_data = data/chimera_multihance_relative.csv

fitness_alg = correlation             # model fitness = 1 - r²
maximum_rule_size = 4                 # motifs span up to 4 slots
maximum_rule_count = 80
population_size = 100
runs = 100
tournament_size = 5
diversity_selection = True
parsimony_pressure = 0.001
matching_mode = substring
enable_gaps = True
```

All the character-level experimental flags (`exp_char_classes`,
`exp_variable_gaps`, `use_cgp`, …) are present but disabled — they are ignored in
token mode.

---

## 7. How to run

```powershell
# 1. Build the tokenized training data from the workbook
python scripts/chimera/build_chimera_dataset.py

# 2. Train a chimera model (uses config.ini by default)
python poet.py

# 3. Rank any-to-any candidates from the trained model
python scripts/chimera/rank_chimera_candidates.py output/chimera/model.csv \
  -c config.ini --backbone B3 --max-swaps 2 --top 50 -o output/chimera/ranked_B3.csv

# End-to-end smoke test (rebuilds data, trains N generations, checks outputs)
python tests/verify_chimera_token_mode.py config.ini 10
```

To run classic peptide-style POET on this branch, point it at the preserved
default: `python poet.py -config configs/default_peptide.ini`.

---

## 8. Modelling assumptions and limitations

- **Additive, epistasis-free by construction.** A construct's score is the sum of
  its matched rule weights. With a `maximum_rule_size` of 4, motifs *can* span up
  to 4 adjacent slots and thus capture some local interaction, but the model is
  fundamentally additive over swaps. This is the right shape for "combine the
  measured single swaps" and will not capture strong long-range epistasis.

- **Only native-position swaps are in the training data.** Every measured
  construct swaps slot *k* with region *k* (tokens where `donor_region ==
  target_slot`). The ranker, however, enumerates **all** any-to-any tokens,
  including cross-position transplants like `S03:B2_07` (region 7 placed in slot
  3). No rule can ever match those tokens because they never occur in training,
  so such tokens contribute nothing and every candidate that differs only in an
  unrecognised cross-position token receives an identical score. In practice the
  model is a reliable predictor for **combinations of the measured single-region
  swaps**, and an uninformed extrapolation for genuine cross-position moves. If
  only same-position multi-swaps are of interest, restrict the candidate set to
  `donor_region == target_slot`.

- **Small dataset.** ~27 constructs drive a model of up to 80 weighted rules.
  Correlation fitness is measured on the training set itself (no held-out split
  in correlation mode), so a high in-sample r² is a ranking heuristic, not a
  validated generalisation estimate. Keep `parsimony_pressure` on and treat the
  ranked list as a prioritisation aid.

- **Per-sheet normalisation.** B2- and B3-backbone constructs are each scaled to
  their own sheet's B2 baseline (raw baselines differ substantially between
  sheets). Fitness values are therefore *relative* within a backbone, which is
  appropriate for ranking but means absolute cross-backbone comparisons should be
  made with care.

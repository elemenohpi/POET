import math
import os
import re
from typing import TYPE_CHECKING, Literal, Optional, cast, overload
from concurrent.futures import ProcessPoolExecutor
from scipy import stats
import numpy as np
import pandas as pd

import pattern_engine as PE
import warnings

warnings.filterwarnings("ignore")

# ── Multiprocessing worker helpers ──────────────────────────────────────────
if TYPE_CHECKING:
    _worker_fitness: Optional["Fitness"] = None
else:
    _worker_fitness = None


def _init_worker(worker_data):
    """Called once per worker process to set up a lightweight Fitness instance."""
    global _worker_fitness
    _worker_fitness = Fitness.__new__(Fitness)
    _worker_fitness.sequences = worker_data["sequences"]
    _worker_fitness.fitness_values = worker_data["fitness_values"]
    _worker_fitness.dataset_size = worker_data["dataset_size"]
    _worker_fitness.config = worker_data["config"]
    _worker_fitness.matching_mode = worker_data["matching_mode"]
    _worker_fitness.enable_gaps = worker_data.get("enable_gaps", True)
    _worker_fitness.parsimony_pressure = worker_data.get("parsimony_pressure", 0.0)
    # Experimental feature flags
    _worker_fitness.exp_char_classes = worker_data.get("exp_char_classes", False)
    _worker_fitness.exp_variable_gaps = worker_data.get("exp_variable_gaps", False)
    _worker_fitness.exp_weighted_positions = worker_data.get(
        "exp_weighted_positions", False
    )
    _worker_fitness.exp_pw_threshold = worker_data.get("exp_pw_threshold", 0.7)
    _worker_fitness.exp_composition = worker_data.get("exp_composition", False)
    _worker_fitness.exp_match_count = worker_data.get("exp_match_count", False)
    _worker_fitness.exp_circular = worker_data.get("exp_circular", False)
    _worker_fitness._use_experimental = worker_data.get("_use_experimental", False)
    _worker_fitness.k = 0


def _worker_measure(individual):
    """Evaluate a single individual in a worker process."""
    assert _worker_fitness is not None
    fitness_val, test_val = _worker_fitness.measureTotal(individual)
    rule_statuses = [(rule.status, rule.match_direction) for rule in individual.rules]
    return fitness_val, test_val, individual.usedRulesCount, rule_statuses


class Fitness:
    def __init__(self, config):
        self.config = config
        try:
            self.learn = pd.read_csv(config["learn_data"])
        except AttributeError:
            print("WARNING: settings.learn_df is not set")

        self.sequences = self.learn.iloc[:, 0].tolist()
        self.fitness_values = self.learn.iloc[:, 1].tolist()
        self.dataset_size = len(self.sequences)

        self.matching_mode = config.get("matching_mode", "substring")
        self.enable_gaps = config.get("enable_gaps", "True").strip().lower() == "true"
        self.k = 0

        # ── Experimental feature flags ──
        _bool = (
            lambda key, default="False": config.get(key, default).strip().lower()
            == "true"
        )
        self.exp_char_classes = _bool("exp_char_classes")
        self.exp_variable_gaps = _bool("exp_variable_gaps")
        self.exp_weighted_positions = _bool("exp_weighted_positions")
        self.exp_pw_threshold = float(config.get("exp_pw_threshold", "0.7"))
        self.exp_composition = _bool("exp_composition")
        self.exp_match_count = _bool("exp_match_count")
        self.exp_circular = _bool("exp_circular")
        self._use_experimental = any(
            [
                self.exp_char_classes,
                self.exp_variable_gaps,
                self.exp_weighted_positions,
                self.exp_composition,
                self.exp_match_count,
                self.exp_circular,
            ]
        )

        self.parsimony_pressure = float(config.get("parsimony_pressure", "0.0"))

        self._executor = None
        self._num_workers = 1

    def model_vs_dataset(self, config, individuals):
        seq_fitness_tuples = list(zip(self.sequences, self.fitness_values))
        individuals_evaluations = []
        values = self.fitness_values

        for individual in individuals:
            self.resetIndividual(individual)
            predictions = []
            responses = []
            for seq_fit_tuple in seq_fitness_tuples:
                error, prediction = self.eval(
                    seq_fit_tuple[0], seq_fit_tuple[1], individual, True
                )
                predictions.append(prediction)
            if config["fitness_alg"] == "correlation":
                align = np.polyfit(predictions, values, 1)
                for prediction in predictions:
                    response = prediction * align[0] + align[1]
                    responses.append(response)
                individuals_evaluations.append(responses)
            elif config["fitness_alg"] == "RMSE":
                individuals_evaluations.append(predictions)
        return seq_fitness_tuples, individuals_evaluations

    @overload
    def eval(
        self,
        sequence: str,
        actualFitness: float,
        individual,
        returnPrediction: Literal[False] = ...,
    ) -> float: ...

    @overload
    def eval(
        self,
        sequence: str,
        actualFitness: float,
        individual,
        returnPrediction: Literal[True],
    ) -> tuple[float, float]: ...

    def eval(self, sequence, actualFitness, individual, returnPrediction=False):
        if self.matching_mode == "regex":
            return self._eval_regex(
                sequence, actualFitness, individual, returnPrediction
            )
        if self._use_experimental:
            return self._eval_substring_experimental(
                sequence, actualFitness, individual, returnPrediction
            )
        return self._eval_substring(
            sequence, actualFitness, individual, returnPrediction
        )

    @staticmethod
    def _match_gap(pattern, substr):
        """Match *pattern* against *substr* treating '_' as a single-char wildcard.

        Both strings must already be the same length.  Returns True when every
        non-wildcard character in *pattern* equals the corresponding character
        in *substr*.
        """
        for pc, sc in zip(pattern, substr):
            if pc != "_" and pc != sc:
                return False
        return True

    def _eval_substring(
        self, sequence, actualFitness, individual, returnPrediction=False
    ):
        seq_len = len(sequence)
        measuredFitness = 0.0
        gaps_enabled = self.enable_gaps

        # Precompute rule info to avoid repeated work in the inner loop
        rule_data = []
        for rule in individual.rules:
            p = rule.pattern
            if not isinstance(p, str) or len(p) == 0:
                continue
            # Skip all-wildcard patterns – they match everything
            # and provide no discriminative power.
            if all(ch == "_" for ch in p):
                continue
            has_gap = gaps_enabled and "_" in p
            rule_data.append((p, p[::-1], len(p), has_gap, rule))

        for pos in range(seq_len):
            remaining = seq_len - pos
            for pattern, rev_pattern, plen, has_gap, rule in rule_data:
                if plen > remaining:
                    continue

                substr = sequence[pos : pos + plen]

                # Use fast exact comparison when there are no gaps
                fwd_match = (
                    self._match_gap(pattern, substr) if has_gap else pattern == substr
                )
                if fwd_match:
                    if rule.status == 0:
                        rule.status = 1
                        rule.match_direction = "forward"
                        individual.usedRulesCount += 1
                    elif rule.match_direction == "reverse":
                        rule.match_direction = "both"
                    measuredFitness += rule.weight
                    break

                rev_match = (
                    self._match_gap(rev_pattern, substr)
                    if has_gap
                    else rev_pattern == substr
                )
                if rev_match:
                    if rule.status == 0:
                        rule.status = 1
                        rule.match_direction = "reverse"
                        individual.usedRulesCount += 1
                    elif rule.match_direction == "forward":
                        rule.match_direction = "both"
                    measuredFitness += rule.weight
                    break

        error = abs(measuredFitness - actualFitness)
        if returnPrediction is True:
            return error, measuredFitness
        return error

    # ── Experimental substring evaluation ──────────────────────────────────

    def _eval_substring_experimental(
        self, sequence, actualFitness, individual, returnPrediction=False
    ):
        """Element-based eval supporting all experimental features."""
        seq = sequence
        seq_len = len(sequence)

        # ── Circular matching: extend sequence to wrap around ──
        max_plen = 0
        if self.exp_circular:
            for rule in individual.rules:
                p = rule.pattern
                if not isinstance(p, str) or len(p) == 0:
                    continue
                elems = PE.parse_pattern(p)
                ml = PE.elements_max_len(elems)
                if ml > max_plen:
                    max_plen = ml
            if max_plen > 1:
                seq = sequence + sequence[: max_plen - 1]
        extended_len = len(seq)

        # ── Parse & prepare rule data ──
        rule_data = []
        for rule in individual.rules:
            p = rule.pattern
            if not isinstance(p, str) or len(p) == 0:
                continue
            elements = PE.parse_pattern(p)
            if PE.is_all_wildcard(elements):
                continue
            rev_elements = PE.reverse_elements(elements)
            min_len = PE.elements_min_len(elements)
            has_vg = PE.has_variable_length(elements)
            pw = (
                rule.position_weights
                if (
                    self.exp_weighted_positions
                    and rule.position_weights
                    and not has_vg
                    and len(rule.position_weights) == len(elements)
                )
                else None
            )
            rev_pw = list(reversed(pw)) if pw else None
            rule_data.append(
                {
                    "elements": elements,
                    "rev_elements": rev_elements,
                    "min_len": min_len,
                    "pw": pw,
                    "rev_pw": rev_pw,
                    "rule": rule,
                    "match_count": 0,
                    "quality_sum": 0.0,
                    "fwd_count": 0,
                    "rev_count": 0,
                }
            )

        # ── Phase 1: find matches ──
        if self.exp_match_count:
            # Match-count mode: check each rule independently at every position
            for rd in rule_data:
                for pos in range(seq_len):
                    remaining = extended_len - pos
                    if rd["min_len"] > remaining:
                        continue
                    q = self._try_match_exp(rd, seq, pos)
                    if q > 0:
                        rd["match_count"] += 1
                        rd["quality_sum"] += q
        else:
            # Position-exclusion mode (original behaviour, first match wins)
            for pos in range(seq_len):
                remaining = extended_len - pos
                for rd in rule_data:
                    if rd["min_len"] > remaining:
                        continue
                    q = self._try_match_exp(rd, seq, pos)
                    if q > 0:
                        rd["match_count"] += 1
                        rd["quality_sum"] += q
                        break

        # ── Phase 2: update rule status ──
        for rd in rule_data:
            rule = rd["rule"]
            if rd["match_count"] == 0:
                continue
            if rule.status == 0:
                rule.status = 1
                individual.usedRulesCount += 1
            if rd["fwd_count"] > 0 and rd["rev_count"] > 0:
                rule.match_direction = "both"
            elif rd["fwd_count"] > 0:
                if rule.match_direction == "reverse":
                    rule.match_direction = "both"
                elif rule.match_direction == "":
                    rule.match_direction = "forward"
            elif rd["rev_count"] > 0:
                if rule.match_direction == "forward":
                    rule.match_direction = "both"
                elif rule.match_direction == "":
                    rule.match_direction = "reverse"

        # ── Phase 3: compute fitness ──
        measuredFitness = 0.0

        for rd in rule_data:
            rule = rd["rule"]
            # Skip grouped rules here — handled below
            if self.exp_composition and rule.group_id != 0:
                continue
            if rd["match_count"] == 0:
                continue

            contribution = rule.weight * rd["match_count"]
            # Position-weight quality scaling
            if rd["pw"] is not None and rd["match_count"] > 0:
                avg_quality = rd["quality_sum"] / rd["match_count"]
                contribution *= avg_quality

            measuredFitness += contribution

        # ── Phase 4: composition groups ──
        if self.exp_composition:
            groups: dict[int, list[dict]] = {}
            for rd in rule_data:
                gid = rd["rule"].group_id
                if gid == 0:
                    continue
                groups.setdefault(gid, []).append(rd)

            for _gid, group in groups.items():
                op = group[0]["rule"].group_op
                matched = [rd for rd in group if rd["match_count"] > 0]

                if op == "and":
                    fires = len(matched) == len(group)
                else:  # "or"
                    fires = len(matched) > 0

                if fires:
                    for rd in matched:
                        c = rd["rule"].weight * rd["match_count"]
                        if rd["pw"] is not None and rd["match_count"] > 0:
                            c *= rd["quality_sum"] / rd["match_count"]
                        measuredFitness += c

        error = abs(measuredFitness - actualFitness)
        if returnPrediction is True:
            return error, measuredFitness
        return error

    def _try_match_exp(self, rd, seq, pos):
        """Try matching at *pos*.  Returns quality > 0 on match, else 0."""
        pw = rd["pw"]
        # Forward
        if pw is not None:
            q = PE.match_at_weighted(
                rd["elements"], seq, pos, pw, self.exp_pw_threshold
            )
        else:
            q = 1.0 if PE.match_at(rd["elements"], seq, pos) else 0.0
        if q > 0:
            rd["fwd_count"] += 1
            return q

        # Reverse
        rev_pw = rd["rev_pw"]
        if rev_pw is not None:
            q = PE.match_at_weighted(
                rd["rev_elements"], seq, pos, rev_pw, self.exp_pw_threshold
            )
        else:
            q = 1.0 if PE.match_at(rd["rev_elements"], seq, pos) else 0.0
        if q > 0:
            rd["rev_count"] += 1
            return q

        return 0.0

    def _eval_regex(self, sequence, actualFitness, individual, returnPrediction=False):
        """Evaluate an individual against a sequence using regex matching."""
        measuredFitness = 0.0

        for rule in individual.rules:
            p = rule.pattern
            if not isinstance(p, str) or len(p) == 0:
                continue

            try:
                compiled = re.compile(p)
            except re.error:
                continue

            match_count = 0
            for match in re.finditer(compiled, sequence):
                match_count += 1

            if match_count > 0:
                if rule.status == 0:
                    rule.status = 1
                    rule.match_direction = "regex"
                    individual.usedRulesCount += 1
                measuredFitness += rule.weight

        error = abs(measuredFitness - actualFitness)
        if returnPrediction is True:
            return error, measuredFitness
        return error

    def resetIndividual(self, individual):
        individual.usedRulesCount = 0
        for rule in individual.rules:
            rule.status = 0
            rule.match_direction = ""

    def measure_dataset(self, individual):
        train_error = 0.0
        for seq, actual in zip(self.sequences, self.fitness_values):
            error = self.eval(seq, actual, individual)
            train_error += error**2

        MSE_train = train_error / self.dataset_size
        RMSE_train = math.sqrt(MSE_train)

        return RMSE_train, 0

    def measureTotal(self, individual):
        self.resetIndividual(individual)
        if self.config["fitness_alg"] == "RMSE":
            chunks: list[float] = [0.0] * 10
            chunk_divisor = self.dataset_size / 10
            for idx, (seq, actual) in enumerate(
                zip(self.sequences, self.fitness_values)
            ):
                j = int(idx / chunk_divisor)
                if j == 10:
                    j = 9
                error = self.eval(seq, actual, individual)
                chunks[j] += error**2

            total = sum(chunks)
            testSize = int(self.dataset_size / 10)
            trainSize = self.dataset_size - testSize

            train_rmse_sum = 0.0
            test_rmse_sum = 0.0
            for i in range(10):
                test_rmse_sum += math.sqrt(chunks[i] / testSize)
                train_rmse_sum += math.sqrt((total - chunks[i]) / trainSize)

            RMSE_train = train_rmse_sum / 10
            RMSE_test = test_rmse_sum / 10
            # Parsimony pressure: penalise model complexity
            if self.parsimony_pressure > 0:
                RMSE_train += self.parsimony_pressure * len(individual.rules)
            return RMSE_train, RMSE_test
        elif self.config["fitness_alg"] == "correlation":
            predictions = []
            for seq, actual in zip(self.sequences, self.fitness_values):
                _, prediction = self.eval(seq, actual, individual, True)
                predictions.append(prediction)
            pearsonr_result = cast(
                float, stats.pearsonr(predictions, self.fitness_values)[0]
            )
            if math.isnan(pearsonr_result):
                return 1, 0
            base_fitness = 1 - pearsonr_result**2
            # Parsimony pressure: penalise model complexity
            if self.parsimony_pressure > 0:
                base_fitness += self.parsimony_pressure * len(individual.rules)
            return base_fitness, 0
        else:
            raise ValueError(
                "Invalid fitness_alg '{}': expected 'correlation' or 'RMSE'".format(
                    self.config["fitness_alg"]
                )
            )

    # ── Parallel evaluation support ───────────────────────────────────────
    def _get_worker_data(self):
        return {
            "sequences": self.sequences,
            "fitness_values": self.fitness_values,
            "dataset_size": self.dataset_size,
            "config": self.config,
            "matching_mode": self.matching_mode,
            "enable_gaps": self.enable_gaps,
            "parsimony_pressure": self.parsimony_pressure,
            "exp_char_classes": self.exp_char_classes,
            "exp_variable_gaps": self.exp_variable_gaps,
            "exp_weighted_positions": self.exp_weighted_positions,
            "exp_pw_threshold": self.exp_pw_threshold,
            "exp_composition": self.exp_composition,
            "exp_match_count": self.exp_match_count,
            "exp_circular": self.exp_circular,
            "_use_experimental": self._use_experimental,
        }

    def start_workers(self, num_workers=None):
        """Create a persistent pool of worker processes."""
        if num_workers is None or num_workers < 1:
            num_workers = os.cpu_count() or 1
        self._num_workers = num_workers
        if self._num_workers > 1:
            self._executor = ProcessPoolExecutor(
                max_workers=self._num_workers,
                initializer=_init_worker,
                initargs=(self._get_worker_data(),),
            )
            print(
                "  Parallel fitness evaluation enabled: {} workers".format(
                    self._num_workers
                )
            )
        else:
            print("  Sequential fitness evaluation (1 worker)")

    def stop_workers(self):
        """Shut down the worker pool."""
        if self._executor:
            self._executor.shutdown(wait=True)
            self._executor = None

    def measureBatch(self, individuals):
        """Evaluate all individuals, using parallel workers if available."""
        if self._executor:
            chunksize = max(1, len(individuals) // (self._num_workers * 4))
            results = list(
                self._executor.map(_worker_measure, individuals, chunksize=chunksize)
            )
            for ind, (f, t, urc, statuses) in zip(individuals, results):
                ind.fitness = f
                ind.test = t
                ind.usedRulesCount = urc
                for rule, (status, direction) in zip(ind.rules, statuses):
                    rule.status = status
                    rule.match_direction = direction
        else:
            for ind in individuals:
                ind.fitness, ind.test = self.measureTotal(ind)

    # ── Gradient-based weight optimisation (Lamarckian step) ─────────────

    def optimize_weights(self, individual, w_min=0.0, w_max=10.0):
        """Set rule weights via bounded least-squares given current patterns.

        Builds a match matrix X[n_sequences, n_rules] where X[i,j] is the
        number of positions where rule j matches in sequence i, then solves
            min ||X w - y||²   s.t.  w_min ≤ w_j ≤ w_max
        using *scipy.optimize.lsq_linear*.
        """
        from scipy.optimize import lsq_linear

        rules = individual.rules
        k = len(rules)
        if k == 0:
            return

        n = self.dataset_size
        X = np.zeros((n, k), dtype=np.float64)

        for j, rule in enumerate(rules):
            p = rule.pattern
            if not isinstance(p, str) or len(p) == 0:
                continue
            for i, seq in enumerate(self.sequences):
                X[i, j] = self._count_matches(seq, rule)

        y = np.array(self.fitness_values, dtype=np.float64)

        # Only optimise columns that have at least one match
        active = X.any(axis=0)
        if not active.any():
            return

        X_active = X[:, active]
        lb = np.full(X_active.shape[1], w_min)
        ub = np.full(X_active.shape[1], w_max)
        result = lsq_linear(X_active, y, bounds=(lb, ub))

        idx = 0
        for j in range(k):
            if active[j]:
                rules[j].weight = round(float(result.x[idx]), 4)
                idx += 1

    def _count_matches(self, sequence, rule):
        """Count match positions for *rule* in *sequence* (independent matching)."""
        p = rule.pattern
        if not isinstance(p, str) or len(p) == 0:
            return 0
        if self.matching_mode == "regex":
            return self._count_matches_regex(sequence, p)
        if self._use_experimental:
            return self._count_matches_experimental(sequence, p)
        return self._count_matches_substring(sequence, p)

    def _count_matches_substring(self, sequence, pattern):
        seq_len = len(sequence)
        plen = len(pattern)
        if plen == 0 or plen > seq_len:
            return 0
        if all(ch == "_" for ch in pattern):
            return 0
        has_gap = self.enable_gaps and "_" in pattern
        rev = pattern[::-1]
        count = 0
        for pos in range(seq_len - plen + 1):
            sub = sequence[pos : pos + plen]
            if (self._match_gap(pattern, sub) if has_gap else pattern == sub):
                count += 1
            elif (self._match_gap(rev, sub) if has_gap else rev == sub):
                count += 1
        return count

    def _count_matches_experimental(self, sequence, pattern):
        elements = PE.parse_pattern(pattern)
        if PE.is_all_wildcard(elements):
            return 0
        rev_elements = PE.reverse_elements(elements)
        min_len = PE.elements_min_len(elements)
        seq_len = len(sequence)
        count = 0
        for pos in range(seq_len):
            remaining = seq_len - pos
            if min_len > remaining:
                break
            if PE.match_at(elements, sequence, pos):
                count += 1
            elif PE.match_at(rev_elements, sequence, pos):
                count += 1
        return count

    def _count_matches_regex(self, sequence, pattern):
        try:
            return len(re.findall(pattern, sequence))
        except re.error:
            return 0

    def predict(self, sequence, individual):
        raise NotImplementedError("Use eval() instead of predict()")

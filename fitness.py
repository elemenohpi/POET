import math
import os
import re
from typing import TYPE_CHECKING, Literal, Optional, cast, overload
from concurrent.futures import ProcessPoolExecutor
from scipy import stats
import numpy as np
import pandas as pd

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
    _worker_fitness.mode = worker_data["mode"]
    _worker_fitness.config = worker_data["config"]
    _worker_fitness.matching_mode = worker_data["matching_mode"]
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

        self.mode = int(config["pattern_mode"])
        self.matching_mode = config.get("matching_mode", "substring")
        self.k = 0

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
        return self._eval_substring(
            sequence, actualFitness, individual, returnPrediction
        )

    def _eval_substring(
        self, sequence, actualFitness, individual, returnPrediction=False
    ):
        seq_len = len(sequence)
        measuredFitness = 0.0
        mode = self.mode

        # Precompute rule info to avoid repeated work in the inner loop
        rule_data = []
        for rule in individual.rules:
            p = rule.pattern
            if not isinstance(p, str) or len(p) == 0:
                continue
            rule_data.append((p, p[::-1], len(p), rule))

        for pos in range(seq_len):
            remaining = seq_len - pos
            for pattern, rev_pattern, plen, rule in rule_data:
                if plen > remaining:
                    continue

                substr = sequence[pos : pos + plen]
                if pattern == substr:
                    if rule.status == 0:
                        rule.status = 1
                        rule.match_direction = "forward"
                        individual.usedRulesCount += 1
                    elif rule.match_direction == "reverse":
                        rule.match_direction = "both"

                    if mode == 0:
                        measuredFitness += rule.weight
                    elif mode == 1:
                        measuredFitness *= rule.weight
                    else:
                        raise ValueError(
                            "Invalid pattern_mode: expected 0 (summation) or 1 (multiplication)"
                        )
                    break
                elif rev_pattern == substr:
                    if rule.status == 0:
                        rule.status = 1
                        rule.match_direction = "reverse"
                        individual.usedRulesCount += 1
                    elif rule.match_direction == "forward":
                        rule.match_direction = "both"

                    if mode == 0:
                        measuredFitness += rule.weight
                    elif mode == 1:
                        measuredFitness *= rule.weight
                    else:
                        raise ValueError(
                            "Invalid pattern_mode: expected 0 (summation) or 1 (multiplication)"
                        )
                    break

        error = abs(measuredFitness - actualFitness)
        if returnPrediction is True:
            return error, measuredFitness
        return error

    def _eval_regex(self, sequence, actualFitness, individual, returnPrediction=False):
        """Evaluate an individual against a sequence using regex matching."""
        measuredFitness = 0.0
        mode = self.mode

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

                if mode == 0:
                    measuredFitness += rule.weight
                elif mode == 1:
                    measuredFitness *= rule.weight
                else:
                    raise ValueError(
                        "Invalid pattern_mode: expected 0 (summation) or 1 (multiplication)"
                    )

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
            return 1 - pearsonr_result**2, 0
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
            "mode": self.mode,
            "config": self.config,
            "matching_mode": self.matching_mode,
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

    def predict(self, sequence, individual):
        raise NotImplementedError("Use eval() instead of predict()")

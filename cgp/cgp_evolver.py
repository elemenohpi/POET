from pathlib import Path

import numpy as np
import numpy.random as random
import pickle
import os
import hashlib
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy

# from dnc.multiparent_wrapper import NeuralCrossoverWrapper
from .cgp_model import CGP
from .fitness_functions import correlation, corr_comp_fitness
from .helper import (
    _get_quartiles,
    pairwise_minkowski_distance,
    get_score,
    get_ssd,
    get_weights,
    clean_values,
    _get_semantic_alignment,
)

# Optional heavy dependencies — only needed for semantic crossover / SMAC tuning
try:
    from Bio.Align import PairwiseAligner, Seq
except ImportError:
    PairwiseAligner = None
    Seq = None

from matplotlib import pyplot as plt

from .cgp_generator import node_to_int

try:
    from ConfigSpace import (
        Categorical,
        Configuration,
        ConfigurationSpace,
        Float,
        Integer,
    )
    from ConfigSpace.conditions import InCondition
    from smac import HyperparameterOptimizationFacade, Scenario
except ImportError:
    ConfigSpace = None
    Scenario = None
    HyperparameterOptimizationFacade = None


class CartesianGP:

    def _missing_key_error(self, key):
        raise KeyError(f"'{key}' is required but was not provided.")

    def __init__(
        self,
        parents=1,
        children=4,
        max_generations=100,
        mutation="Point",
        selection="Elite",
        xover=None,
        fixed_length=True,
        fitness_function="correlation_complexity",
        model_parameters=None,
        function_bank=None,
        solution_threshold=0.005,
        checkpoint_filename="checkpoint.pkl",
        seed=42,
        dnc_hp: dict = None,
        **kwargs,
    ):
        # hyperparameter tuning
        # Basic attributes
        self.model_keys = None
        np.random.seed(seed)
        self.model_key_map = None
        self.child_id_counter = 0
        self.population = np.empty(parents + children, dtype=object)
        self.fitnesses = np.full(parents + children, np.inf, dtype=np.float64)
        self.corrs = np.full(parents + children, np.inf, dtype=np.float64)
        self.fitnesses_test = np.full(parents + children, np.inf, dtype=np.float64)
        self.corr_test = np.full(parents + children, np.inf, dtype=np.float64)
        self.best_model = None
        self.function_bank = function_bank
        self.x = None
        self.x_test = None
        self.y = None
        self.y_test = None
        self.solution_threshold = solution_threshold
        self.max_p = parents
        self.max_c = children
        self.max_g = max_generations
        self.original_max_g = max_generations
        self.fixed_length = fixed_length
        self.model_kwargs = model_parameters or {}
        self.ckpt_filename = checkpoint_filename
        self.current_generation = 0
        self.kwargs = kwargs
        self.first_submission = True
        self.dnc_hyperparams = None
        self.dnc = None
        # Config values with defaults
        self.mutation_type = mutation.lower()
        self.selection_type = selection.lower()
        self.xover_type = xover.lower() if xover else None
        self.ff_string = fitness_function
        self.semantic = "semantic" in self.xover_type if xover else False
        self.aligned = "aligned" in self.xover_type if xover else False
        self.homologous = "homologous" in self.xover_type if xover else False

        self.pareto_layers = []
        self.pareto_elite = []

        # Safe defaults
        self.tournament_size = int(kwargs.get("tournament_size", 2))
        self.n_elites = int(kwargs.get("n_elites", 0))
        self.n_points = int(kwargs.get("n_points", 1))
        self.tournament_diversity = kwargs.get("tournament_diversity", True)
        self.one_d = kwargs.get("one_dimensional_xover", False)

        # Mutation fallback
        self.mutation_can_make_children = self.max_p < 2 or kwargs.get(
            "asexual_reproduction", False
        )

        # Fitness function
        self.fitness_function = {
            "correlation": correlation,
            "correlation_complexity": corr_comp_fitness,
        }.get(fitness_function.lower())
        if not self.fitness_function:
            raise ValueError(f"Invalid fitness function: {fitness_function}")

        # Selection methods
        self.selection_methods = {
            "elite": self.elite_selection,
            "paretoelite": self.pareto_elite_selection,
            "tournament": self.tournament_selection,
            "paretotournament": self.pareto_tournament_selection,
            "elite_tournament": self.elite_tournament_selection,
            "competent_tournament": self.competent_tournament_selection,
        }
        if self.selection_type not in self.selection_methods:
            raise ValueError(f"Invalid selection type: {self.selection_type}")
        self.selection = self.selection_methods[self.selection_type]

        # Crossover setup
        self.xover_methods = {
            "none": None,
            "n_point": self._n_point_xover,
            "uniform": self._uniform_xover,
            "semantic_uniform": self._uniform_xover,
            # 'dnc_semantic_uniform': self._uniform_xover,
            # 'dnc_uniform': self._uniform_xover,
            "aligned_semantic_uniform": self._uniform_xover,
            "aligned_homologous_semantic_uniform": self._uniform_xover,
            "homologous_semantic_uniform": self._uniform_xover,
            "semantic_n_point": self._n_point_xover,
            # 'dnc_semantic_n_point': self._n_point_xover,
            # 'dnc_n_point': self._n_point_xover,
            "homologous_semantic_n_point": self._n_point_xover,
            "aligned_homologous_semantic_n_point": self._n_point_xover,
            "aligned_semantic_n_point": self._n_point_xover,
            "subgraph": self._subgraph_xover,
        }
        if self.xover_type:
            if self.xover_type not in self.xover_methods:
                raise ValueError(f"Invalid crossover type: {self.xover_type}")
            # if 'dnc' in self.xover_type:
            #    self.dnc_hyperparams = dnc_hp
            #    if 'uniform' in self.xover_type:
            #        dnc_xover = 'uniform'
            #    elif 'n_point' in self.xover_type:
            #        dnc_xover = 'n_point'
            #    else:
            #        raise KeyError(f"{self.xover_type} is not valid for use with Deep Neural Crossover.")
            #    self.dnc = NeuralCrossoverWrapper(**self.dnc_hyperparams, crossover_type=dnc_xover,
            #                                      n_points=self.n_points)
            self.xover = self.xover_methods[self.xover_type]
        else:
            self.xover = None

        # one dimensional xover is incompatible with semantic and subgraph methods:
        print(f"1d crossover: {self.one_d}")
        if self.one_d and (
            "semantic" in self.xover_type or "subgraph" in self.xover_type
        ):
            raise ValueError(
                f"{self.xover_type} crossover is incompatible with one-dimensional crossover."
            )

        # Validate crossover points
        if self.xover_type and "n_point" in self.xover_type:
            max_size = self.model_kwargs.get("max_size", 10)
            if self.n_points < 1 or self.n_points > max_size // 2:
                raise ValueError(f"Invalid n_points: {self.n_points}")

        # Metrics and tracking
        self.metrics = np.zeros((self.max_g + 1, 32), dtype=np.float64)
        self.xover_index = {
            cat: np.zeros((self.max_g, self.max_p))
            for cat in ["deleterious", "neutral", "beneficial"]
        }
        self.mut_index = np.zeros((self.max_g, self.max_p))

    @staticmethod
    def hash_model(m):
        return hashlib.md5(m.tobytes()).hexdigest()

    def save_checkpoint(self, filename="cgp_checkpoint.pkl", generation=None):
        temp_file = filename[0:-4] + ".tmp"
        with open(temp_file, "wb") as f:
            pickle.dump(self.__dict__, f)
        os.replace(temp_file, filename)  # atomic rename to avoid partial writes
        print(
            f"Checkpoint saved at generation {generation or self.current_generation} in {filename}"
        )

    @classmethod
    def load_checkpoint(cls, filename="cgp_checkpoint.pkl"):
        def _try_load(path):
            with open(path, "rb") as f:
                return pickle.load(f)

        try:
            data = _try_load(filename)
        except EOFError:
            print(
                f"ΓÜá∩╕Å Warning: Failed to load checkpoint '{filename}' (EOFError). Trying backup..."
            )
            tmp_file = filename[0:-4] + ".tmp"
            if os.path.exists(tmp_file):
                try:
                    data = _try_load(tmp_file)
                    print("Γ£à Loaded from backup:", tmp_file)
                except Exception as e:
                    raise RuntimeError(
                        f"Γ¥î Failed to load from both '{filename}' and backup: {e}"
                    )
            else:
                raise RuntimeError(
                    f"Γ¥î Checkpoint corrupted and no backup found: {filename}"
                )

        obj = cls.__new__(cls)  # Don't call __init__!
        obj.__dict__.update(data)

        # Restore method bindings after unpickling
        obj.selection_methods = {
            "elite": obj.elite_selection,
            "paretoelite": obj.pareto_elite_selection,
            "tournament": obj.tournament_selection,
            "paretotournament": obj.pareto_tournament_selection,
            "elite_tournament": obj.elite_tournament_selection,
            "competent_tournament": obj.competent_tournament_selection,
        }

        obj.xover_methods = {
            "none": None,
            "n_point": obj._n_point_xover,
            "uniform": obj._uniform_xover,
            "semantic_uniform": obj._uniform_xover,
            "aligned_semantic_uniform": obj._uniform_xover,
            "aligned_homologous_semantic_uniform": obj._uniform_xover,
            "homologous_semantic_uniform": obj._uniform_xover,
            "semantic_n_point": obj._n_point_xover,
            "homologous_semantic_n_point": obj._n_point_xover,
            "aligned_homologous_semantic_n_point": obj._n_point_xover,
            "aligned_semantic_n_point": obj._n_point_xover,
            "subgraph": obj._subgraph_xover,
        }

        obj.selection = obj.selection_methods.get(obj.selection_type)
        obj.xover = obj.xover_methods.get(obj.xover_type)

        obj.first_submission = False

        print(f"Checkpoint loaded at generation {obj.current_generation}")
        return obj

    def initialize_population(self):
        """Initialize population in parallel."""
        with ThreadPoolExecutor() as executor:
            future_models = {
                i: executor.submit(
                    CGP,
                    fixed_length=self.fixed_length,
                    fitness_function=self.ff_string,
                    mutation_type=self.mutation_type,
                    function_bank=self.function_bank,
                    **self.model_kwargs,
                )
                for i in range(self.max_p)
            }

        for i, future in future_models.items():
            ind = future.result()

            # Γ£à Force re-initialize xover_index for consistency in 1D
            if self.one_d:
                n_zeros = (ind.arity + 1) * ind.max_size + ind.outputs
                ind.xover_index = np.zeros(n_zeros)
                print(
                    f"≡ƒôª Initialized xover_index for individual {i}: {ind.xover_index.shape}"
                )
            else:
                ind.xover_index = np.zeros(ind.max_size + ind.outputs)

            self.population[i] = ind
        self.model_keys = deepcopy(self.population[0].model_keys)

    def elite_selection(self, models=None, n_elites=None, indices=False):
        """Select top-n elite individuals based on their fitness attribute."""
        n_elites = n_elites or self.max_p
        # print(self.population)
        if models is None:
            models = self.population
        # Combine models and their indices
        valid_pairs = [(i, m) for i, m in enumerate(self.population) if m is not None]
        np.random.shuffle(valid_pairs)

        # Sort pairs by model fitness
        sorted_pairs = sorted(valid_pairs, key=lambda x: x[1].fitness)
        # Unpack sorted indices and models
        sorted_indices, sorted_models = zip(*sorted_pairs) if sorted_pairs else ([], [])
        # print("Fitnesses:", self.fitnesses)
        # print("Elite indices:", elite_indices)
        if indices:
            return [deepcopy(m) for m in sorted_models[:n_elites]], sorted_indices[
                :n_elites
            ]
        else:
            return [deepcopy(m) for m in sorted_models[:n_elites]]

    """
    def get_pareto_front(self, points):
        is_efficient = np.arange(points.shape[0])
        n_points = points.shape[0]
        next_point_index = 0
        while next_point_index < len(points):
            nondominated_point_mask = np.any(points < points[next_point_index], axis=1)
            nondominated_point_mask[next_point_index] = True
            is_efficient = is_efficient[nondominated_point_mask]
            points = points[nondominated_point_mask]
            next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
        return is_efficient
    """

    def get_pareto_front(self, points: np.ndarray) -> np.ndarray:
        """
        points: shape (n, 2) with columns [correlation, complexity]
        returns: array of indices of Pareto-optimal points (minimization)
        """
        if points.size == 0:
            return np.array([], dtype=int)

        # sort by first objective (correlation)
        order = np.argsort(points[:, 0])
        best_complexity = np.inf
        front_indices = []

        for idx in order:
            comp = points[idx, 1]
            if comp < best_complexity:
                best_complexity = comp
                front_indices.append(idx)

        return np.array(front_indices, dtype=int)

    def pareto_elite_selection(self, models=None, n_elites=None, return_indices=False):
        """
        Select elites based on 2D objectives (correlation, complexity).

        - correlation dominates (primary objective, minimize)
        - complexity is secondary (tie-break / small pressure)
        - always returns exactly n_elites models (or indices if return_indices=True)
        """
        if models is None:
            models = self.population

        # remove Nones but remember their indices
        valid_pairs = [(i, m) for i, m in enumerate(models) if m is not None]
        if not valid_pairs:
            return [] if not return_indices else np.array([], dtype=int)

        orig_indices, models_valid = zip(*valid_pairs)
        orig_indices = np.array(orig_indices)
        models_valid = list(models_valid)

        n_elites = n_elites or self.n_elites

        # build points array: [corr, complexity]
        points = np.array(
            [[m.correlation, m.complexity] for m in models_valid], dtype=float
        )

        # 1) find Pareto front indices (relative to models_valid)
        front_idx = self.get_pareto_front(points)

        # 2) rank points on the front lexicographically: corr first, then complexity
        front_points = points[front_idx]

        corr = front_points[:, 0]
        comp = front_points[:, 1]

        eps_corr = 1e-3  # ΓåÉ tune this

        best_corr = np.min(corr)
        corr_bucket = np.floor((corr - best_corr) / eps_corr)

        order_front = np.lexsort((comp, corr_bucket))
        front_sorted_idx = front_idx[order_front]

        # 3) pick as many elites as possible from the front
        chosen_idx = list(front_sorted_idx[:n_elites])

        # 4) if still need more elites, fill from remaining models by lexicographic order
        if len(chosen_idx) < n_elites:
            remaining_mask = np.ones(len(models_valid), dtype=bool)
            remaining_mask[chosen_idx] = False
            remaining_points = points[remaining_mask]

            order_all = np.lexsort((remaining_points[:, 1], remaining_points[:, 0]))
            remaining_indices = np.where(remaining_mask)[0][order_all]

            needed = n_elites - len(chosen_idx)
            chosen_idx.extend(list(remaining_indices[:needed]))

        chosen_idx = np.array(chosen_idx, dtype=int)

        if return_indices:
            # map back to original indices of self.population
            return orig_indices[chosen_idx]

        # otherwise return deep-copied models
        elites = [deepcopy(models_valid[i]) for i in chosen_idx]
        return elites

    def tournament_selection(self, n_to_select=None):
        """Performs tournament selection with optional diversity enforcement."""
        n_to_select = n_to_select or self.max_p
        new_population = np.empty(n_to_select, dtype=object)

        # Filter out None values from population
        t_pop = self.population[self.population is not None][0]

        available_indices = list(range(len(t_pop)))
        remaining_slots = len(new_population)

        new_population = self.t_select(
            available_indices, new_population, remaining_slots, pareto=False
        )

        return new_population

    def pareto_tournament_selection(self, n_to_select=None):
        n_to_select = n_to_select or self.max_p
        new_population = np.empty(n_to_select, dtype=object)

        # Filter out None values from population
        available_indices = np.array(
            [i for i, m in enumerate(self.population) if m is not None]
        )

        # available_indices = list(range(len(t_pop)))
        remaining_slots = len(new_population)

        new_population = self.t_select(
            available_indices, new_population, remaining_slots, pareto=True
        )

        return new_population

    def elite_tournament_selection(self):
        """Combines elite selection with tournament selection, ensuring diversity enforcement."""

        # Step 1: Select elite individuals (from self.population)
        elite_population, elite_indices = self.elite_selection(
            n_elites=self.n_elites, indices=True
        )
        # Step 2: Prepare for tournament selection
        remaining_slots = self.max_p - len(elite_population)
        # Create available index pool from self.population
        available_indices = [
            i
            for i, ind in enumerate(self.population)
            if i not in elite_indices and ind is not None
        ]
        # Step 3: Perform tournament selection
        new_population = np.empty(remaining_slots, dtype=object)
        available_indices, new_population = self.t_select(
            available_indices, new_population, remaining_slots
        )

        # Step 4: Combine elite and tournament-selected individuals
        final_population = np.concatenate((elite_population, new_population))
        return final_population

    def t_select(
        self, available_indices, new_population, remaining_slots, pareto=False
    ):
        for i in range(remaining_slots):
            if (
                self.tournament_diversity
                and len(available_indices) < self.tournament_size
            ):
                contestants_indices = available_indices
            else:
                contestants_indices = np.random.choice(
                    available_indices, size=self.tournament_size, replace=False
                )

            # Select best contestant
            if not pareto:
                best_index = min(
                    contestants_indices, key=lambda idx: self.population[idx].fitness
                )
            else:
                selected_models = [self.population[i] for i in contestants_indices]

                best_contestant = self.pareto_elite_selection(
                    models=selected_models, return_indices=True, n_elites=1
                )
                best_index = contestants_indices[best_contestant[0]]
            new_population[i] = deepcopy(self.population[best_index])

            # Remove selected individual to enforce diversity
            if self.tournament_diversity:
                available_indices = available_indices[available_indices != best_index]

        return new_population

    def _compute_semantics(self):
        """
        Compute semantics for all individuals in the CGP population at once.
        Uses vectorized operations for better efficiency.

        Returns:
            np.ndarray: A matrix where each row is an individual's output.
        """
        valid_individuals = [ind for ind in self.population if ind is not None]

        if not valid_individuals:  # Handle empty population case
            return np.empty((0, len(self.x)))

        # Compute outputs for all individuals in one go
        semantics = np.vstack([ind(self.x).flatten() for ind in valid_individuals])

        return semantics

    def competent_tournament_selection(self, n_to_select=None):
        """
        Perform competent tournament selection with semantic distance-based scoring,
        using optimized vectorized operations.

        Returns:
            np.ndarray: Array of selected parent CGP models.
        """
        n_to_select = n_to_select or self.max_p

        # Filter out None individuals and get their indices
        t_pop = np.array([ind for ind in self.population if ind is not None])
        parent_indices = np.arange(len(t_pop))

        # Compute all parent semantics at once
        parent_semantics = self._compute_semantics()

        # Compute target vector (ground truth)
        def _normalize(v):
            return (v - np.mean(v)) / (np.std(v) + 1e-8)

        parent_semantics = np.apply_along_axis(_normalize, 1, parent_semantics)
        target = _normalize(np.ravel(self.y))
        # Compute distances to target using vectorized Minkowski distance (p=2)
        target_distances = np.linalg.norm(parent_semantics - target, axis=1)

        # Prepare selected indices and remaining indices
        selected_indices = []
        remaining = list(parent_indices)

        t_size = min(self.tournament_size, len(remaining))
        enforce_diversity = self.tournament_diversity

        while len(selected_indices) < n_to_select and len(remaining) >= t_size:

            # Ensure we have enough candidates to sample
            if len(remaining) < t_size:
                break

            # Sample contestants
            contestants = np.random.choice(remaining, size=t_size, replace=False)

            # Select first parent (lowest target distance)
            first_index = contestants[np.argmin(target_distances[contestants])]

            first_sem = parent_semantics[first_index]
            first_target_dist = target_distances[first_index]

            # Compute semantic distances between first parent and all contestants
            sem_distances = np.linalg.norm(
                parent_semantics[contestants] - first_sem, axis=1
            )

            # Compute composite scores
            scores = {
                idx: get_score(
                    first_target_dist,
                    pairwise_minkowski_distance(first_sem, parent_semantics[idx], p=2),
                    target_distances[idx],
                )
                for idx in contestants
            }

            # Select second parent (minimum composite score)
            second_index = min(scores, key=scores.get)

            # Add both parents to selected list
            selected_indices.extend([first_index, second_index])

            # Enforce diversity
            if enforce_diversity:
                remaining = [
                    i for i in remaining if i not in (first_index, second_index)
                ]

        return self.population[np.array(selected_indices[:n_to_select])]

    def crossover(self, parents, xover_rate, gen):
        """Perform crossover to generate children."""
        children = []
        parent_pairs = [(parents[i], parents[i + 1]) for i in range(0, len(parents), 2)]

        while len(children) < self.max_c:
            for i, (p1, p2) in enumerate(parent_pairs):
                if np.random.rand() < xover_rate:
                    if self.xover_type == "subgraph":
                        c1 = self.xover(p1, p2, gen)
                        c2 = self.xover(p2, p1, gen)
                    # elif 'dnc' in self.xover_type:
                    #    parent_pairs = [(parents[i], parents[i + 1]) for i in range(0, len(parents), 2)]
                    #    semantic_pairs = [(clean_values(p[0], self.x, include_output=True),
                    #                       clean_values(p[1], self.x, include_output=True)) for p in
                    #                      parent_pairs] if self.semantic else None
                    #    offspring_pairs, _ = self.dnc.cross_pairs(parent_pairs, self.x, self.y,
                    #                                              semantic_pairs)  # returns children, updates training

                    #    # Flatten pairs into a single list of individuals
                    #    new_children = [CGP(model=c, model_keys=self.model_keys, fixed_length=self.fixed_length,
                    #                        fitness_function=self.ff_string,
                    #                        mutation_type=self.mutation_type) for pair in offspring_pairs[:2] for c in
                    #                    pair]
                    #    c1, c2 = new_children[:2]
                    else:
                        c1, c2 = self.xover(p1, p2, gen)
                else:
                    c1, c2 = deepcopy(p1), deepcopy(p2)

                # Get parent keys
                p1_key = getattr(p1, "child_keys", f"Model_{2 * i:03d}_g{gen - 1}")
                p2_key = getattr(p2, "child_keys", f"Model_{2 * i + 1:03d}_g{gen - 1}")

                # Assign keys to both children
                for child in [c1, c2]:
                    child_key = f"Child_{self.child_id_counter:03d}_g{gen}"
                    child.set_parent_key([p1_key, p2_key])
                    child.set_child_key(child_key)

                    if hasattr(self, "model_key_map"):
                        self.model_key_map[child_key] = child

                    self.child_id_counter += 1
                    children.append(child)

                    if len(children) >= self.max_c:
                        break
                if len(children) >= self.max_c:
                    break
        return np.array(children, dtype=object)

    def flatten_parent(self, parent):
        """
        Flattens a CGP individual's function/output nodes into a 1D array.

        Returns:
            flat: np.ndarray of [Operator, Operand0, ..., Operand{arity-1}, ...]
            node_types: np.ndarray of corresponding NodeTypes (to distinguish Function vs Output)
        """
        mask = parent.model[:, self.model_keys["NodeType"]] == node_to_int("Function")
        nodes = parent.model[mask]
        arity = parent.arity  # number of operands per function/output node
        flat = []
        for node in nodes:
            flat.append(node[self.model_keys["Operator"]])
            for j in range(arity):
                flat.append(
                    int(node[self.model_keys[f"Operand{j}"]])
                )  # ensure operands are stored as int
        return np.array(flat, dtype=np.int64), nodes[:, self.model_keys["NodeType"]]

    def unflatten_model(self, flattened_model, node_types, arity):
        """
        Reconstructs a structured NumPy model array from a flat vector and node types.

        Args:
            flattened_model (np.ndarray): flat vector of alternating [Operator, Operand0...n]
            node_types (np.ndarray): array of NodeType strings ('Function' or 'Output')
            arity (int): number of operands per node

        Returns:
            np.ndarray: structured model array of nodes (Function + Output only)
        """
        num_nodes = len(node_types)
        stride = 1 + arity  # number of fields per node
        model = np.zeros((num_nodes, len(self.model_keys)), dtype=np.int64)

        for i in range(num_nodes):
            offset = i * stride
            model[i, self.model_keys["NodeType"]] = node_types[i]
            model[i, self.model_keys["Operator"]] = int(flattened_model[offset])
            for j in range(arity):
                operand_value = flattened_model[offset + j]
                model[i][self.model_keys[f"Operand{j}"]] = int(operand_value)

            model[i, self.model_keys["Value"]] = 0.0
            model[i, self.model_keys["Active"]] = 0.0
        return model

    def _n_point_xover(self, p1, p2, gen, **kwargs):
        """Performs n-point crossover, supporting semantic and homologous crossover."""

        def get_crossover_points(length, offset=0, weights=None):
            indices = np.arange(offset, length)
            if weights is not None and weights.sum() > 0:
                if len(indices) > len(weights):
                    indices = indices[: len(weights)]
                return np.sort(
                    np.random.choice(
                        indices, size=self.n_points, replace=False, p=weights
                    )
                )
            return np.sort(np.random.choice(indices, size=self.n_points, replace=False))

        if self.one_d:
            flat_p1, types_p1 = self.flatten_parent(p1)
            flat_p2, types_p2 = self.flatten_parent(p2)
            xover_length = len(flat_p1)

            len1 = len(types_p1)

            xover_points_p1 = get_crossover_points(len1)

            # Optional: track stats here
            for x in xover_points_p1:
                p1.xover_index[x] += 1
            for x in xover_points_p1:
                p2.xover_index[x] += 1

            # Split and recombine
            stride = 1 + p1.arity
            xover_points_p1 = xover_points_p1[xover_points_p1 < len(types_p1)]

            parts1 = np.split(flat_p1, xover_points_p1)
            parts2 = np.split(flat_p2, xover_points_p1)
            types_parts1 = np.split(types_p1, xover_points_p1)
            types_parts2 = np.split(types_p2, xover_points_p1)

            child1_flat = np.concatenate(parts1[::2] + parts2[1::2])
            child2_flat = np.concatenate(parts2[::2] + parts1[1::2])
            child1_types = np.concatenate(types_parts1[::2] + types_parts2[1::2])
            child2_types = np.concatenate(types_parts2[::2] + types_parts1[1::2])

            # Γ£à Make sure unflatten won't crash
            assert (
                len(child1_flat) == len(child1_types) * stride
            ), f"child1_flat={len(child1_flat)} vs types={len(child1_types)}"
            assert (
                len(child2_flat) == len(child2_types) * stride
            ), f"child2_flat={len(child2_flat)} vs types={len(child2_types)}"

            # Reconstruct models
            child1_model = self.unflatten_model(
                child1_flat, child1_types, arity=p1.arity
            )
            child2_model = self.unflatten_model(
                child2_flat, child2_types, arity=p2.arity
            )
            o1_nodes = [
                node
                for node in p1.model
                if node[self.model_keys["NodeType"]] == node_to_int("Output")
            ]
            o2_nodes = [
                node
                for node in p2.model
                if node[self.model_keys["NodeType"]] == node_to_int("Output")
            ]

            # Add back Input and Constant nodes
            def insert_io_nodes(original, body, o_nodes):
                i_nodes = [
                    node
                    for node in original.model
                    if node[self.model_keys["NodeType"]]
                    in map(node_to_int, ["Input", "Constant"])
                ]
                return np.array(
                    i_nodes + list(body) + o_nodes, dtype=original.model.dtype
                )

            full_model_1 = insert_io_nodes(p1, child1_model, o2_nodes)
            full_model_2 = insert_io_nodes(p2, child2_model, o1_nodes)

            c1 = CGP(
                model=full_model_1,
                model_keys=self.model_keys,
                fixed_length=self.fixed_length,
                fitness_function=self.ff_string,
                mutation_type=self.mutation_type,
                xover_length=xover_length,
            )
            c2 = CGP(
                model=full_model_2,
                model_keys=self.model_keys,
                fixed_length=self.fixed_length,
                fitness_function=self.ff_string,
                mutation_type=self.mutation_type,
                xover_length=xover_length,
            )
            return c1, c2

        # Else: fixed-length 2D crossover
        fb_node = min(p1.first_body_node, p2.first_body_node)
        weights = None
        if self.semantic:
            if self.aligned:
                weights = _get_semantic_alignment(p1, p2, self.x)
            else:
                vmat_1, vmat_2 = clean_values(p1, self.x), clean_values(p2, self.x)
                weights = get_weights(get_ssd(vmat_1, vmat_2))
            if self.homologous and not np.all(weights == weights[0]):
                weights = (weights.max() - weights) / (
                    weights.max() - weights.min() + 1e-8
                )
                weights /= weights.sum()
        xover_points_p1 = get_crossover_points(
            len(p1.model), p1.first_body_node, weights
        )

        p1.xover_index[xover_points_p1 - fb_node] += 1
        p2.xover_index[xover_points_p1 - fb_node] += 1

        parts1 = np.split(p1.model, xover_points_p1)
        parts2 = np.split(p2.model, xover_points_p1)
        child1 = np.concatenate(parts1[::2] + parts2[1::2])
        child2 = np.concatenate(parts2[::2] + parts1[1::2])

        c1 = CGP(
            model=child1,
            model_keys=self.model_keys,
            fixed_length=self.fixed_length,
            fitness_function=self.ff_string,
            mutation_type=self.mutation_type,
        )
        c2 = CGP(
            model=child2,
            model_keys=self.model_keys,
            fixed_length=self.fixed_length,
            fitness_function=self.ff_string,
            mutation_type=self.mutation_type,
        )
        return c1, c2

    def _uniform_xover(self, p1, p2, gen, weights: np.ndarray | list = None, **kwargs):
        """Performs uniform crossover, supporting semantic and homologous crossover."""

        if self.one_d:
            flat1, types1 = self.flatten_parent(p1)
            flat2, types2 = self.flatten_parent(p2)

            xover_length = len(flat1)

            assert len(flat1) == len(
                flat2
            ), "Flattened parents must have the same length for uniform crossover."

            # Create swap mask
            swap_mask = np.random.rand(len(flat1)) < 0.5
            child1_flat = flat1.copy()
            child2_flat = flat2.copy()
            child1_flat[swap_mask] = flat2[swap_mask]
            child2_flat[swap_mask] = flat1[swap_mask]
            o1_nodes = [
                node
                for node in p1.model
                if node[self.model_keys["NodeType"]] == node_to_int("Output")
            ]
            o2_nodes = [
                node
                for node in p2.model
                if node[self.model_keys["NodeType"]] == node_to_int("Output")
            ]

            o_nodes = []
            o_nodes_complement = []

            o_count = -1
            for n1, n2 in zip(o1_nodes, o2_nodes):
                if np.random.rand() > 0.5:
                    o_nodes.append(n1)
                    p1.xover_index[o_count - p1.arity + 1] += 1
                    o_nodes_complement.append(n2)
                else:
                    o_nodes.append(n2)
                    p2.xover_index[o_count - p2.arity + 1] += 1
                    o_nodes_complement.append(n1)
                o_count -= 1

                # Track swapped indices
            for idx in np.where(swap_mask)[0]:
                assert np.max(np.where(swap_mask)) < len(
                    p1.xover_index
                ), "swap index out of bounds!"

                p1.xover_index[idx] += 1
                p2.xover_index[idx] += 1

            # Reconstruct models
            child1_body = self.unflatten_model(child1_flat, types1, arity=p1.arity)
            child2_body = self.unflatten_model(child2_flat, types2, arity=p2.arity)

            # Add back input + constant nodes
            def insert_io_nodes(original, body, o_nodes):
                i_nodes = [
                    node
                    for node in original.model
                    if node[self.model_keys["NodeType"]]
                    in map(node_to_int, ["Input", "Constant"])
                ]
                return np.array(
                    i_nodes + list(body) + o_nodes, dtype=original.model.dtype
                )

            full_model_1 = insert_io_nodes(p1, child1_body, o_nodes_complement)
            full_model_2 = insert_io_nodes(p2, child2_body, o_nodes)

            c1 = CGP(
                model=full_model_1,
                model_keys=self.model_keys,
                fixed_length=self.fixed_length,
                fitness_function=self.ff_string,
                mutation_type=self.mutation_type,
                xover_length=xover_length,
            )

            c2 = CGP(
                model=full_model_2,
                model_keys=self.model_keys,
                fixed_length=self.fixed_length,
                fitness_function=self.ff_string,
                mutation_type=self.mutation_type,
                xover_length=xover_length,
            )

            return c1, c2

        # Standard 2D structured uniform crossover
        n_outputs = 0
        if self.semantic:
            if self.aligned:
                weights = _get_semantic_alignment(p1, p2, self.x)
            else:
                vmat_1, vmat_2 = clean_values(p1, self.x), clean_values(p2, self.x)
                weights = get_weights(get_ssd(vmat_1, vmat_2), epsilon=0.001)

            if self.homologous and not np.all(weights == weights[0]):
                max_weight, min_weight = weights.max(), weights.min()
                weights = (max_weight - weights) / (max_weight - min_weight + 1e-8)
                weights /= weights.sum()

            n_outputs = p1.outputs

        fb_node = min(p1.first_body_node, p2.first_body_node)
        assert len(p1.model) == len(
            p2.model
        ), "Parents in Uniform Xover must have the same length."
        possible_indices = np.arange(fb_node, len(p1.model) - n_outputs)

        if weights is not None:
            mask = weights > 1e-8
            filtered_indices = possible_indices[mask]
            filtered_weights = weights[mask]
            n_swap = min(len(filtered_indices), len(possible_indices) // 2)

            if n_swap > 0:
                swapped_indices = np.random.choice(
                    filtered_indices,
                    size=n_swap,
                    replace=False,
                    p=filtered_weights / filtered_weights.sum(),
                )
            else:
                swapped_indices = np.array([], dtype=int)
        else:
            swapped_indices = np.random.choice(
                possible_indices, size=len(possible_indices) // 2, replace=False
            )

        # Perform crossover
        c1_model, c2_model = p1.model.copy(), p2.model.copy()
        c1_model[swapped_indices] = p2.model[swapped_indices]
        c2_model[swapped_indices] = p1.model[swapped_indices]

        # Track crossover stats
        p1.xover_index[swapped_indices - fb_node] += 1
        p2.xover_index[swapped_indices - fb_node] += 1

        # Final CGP children
        c1 = CGP(
            model=c1_model,
            model_keys=self.model_keys,
            fixed_length=self.fixed_length,
            fitness_function=self.ff_string,
            mutation_type=self.mutation_type,
        )
        c2 = CGP(
            model=c2_model,
            model_keys=self.model_keys,
            fixed_length=self.fixed_length,
            fitness_function=self.ff_string,
            mutation_type=self.mutation_type,
        )
        return c1, c2

    def _subgraph_xover(self, p1: CGP, p2: CGP, gen: int, **kwargs):
        """Performs subgraph crossover using NumPy arrays."""

        def random_node_number(n_i, I=None, n_f=None, m=None):
            """Selects a random node number with constraints."""
            n_r = []  # List of valid random choices

            if n_f is not None:
                if m is not None:
                    n_m = n_f[n_f <= m]
                    if len(n_m) == 0:
                        n_r.append(np.random.randint(0, n_i))
                    else:
                        n_r.append(np.random.choice(n_m))
                else:
                    n_r.append(np.random.choice(n_f))

            if I is not None:
                n_r.append(np.random.choice(I))

            return np.random.choice(n_r)

        def determine_crossover_point(m1, m2):
            """Determines a crossover point for two models."""
            a, b, c, d = min(m1), max(m1), min(m2), max(m2)
            if a >= b:
                b += 1  # Ensure valid range
            if c >= d:
                d += 1
            cp1, cp2 = np.random.randint(a, b), np.random.randint(c, d)
            return min(cp1, cp2)

        def neighborhood_connect(nf, nb, model):
            """Ensures neighborhood connectivity in the new model."""
            model[nb, self.model_keys["Operand0"]] = nf
            return model

        def random_active_connect(n_i, n_a, c_p, model):
            """Reconnects inactive nodes in the new model."""
            operand_indices = [f"Operand{i}" for i in range(p1.arity)]
            input_nodes = np.where(
                (model[self.model_keys["NodeType"]] == node_to_int("Constant"))
                | (model[self.model_keys["NodeType"]] == node_to_int("Input"))
            )[0]

            for n in n_a:
                if n > c_p:
                    for operand in [self.model_keys[op] for op in operand_indices]:
                        if model[n, operand] not in n_a:
                            model[n, operand] = random_node_number(
                                n_i, I=input_nodes, n_f=n_a, m=c_p
                            )

            output_nodes = np.where(
                model[self.model_keys["NodeType"]] == node_to_int("Output")
            )[0]
            for idx in output_nodes:
                if model[idx, self.model_keys["Operand0"]] not in n_a:
                    model[idx, self.model_keys["Operand0"]] = random_node_number(
                        n_i, I=input_nodes, n_f=n_a
                    )

            return model

        def ensure_active_nodes(parent, n_max_mutations=64):
            """Ensures a model has active nodes before crossover."""
            active_nodes = np.array(list(parent.get_active_nodes()))
            while (
                len(active_nodes) < 1
            ):  # If there are no active nodes, mutate output until one is found
                parent.mutate_output()
                generic = np.zeros((1, parent.inputs))
                parent.fit(generic, generic)
                active_nodes = np.array(list(parent.get_active_nodes()))
            return parent.model.copy(), active_nodes

        # Ensure active nodes for both parents
        g1, m1 = ensure_active_nodes(p1)
        g2, m2 = ensure_active_nodes(p2)

        # Determine crossover point
        xover_point = determine_crossover_point(m1, m2)
        if xover_point <= 0:
            return CGP(
                model=g1,
                model_keys=self.model_keys,
                fixed_length=self.fixed_length,
                fitness_function=self.ff_string,
                mutation_type=self.mutation_type,
            )

        # Create new model by swapping sections
        g0 = np.concatenate(
            (g1[:xover_point], g2[xover_point:])
        )  # NumPy-based slicing and concatenation

        # Determine first function node for reference
        fb_node = min(p1.first_body_node, p2.first_body_node)

        # Extract active nodes surrounding the crossover point
        n_a1 = m1[m1 <= xover_point]
        n_a2 = m2[m2 > xover_point]

        if (
            len(n_a1) > 0 and len(n_a2) > 0
        ):  # Ensure both lists contain active function nodes
            n_f, n_b = n_a1[-1], n_a2[0]  # Identify boundary nodes
            g0 = neighborhood_connect(n_f, n_b, g0)  # Ensure neighborhood connectivity

        # Merge active nodes and reconnect
        n_a = np.concatenate((n_a1, n_a2))
        if len(n_a) > 0:
            g0 = random_active_connect(
                (p1.inputs + len(p1.constants)), n_a, xover_point, g0
            )

        # Create and return new CGP instance
        g0 = CGP(
            model=g0,
            model_keys=self.model_keys,
            fixed_length=self.fixed_length,
            fitness_function=self.ff_string,
            mutation_type=self.mutation_type,
        )
        g0.xover_index[xover_point - fb_node] += 1
        return g0

    def _mutate(self, models, gen, mutation_rate, verbose=False):
        if self.mutation_can_make_children:
            children = []
            for m, model in enumerate(models):
                for _ in range(self.max_c // len(models)):  # N children per parent
                    if np.random.rand() < mutation_rate or self.max_p == 1:
                        child = model.mutate(verbose)

                        # Sanity checks to catch aliasing bugs
                        assert (
                            child.id != model.id
                        ), "Mutated child has same ID as parent"
                        assert child is not model, "Child is not a distinct instance"
                        assert id(child.model) != id(
                            model.model
                        ), "Structured array not deeply copied"

                        child.fit(self.x, self.y, mutable=False)

                        child_key = f"Child_{self.child_id_counter:03d}_g{gen}"
                        parent_key = getattr(
                            model, "child_keys", f"Model_{m:03d}_g{gen - 1}"
                        )

                        child.set_parent_key([parent_key])
                        child.set_child_key(child_key)

                        if hasattr(self, "model_key_map"):
                            self.model_key_map[child_key] = child

                        self.child_id_counter += 1
                        children.append(child)

                    if len(children) >= self.max_c:
                        break

            return np.array(children, dtype=object)
        else:  # In-place mutation
            for m, model in enumerate(models):
                if np.random.rand() < mutation_rate:
                    mutated = model.mutate(verbose)
                    mutated.fit(self.x, self.y)

                    child_key = f"Child_{self.child_id_counter:03d}_g{gen}"
                    parent_key = getattr(
                        model, "child_keys", f"Model_{m:03d}_g{gen - 1}"
                    )
                    mutated.set_parent_key([parent_key])
                    mutated.set_child_key(child_key)

                    if hasattr(self, "model_key_map"):
                        self.model_key_map[child_key] = mutated

                    self.child_id_counter += 1
                    models[m] = mutated  # Overwrite in-place

            return models

    def _get_fitnesses(self, mode="train", pop_list=None, mutable=True):
        """Compute fitnesses for all models."""

        if mode == "train":
            x = self.x
            y = self.y
            f_list = self.fitnesses
            corr_list = self.corrs
        elif mode == "test":
            x = self.x_test
            y = self.y_test
            f_list = self.fitnesses_test
            corr_list = self.corr_test
        if pop_list is not None:
            f_list = deepcopy(self.fitnesses)
        if pop_list is None:
            pop_list = self.population
        for i in range(len(pop_list)):
            if pop_list[i] is not None:
                (
                    corr_list[i],
                    _,
                    f_list[i],
                ) = pop_list[
                    i
                ].fit(x, y, mutable=mutable)
        return np.array(f_list), np.array(corr_list)

    def _group_parents_and_children(self):
        individuals_with_parents = [
            model
            for model in self.population
            if model is not None
            and hasattr(model, "parent_keys")
            and model.parent_keys is not None
        ]

        parent_child_map = {}

        for child in individuals_with_parents:
            parent_pair = tuple(child.parent_keys)
            parent_child_map.setdefault(parent_pair, []).append(child)

        parent_child_groups = {
            parent_pair: (
                [
                    self.model_key_map[p_key]
                    for p_key in parent_pair
                    if p_key in self.model_key_map
                ],
                children,
            )
            for parent_pair, children in parent_child_map.items()
        }

        return parent_child_groups

    def _get_similarity_score(self, model1_obj, model2_obj):
        """
        Computes a similarity score between two models based on structural alignment.

        Args:
            model1_obj: First CGP model object (must have .model as 2D array).
            model2_obj: Second CGP model object (must have .model as 2D array).

        Returns:
            float: A similarity score (higher = more similar).
        """
        m1 = deepcopy(model1_obj.model)
        print(m1)
        m2 = deepcopy(model2_obj.model)

        # Validate dimensionality
        if not (isinstance(m1, np.ndarray) and m1.ndim == 2):
            print(f"Model1 not 2D: shape={getattr(m1, 'shape', None)}")
            return 0.0
        if not (isinstance(m2, np.ndarray) and m2.ndim == 2):
            print(f"Model2 not 2D: shape={getattr(m2, 'shape', None)}")
            return 0.0

        def _map_functions(mo1, mo2):
            unique_functions = np.unique(np.concatenate((mo1[:, 1], mo2[:, 1])))
            function_map = {func: i for i, func in enumerate(unique_functions)}
            mo1[:, 1] = np.vectorize(function_map.get)(mo1[:, 1])
            mo2[:, 1] = np.vectorize(function_map.get)(mo2[:, 1])
            return mo1, mo2

        try:
            m1, m2 = _map_functions(m1, m2)
            sequence1 = m1[:, 1:].flatten().astype(str)
            sequence2 = m2[:, 1:].flatten().astype(str)
        except Exception as e:
            print("Error during sequence mapping or flattening:", e)
            return 0.0

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -2

        try:
            seq1 = Seq("".join(sequence1))
            seq2 = Seq("".join(sequence2))
            score = aligner.score(seq1, seq2)
            return score if np.isfinite(score) else 0.0
        except Exception as e:
            print("Alignment error:", e)
            return 0.0

    def _analyze_similarity(self):
        """
        Computes the similarity between parents and their best offspring.
        Uses structural similarity to compare genetic representations.
        """
        parent_child_groups = self._group_parents_and_children()
        similarity_scores = []

        for parent_pair, (parents, children) in parent_child_groups.items():
            # Filter out None or invalid entries
            if any("g-1" in p for p in parent_pair):
                continue
            parents = [p for p in parents if p is not None]
            children = [c for c in children if c is not None]
            if not parents or not children:
                continue

            best_parent = min(parents, key=lambda ind: getattr(ind, "fitness", np.inf))
            best_child = min(children, key=lambda ind: getattr(ind, "fitness", np.inf))

            score = self._get_similarity_score(best_parent, best_child)
            similarity_scores.append((parent_pair, score))

        return similarity_scores

    def _record_metrics(self, gen: int):
        """
        Records key performance metrics for each generation, including fitness statistics
        and similarity measurements.
        """
        # Extract fitness values & active node counts
        fit_list = [
            self.fitnesses[i]
            for i in range(len(self.population))
            if self.population[i] is not None
        ]
        fit_test_list = [
            self.fitnesses_test[i]
            for i in range(len(self.population))
            if self.population[i] is not None
        ]

        corr_list = [m.correlation for m in self.population if m is not None]
        corr_test_list = [
            self.corr_test[i]
            for i in range(len(self.population))
            if self.population[i] is not None
        ]
        comp_list = [m.complexity for m in self.population if m is not None]

        active_nodes_list = []
        for p in range(len(self.population)):
            if self.population[p] is not None:
                active_nodes_list.append(self.population[p].count_active_nodes())
        active_nodes_list = np.atleast_1d(active_nodes_list)

        # Compute similarity scores
        # similarity_scores = self._analyze_similarity()
        # similarity_values = np.array([score for _, score in similarity_scores])

        # Find the best model by fitness (lower is better)
        best_model_index = np.argmin(fit_list)
        best_test_model_index = np.argmin(fit_test_list)
        self.best_model = self.population[best_model_index]

        # Compute quartile statistics efficiently
        fit_statistics = _get_quartiles(fit_list)
        fit_test_statistics = _get_quartiles(fit_test_list)

        corr_statistics = _get_quartiles(corr_list)
        corr_test_statistics = _get_quartiles(corr_test_list)
        comp_statistics = _get_quartiles((comp_list))

        active_nodes_statistics = _get_quartiles(active_nodes_list)
        semantic_diversity = np.nanstd(fit_list)
        """
        # Compute similarity quartiles (handle empty case)
        if len(similarity_values) > 0:
            similarity_quartiles = _get_quartiles(similarity_values)
        else:
            similarity_quartiles = [np.nan] * 5  # Default to NaN if no similarity data
        """
        # Store metrics in the NumPy structured array
        self.metrics[gen] = (
            *corr_statistics,
            *corr_test_statistics,
            *comp_statistics,
            *fit_statistics,
            *fit_test_statistics,
            self.best_model.count_active_nodes(),
            *active_nodes_statistics,
            semantic_diversity,  # Semantic diversity
        )

    def save_metrics(self, path=None):
        path = path if path is not None else "."
        print(f"{path}/statistics.csv")
        # Save the metrics DataFrame only once
        np.savetxt(f"{path}/statistics.csv", self.metrics, delimiter=",")

        # Save the xover_index categories efficiently
        for cat in ["deleterious", "neutral", "beneficial"]:
            np.savetxt(
                f"{path}/xover_density_{cat}.csv",
                self.xover_index[cat].astype(np.int32),
                delimiter=",",
            )

    def _report_generation(self, g: int):
        # Efficient logging instead of multiple print calls
        best_ind = self.population[np.argmin(self.fitnesses)]
        print(f"Generation {g}")
        print(
            f"Best Fitness: {np.min(self.fitnesses)}:\tCorrelation: {best_ind.correlation}\tComplexity: {best_ind.complexity}"
        )
        print("################")

    def _compare_child_parents(self):
        parent_child_groups = self._group_parents_and_children()

        for parent_pair, (parents, children) in parent_child_groups.items():
            if not parents or not children:
                continue

            # Extract fitness values from parents for quicker access
            parent_fitness = [parent.fitness for parent in parents]

            for child in children:
                if child.better_than_parents is None:
                    # Check if the child's fitness is better or worse than any parent
                    child_fitness = child.fitness
                    if any(child_fitness > f for f in parent_fitness):
                        child.better_than_parents = "deleterious"
                    elif any(child_fitness < f for f in parent_fitness):
                        child.better_than_parents = "beneficial"

    def _box_distribution(self, gen):
        """Accumulates crossover distribution statistics by type, enforcing fixed length."""
        individuals_with_parents = [
            ind
            for ind in self.population
            if ind is not None
            and ind.parent_keys is not None
            and ind.better_than_parents is not None
        ]

        xover_index = self.xover_index
        expected_len = xover_index["beneficial"].shape[1]

        for ind in individuals_with_parents:
            if ind.xover_index.shape[0] != expected_len:
                raise ValueError(
                    f"xover_index length mismatch at generation {gen}: "
                    f"expected {expected_len}, got {ind.xover_index.shape[0]}"
                )

            if ind.better_than_parents == "beneficial":
                xover_index["beneficial"][gen - 1, :] += ind.xover_index
            elif ind.better_than_parents == "deleterious":
                xover_index["deleterious"][gen - 1, :] += ind.xover_index
            else:
                xover_index["neutral"][gen - 1, :] += ind.xover_index

            ind.xover_index.fill(0)

    def set_max_gens(self, gens):
        self.max_g = gens

    def expand_generations_if_needed(self, new_max_g: int):
        if new_max_g <= self.original_max_g:
            return
        pad = new_max_g - self.original_max_g
        self.metrics = np.pad(self.metrics, ((0, pad), (0, 0)), mode="constant")
        self.mut_index = np.pad(self.mut_index, ((0, pad), (0, 0)), mode="constant")
        for k in self.xover_index:
            self.xover_index[k] = np.pad(
                self.xover_index[k], ((0, pad), (0, 0)), mode="constant"
            )

    def initialize_xover_index(self):
        """Initialize self.xover_index dynamically based on self.one_d mode."""
        """
        if self.one_d:
            flat, _ = self.flatten_parent(self.population[0])
            xover_len = len(flat)+self.population[0].outputs
            print(xover_len)
        else:
            xover_len = self.population[0].max_size + self.population[0].outputs
        """
        xover_len = len(self.population[0].xover_index)
        self.xover_index = {
            cat: np.zeros((self.max_g, xover_len))
            for cat in ["deleterious", "neutral", "beneficial"]
        }

    def _reinsert_elites(self, protected_parents):
        """Reinsert top n_elites based on fitness, ensuring integrity."""
        elite_indices = np.argsort([p.fitness for p in protected_parents])[
            : self.n_elites
        ]
        elites = [deepcopy(protected_parents[i]) for i in elite_indices]
        corrected_elites = []

        for elite in elites:
            original_fitness = elite.fitness
            elite_copy = deepcopy(elite)
            elite_copy.slope = elite.slope
            elite_copy.intercept = elite.intercept

            recomputed = elite_copy.fit(self.x, self.y, mutable=False)
            diff = abs(original_fitness - recomputed)
            if diff > 1e-5:
                print(f"ΓÜá∩╕Å Minor mismatch ({diff:.2e}) ΓÇö tolerating.")
            elif diff > 1e-2:
                raise RuntimeError(
                    f"ΓÜá∩╕Å Mismatch in elite fitness ΓÇö overwriting stored fitness: {original_fitness} ΓåÆ {recomputed}"
                )

            elite_copy.fitness = recomputed

            # print(f"[ELITISM] Re-inserting elite ID: {elite.id} (Fitness: {elite_copy.fitness})")
            corrected_elites.append(elite_copy)

        # Replace the worst individuals with elites
        worst_indices = np.argsort([ind.fitness for ind in self.population])[
            -self.n_elites :
        ]
        for idx, elite in zip(worst_indices, corrected_elites):
            self.population[idx] = elite
            self.fitnesses[idx] = elite.fitness

    def fit(
        self,
        train_x: np.ndarray,
        test_x: np.ndarray,
        train_y: np.ndarray,
        test_y: np.ndarray,
        step_size: int = None,
        xover_rate: float = 0.5,
        mutation_rate: float = 0.5,
    ):
        """
        Trains the Cartesian Genetic Programming model using evolutionary techniques.

        Args:
            train_x (np.ndarray): Training input data.
            train_y (np.ndarray): Training target data.
            step_size (int, optional): Frequency of reporting progress.
            xover_rate (float): Probability of performing crossover.
            mutation_rate (float): Probability of performing mutation.

        Returns:
            CGP: The best evolved model.
        """
        self.x, self.x_test, self.y, self.y_test = train_x, test_x, train_y, test_y
        best_fitness_train = []
        best_fitness_test = []
        # Sanity checks
        if len(self.x) < 1:
            raise ValueError("Must have at least one input value.")
        if len(self.y) != len(self.x):
            raise ValueError("Must have a 1:1 mapping for input set to output values.")
        if step_size is not None and not isinstance(step_size, int):
            raise TypeError("Step size must be either of type `int` or `None`.")
        if self.first_submission:  # first time setup
            print("First Time Setup")
            print(f"1D Xover {self.one_d}")
            self.initialize_population()

            self.initialize_xover_index()
            self.model_key_map = {}  # Add this at the beginning of fit()

            self._get_fitnesses(mode="train")
            self._get_fitnesses(mode="test")

            # Metrics and tracking
            self.metrics = np.zeros((self.max_g + 1, 32), dtype=np.float64)
            self.xover_index = {
                cat: np.zeros((self.max_g, self.max_p))
                for cat in ["deleterious", "neutral", "beneficial"]
            }
            self.mut_index = np.zeros((self.max_g, self.max_p))

            self._record_metrics(0)
            self._report_generation(0)

            # Setup mutation and crossover tracking
            genes_per_instruction = self.population[0].arity + 1  # for the operator
            model_size = len(self.xover_index)
            if self.one_d:
                self.mut_index = np.zeros((self.max_g, model_size))
            else:
                self.mut_index = np.zeros(
                    (self.max_g, model_size * genes_per_instruction)
                )

        else:
            self.expand_generations_if_needed(self.max_g)

        self._get_fitnesses(mode="test", mutable=False)
        self._get_fitnesses(mode="train", mutable=True)

        n_elites = (
            self.n_elites if hasattr(self, "n_elites") else 1
        )  # or pass as parameter
        # Generation 0: Just record metrics, no elite reinsertion
        self.current_generation += (
            1 if self.current_generation else self.current_generation
        )
        self._record_metrics(self.current_generation)
        self._report_generation(self.current_generation)
        elite_prev = self.elite_selection(n_elites=n_elites)

        for gen in range(self.current_generation + 1, self.max_g + 1):
            self.current_generation = gen
            # **Parent Selection**
            # self._get_fitnesses(mutable=False, mode='test')
            # self._get_fitnesses(mutable=True, mode='train')

            selected_parents = [deepcopy(p) for p in self.selection()]

            # **Crossover to Generate Children**
            if self.xover:
                children = self.crossover(selected_parents, xover_rate, gen)
                child_fitnesses = self._get_fitnesses(
                    pop_list=children, mutable=False, mode="train"
                )
                self._compare_child_parents()

            else:
                children = deepcopy(
                    selected_parents
                )  # No crossover, pass parents as is

            if not isinstance(selected_parents, (list, np.ndarray)):
                selected_parents = [selected_parents]

            #    # Clone selected parents
            cloned_parents = [deepcopy(p) for p in selected_parents]
            protected_parents = [deepcopy(p) for p in cloned_parents]

            for i, (orig, protected) in enumerate(
                zip(selected_parents, protected_parents)
            ):
                assert not np.shares_memory(
                    orig.model, protected.model
                ), f"Memory shared at index {i}"

            mutated_children = self._mutate(children, gen, mutation_rate)

            if self.mutation_can_make_children:
                assert all(
                    c.id != p.id for c in mutated_children for p in selected_parents
                ), "Mutation may be in-place!"

            # Ensure both are proper lists of CGP instances
            if isinstance(protected_parents, CGP):
                protected_parents = [protected_parents]
            elif isinstance(protected_parents, np.ndarray):
                protected_parents = list(protected_parents)

            if isinstance(mutated_children, CGP):
                mutated_children = [mutated_children]
            elif isinstance(mutated_children, np.ndarray):
                mutated_children = list(mutated_children)

            # Final sanity check
            assert all(
                isinstance(p, CGP) for p in protected_parents
            ), "Non-CGP in protected_parents"
            assert all(
                isinstance(c, CGP) for c in mutated_children
            ), "Non-CGP in mutated_children"
            assert all(
                p is not None for p in protected_parents
            ), "protected_parents contains None"
            assert all(
                c is not None for c in mutated_children
            ), "mutated_children contains None"

            self.population = protected_parents + mutated_children

            self._get_fitnesses(mutable=False, mode="test")
            self._get_fitnesses(mutable=True, mode="train")
            """
            # Γ£à Step 2: Save elite for next generation
            if elite_prev and 'elite' in self.selection_type:
                fitnesses = [m.fitness for m in self.population if m is not None]
                worst_indices = np.argsort(fitnesses)[-n_elites:]
                for i, idx in enumerate(worst_indices):
                    print(i)
                    print(elite_prev)
                    self.population[idx] = deepcopy(elite_prev[i])

                # print(f"[Gen {gen}] Reinserted {n_elites} elite(s) with fitnesses:",
                #       [f"{e.fitness:.4e}" for e in elite_prev])

            elite_prev = self.elite_selection(n_elites=n_elites)
            # for e in elite_prev:
            #    print(f"Elite fitness: {e.fitness}")
            """

            assert all(
                isinstance(e, CGP) for e in elite_prev
            ), "Elite_selection returned non-CGPs"
            assert all(
                hasattr(e, "fitness") and e.fitness is not None for e in elite_prev
            ), "Elite has no fitness"

            true_elite_train = min(self.population, key=lambda x: x.fitness)
            # true_elite_test = self.population[np.argmin(self.fitnesses_test)]
            best_fitness_train.append(true_elite_train.fitness)
            best_fitness_test.append(
                self.fitnesses_test[np.argmin(self.fitnesses_test)]
            )

            self._record_metrics(gen)

            if step_size and gen % step_size == 0:
                self._report_generation(gen)
                self.save_checkpoint(filename=self.ckpt_filename, generation=gen)

        return (
            self.population[np.argmin(self.fitnesses)],
            self.population[np.argmin(self.fitnesses_test)],
        )

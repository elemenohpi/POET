import numpy as np
import uuid
import hashlib
from .cgp_generator import generate_model, node_to_int
from .cgp_operators import add, sub, mul, div
from .fitness_functions import correlation, align, corr_comp_fitness
from copy import deepcopy

class CGP:

    def __init__(self, model=None, model_keys=None, fixed_length=True, fitness_function='Correlation',
                 mutation_type='Point',
                 parent_keys=None, xover_length=None, **kwargs):
        self.correlation = 1.0
        self.complexity = 1.0
        self.id = uuid.uuid4()
        self.mutation = None
        self.slope = None
        self.intercept = None
        self.fitness = None
        self.fixed_length = fixed_length
        self.parent_keys = parent_keys
        self.child_keys = None
        self.better_than_parents = None

        # Assign fitness function
        if fitness_function.lower() == 'correlation':
            self.fitness_function = correlation
        elif fitness_function.lower() == 'correlation_complexity':
            self.fitness_function = corr_comp_fitness
        else:
            raise ValueError(f'Invalid Fitness Function: {fitness_function}')

        # Function bank setup
        self.function_bank = kwargs.get('function_bank', {'add': add, 'sub': sub, 'mul': mul, 'div': div})
        self.function_bank = {mapping: key for mapping, key in enumerate(self.function_bank.values())}
        self.n_operations = len(self.function_bank)
        assert self.n_operations >= 1, 'At least one operator is required.'

        # Initialize model
        if model is not None and model_keys is not None:
            self.model = model
            self.model_keys = model_keys
            self._initialize_from_model()
        elif (model is not None and model_keys is None) or (model is None and model_keys is not None):
            raise RuntimeError('Must specify both model and model_keys.')
        else:
            self._initialize_from_kwargs(kwargs)

        self.first_body_node = self.inputs + np.sum(
            self.model[:, self.model_keys['NodeType']] == node_to_int('Constant'))
        self.last_body_node = self.first_body_node + self.max_size - 1
        # Γ£à Initialize xover_index correctly
        if xover_length is not None:
            self.xover_index = np.zeros(xover_length)
        else:
            self.xover_index = np.zeros(self.max_size + self.outputs)

        self._choose_mutation(mutation_type)

    def _initialize_from_model(self):
        """Initialize attributes from a pre-built model."""
        self.constants = self.model[:, self.model_keys['Value']][
            self.model[:, self.model_keys['NodeType']] == node_to_int('Constant')
            ]
        self.inputs = np.sum(self.model[:, self.model_keys['NodeType']] == node_to_int('Input'))
        self.outputs = np.sum(self.model[:, self.model_keys['NodeType']] == node_to_int('Output'))
        self.max_size = np.sum(self.model[:, self.model_keys['NodeType']] == node_to_int('Function'))
        self.n_operations = len(set(self.model[:, self.model_keys['Operator']]))
        self.arity = sum(name.startswith('Operand') for name in self.model_keys)

    def _initialize_from_kwargs(self, kwargs):
        """Initialize attributes from keyword arguments."""
        # self.constants = np.array(kwargs.get('constants', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        self.constants = np.atleast_1d(kwargs.get('constants', [1]))
        self.inputs = kwargs.get('inputs', 1)
        self.outputs = kwargs.get('outputs', 1)
        self.arity = kwargs.get('arity', 2)
        self.max_size = kwargs.get('max_size', 16)

        assert self.inputs >= 1, 'There must be at least one input feature.'
        assert self.outputs >= 1, 'There must be at least one output.'
        assert self.arity >= 1, 'Operators must take at least one input.'
        assert self.max_size >= 1, 'There must be at least one instruction.'

        # Generate the model using NumPy instead of Pandas
        self.model, self.model_keys = generate_model(
            self.max_size, self.inputs, self.constants, self.arity, self.outputs,
            self.n_operations, self.function_bank, self.fixed_length
        )

    def __call__(self, data, mutable=True):
        data = np.atleast_2d(data)
        self.visited = set()
        if not mutable:
            model = np.copy(self.model)
            assert not np.shares_memory(model, self.model), "Γ¥î Model and self.model share memory!"
        else:
            model = self.model

        outputs = np.array([
            self._compute_single_input(datum, model, mutable) for datum in data
        ])

        if self.slope is not None and self.intercept is not None:
            outputs = outputs * self.slope + self.intercept
        return outputs

    def _get_node_value(self, model, operand, mutable, visited=None):
        """Compute the value of a node, with cycle detection."""
        if self.visited is None:
            self.visited = set()

        operand = int(operand)

        if operand in visited:
            raise RuntimeError(f"Cycle detected at node {operand} ΓÇö already visited")


        try:
            node = model[int(operand)]
        except IndexError as e:
            print(f'_get_node_value(): {e}')
            print(f'operand: {operand}')
            print(model)
            print(f'model size: {model.shape}')
            exit()

        node_type = node[self.model_keys['NodeType']]

        if node_type in map(node_to_int, ['Input', 'Constant']):
            return node[self.model_keys['Value']]

        elif node_type == node_to_int('Function'):
            self.visited.add(operand)
            operand_values = np.array([
                self._get_node_value(model, node[self.model_keys[f'Operand{i}']], mutable, visited.copy())
                for i in range(self.arity)
            ])
            operator = self.function_bank[node[self.model_keys["Operator"]]]
            result = operator(*operand_values)
            if not np.isfinite(result):
                print(f'Warning: result of operation on {operand_values} is infinite or invalid. Returning 0.0')
                result = 0.0

            if mutable:
                try:
                    model[operand, self.model_keys['Value']] = result
                except OverflowError as e:
                    print(
                        f'Cannot cast {result}\noperand: {operand}\noperand values: {operand_values}\noperator: {operator}\nReturning np.inf')
            model[operand, self.model_keys['Active']] = 1

            return result

        else:
            print(node)
            raise ValueError(f"Invalid node type: {node_type}")

    def _run(self, model, mutable):
        """Run the model and compute output values."""
        model[:, self.model_keys['Active']] = 0  # Reset active flags
        if mutable:
            model[model[:, self.model_keys['NodeType']] == 'Function', self.model_keys['Value']] = 0

        output_indices = list(range(len(model) - self.outputs, len(model)))

        output_values = np.empty(self.outputs)

        for i, idx in enumerate(output_indices):
            output_values[i] = self._get_node_value(model, model[idx, self.model_keys['Operand0']], mutable,
                                                    visited=set())

        if mutable:
            model[output_indices, self.model_keys['Value']] = output_values

        return output_values

    def _compute_single_input(self, datum, model, mutable):
        input_indices = list(range(self.inputs))
        if not mutable:
            model = model.copy()
            model[:, self.model_keys['Value']] = model[:, self.model_keys['Value']].copy()
        model[input_indices, self.model_keys['Value']] = datum
        return self._run(model, mutable)

    def fit(self, data, ground_truth, mutable=True):
        predictions = self.__call__(data, mutable=mutable)

        # Compare raw unaligned outputs
        # pred1 = self.__call__(data, mutable=False)
        # pred2 = self.__call__(data, mutable=False)
        # assert np.allclose(pred1, pred2), "Γ¥î Model output changed between evaluations!"

        n_active_nodes = self.count_active_nodes()

        # Compute fitness using raw predictions
        self.correlation, self.complexity, self.fitness = self.fitness_function(predictions, ground_truth, n_active_nodes, float(self.max_size))
        #fitness_check = self.fitness_function(predictions, ground_truth)
        #if not np.isclose(self.fitness, fitness_check, atol=1e-8):
        #    raise RuntimeError(f"Fitness changed on re-computation: {self.fitness} vs {fitness_check}")

        # Now align (only for prediction)
        if self.fitness_function == correlation:
            self.slope, self.intercept = align(predictions, ground_truth)

        return self.correlation, self.complexity, self.fitness

    def get_active_nodes(self):
        return self.visited

    def count_active_nodes(self):
        return len(self.visited)

    def _choose_mutation(self, mutation_type):
        """Assign mutation function."""
        mutation_type = mutation_type.lower()
        mutations = {'point': self._point_mutation, 'full': self._full_mutation}
        if mutation_type not in mutations:
            raise ValueError(f'{mutation_type} is an invalid mutation operator.')
        self.mutation = mutations[mutation_type]

    def mutate_output(self):
        active_indices = np.where((self.model[:, self.model_keys['NodeType']] == node_to_int('Output')))[0]
        mutation_index = np.random.choice(active_indices)
        new_node = [node_to_int('Output'), 0, 0, np.random.randint(0, self.first_body_node + self.max_size),
                    *[0 for _ in range(self.arity - 1)], 0]
        self.model[mutation_index] = deepcopy(new_node)

    def _full_mutation(self, verbose=True):
        active_indices = np.where((self.model[:, self.model_keys['NodeType']] == node_to_int('Function')) | (
                self.model[:, self.model_keys['NodeType']] == node_to_int('Output')))[0]

        # Select a random function node
        mutation_index = np.random.choice(active_indices)
        if verbose:
            print(f'Mutating at index {mutation_index}')
        if self.model[mutation_index, self.model_keys['NodeType']] == node_to_int('Function'):
            new_node = [node_to_int('Function'), 0, np.random.choice(list(self.function_bank.keys())),
                        *[np.random.randint(0, mutation_index) for _ in range(self.arity)], 1]
            if verbose:
                print(f'Replacing {self.model[mutation_index]} at index {mutation_index} to {new_node}')
            self.model[mutation_index] = deepcopy(new_node)
        else:
            new_node = [node_to_int('Output'), 0, 0, np.random.randint(0, self.first_body_node + self.max_size),
                        *[0 for _ in range(self.arity - 1)], 0]
            if verbose:
                print(f'Replacing {self.model[mutation_index]} at index {mutation_index} to {new_node}')
            self.model[mutation_index] = deepcopy(new_node)

    def _point_mutation(self, verbose=False):
        # active_indices = np.where((self.model['NodeType'] == 'Function') & (self.model['Active'] == 1))[0]
        active_indices = np.where((self.model[:, self.model_keys['NodeType']] == node_to_int('Function')) | (
                self.model[:, self.model_keys['NodeType']] == node_to_int('Output')))[0]

        # if len(active_indices) == 0:
        #    # Fallback: Mutate any function node if no active ones exist
        #    active_indices = np.where((self.model[:, self.model_keys['NodeType']] == node_to_int('Function)) | (self.model['NodeType'] == 'Output'))[0]

        # Select a random function node
        mutation_index = np.random.choice(active_indices)
        if verbose:
            print(f'Mutating at index {mutation_index}')
        if self.model[mutation_index, self.model_keys['NodeType']] == node_to_int('Function'):
            mutation_column = np.random.choice(
                [self.model_keys['Operator']] + [self.model_keys[f'Operand{i}'] for i in range(self.arity)]
            )
            if mutation_column == self.model_keys['Operator']:
                current_op = self.model[mutation_index, self.model_keys['Operator']]
                ops = list(self.function_bank.keys())
                new_operator = np.random.choice(ops)
                while new_operator == current_op:
                    new_operator = np.random.choice(ops)
                if verbose:
                    print(f'Mutating column {mutation_column} to {new_operator}')
                self.model[mutation_index, self.model_keys['Operator']] = new_operator
            else:
                # Mutate operand, ensuring a different value
                new_operand = np.random.randint(0, mutation_index)
                while new_operand == self.model[mutation_index, mutation_column]:
                    new_operand = np.random.randint(0, mutation_index)
                if verbose:
                    print(f'Mutating column {mutation_column} to {new_operand}')
                self.model[mutation_index, mutation_column] = new_operand
        else:  # mutate output node
            old_operand = self.model[mutation_index, self.model_keys['Operand0']]
            new_operand = old_operand
            while old_operand == new_operand:
                new_operand = np.random.randint(0, mutation_index)
            if verbose:
                print(f'Mutating output to {new_operand}')
            self.model[mutation_index, self.model_keys['Operand0']] = new_operand

    def mutate(self, verbose=False):
        """Return a mutated copy of the individual."""
        clone = deepcopy(self)
        clone.id = uuid.uuid4()
        clone.model = self.model.copy()
        if clone.mutation is None:
            raise ValueError("Mutation function not set. Call _choose_mutation() first.")
        clone.mutation(verbose)
        return clone

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            if key == 'model':
                result.model = np.copy(value)
            elif key == 'model_keys':
                result.model_keys = deepcopy(value, memo)
            else:
                setattr(result, key, deepcopy(value, memo))
        result.fitness = self.fitness
        return result

    def print_parameters(self):
        """Print key parameters of the CGP model."""
        print(f"Constants: {self.constants}")
        print(f"Inputs: {self.inputs}")
        print(f"Outputs: {self.outputs}")
        print(f"Arity: {self.arity}")
        print(f"Max Instructions: {self.max_size}")
        print(f"Function Bank: {[func.__name__ for func in self.function_bank]}")
        print(f"Number of Functions: {self.n_operations}")

    def print_model(self):
        print(self.model)

    def set_parent_key(self, key):
        self.parent_keys = key

    def set_child_key(self, key):
        self.child_keys = key

    @staticmethod
    def _hash_model(model: np.ndarray) -> str:
        """Return a hash of the model's raw bytes for integrity checking."""
        return hashlib.md5(model.tobytes()).hexdigest()

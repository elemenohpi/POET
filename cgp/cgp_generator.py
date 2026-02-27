import numpy as np
from .cgp_operators import *


def int_to_node(i):
    l = ['Input', 'Constant', 'Function', 'Output']
    return l[i]


def node_to_int(node):
    l = ['Input', 'Constant', 'Function', 'Output']
    try:
        return l.index(node)
    except ValueError:
        raise ValueError(
            f'Tried to convert node {node} to an integer. The only valid types are \"Input\", \"Constant\", \"Function\", and \"Output\"')
    except RecursionError:
        raise RecursionError(f'Recursion Error in cgp_generator.py::node_to_int\tTried to call {node}')


def generate_model(max_size: int, inputs: int, constants: list | np.ndarray, arity: int, outputs: int,
                   n_operations: int, function_bank: dict, fixed_length: bool = True):
    """
    Generates a CGP model using NumPy arrays instead of pandas.

    Args:
        max_size (int): Maximum number of body nodes.
        inputs (int): Number of variable inputs.
        constants (list | np.ndarray): List of constants.
        arity (int): Number of arguments per function node.
        outputs (int): Number of outputs.
        n_operations (int): Number of operations in the function bank.
        function_bank (tuple): Tuple of available operations.
        fixed_length (bool): If False, the model size can be <= max_size.

    Returns:
        np.ndarray: NumPy array representing the model.
        list: Names of Model Columns

    Structured numpy arrays don't work because of copy issues
    """

    constants = np.array(constants) if isinstance(constants, list) else constants

    # Define model keys
    model_keys = ['NodeType', 'Value', 'Operator', *[f'Operand{i}' for i in range(arity)], 'Active']
    num_keys = len(model_keys)
    model_keys = {key: i for i, key in enumerate(model_keys)}

    # Input and Constant Nodes
    num_constants = len(constants)
    num_inputs = inputs + num_constants
    input_nodes = np.zeros((num_inputs, num_keys))
    input_nodes[:inputs, model_keys['NodeType']] = node_to_int('Input') #should be 0
    input_nodes[inputs:, model_keys['NodeType']] = node_to_int('Constant')
    input_nodes[inputs:, model_keys['Value']] = constants  # Assign constant values

    # Determine the actual number of function nodes
    model_size = max_size if fixed_length else np.random.randint(1, max_size)
    first_body_node = num_inputs
    last_body_node = first_body_node + model_size

    # Function Nodes
    body_nodes = np.zeros((model_size, num_keys))
    body_nodes[:, model_keys['NodeType']] = node_to_int('Function')
    body_nodes[:, model_keys['Operator']] = np.random.randint(0, len(function_bank), model_size)
    # Generate operands correctly
    for i in range(arity):
        body_nodes[:, model_keys[f'Operand{i}']] = np.random.randint(0, first_body_node + np.arange(model_size), size=model_size)

    # Output Nodes
    output_nodes = np.zeros((outputs, num_keys))
    output_nodes[:, model_keys['NodeType']] = node_to_int('Output')
    output_nodes[:, model_keys['Operand0']] = np.random.randint(0, last_body_node, outputs)

    # Combine all nodes into a single structured NumPy array
    model = np.concatenate((input_nodes, body_nodes, output_nodes), axis=0)

    return model, model_keys

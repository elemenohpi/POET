import numpy as np
from .cgp_model import CGP
from .cgp_generator import node_to_int
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def _get_quartiles(data):
    data = np.where(data == np.inf, 1, data)
    return np.quantile(data, [0, 0.25, 0.5, 0.75, 1])


def _get_semantic_alignment(model1: CGP, model2: CGP, x_train, cumulative: bool = False):
    """
    Compute semantic alignment weights between two CGP models.

    This function calculates a set of weights representing the semantic similarity between 
    corresponding nodes of two Cartesian Genetic Programming (CGP) models. The alignment 
    can be based either on the final output semantics (cumulative) or on intermediate 
    node-level semantics.

    Parameters:
        model1 (CGP): The first CGP model.
        model2 (CGP): The second CGP model.
        cumulative (bool, optional): If True, compute alignment based on final output semantics 
                                     (1D vector). If False, use intermediate semantics for all 
                                     active nodes (2D matrix). Defaults to False.

    Returns:
        np.ndarray: An array of alignment weights indicating semantic similarity between nodes.
    """

    # Get internal semantics (n_nodes, n_samples)
    if cumulative:
        semantics_1 = np.array([model1(x) for x in x_train])  # shape: (n_samples,) -> will need reshape
        semantics_2 = np.array([model2(x) for x in x_train])
        semantics_1 = semantics_1.T.reshape(1, -1)  # shape: (1, n_samples)
        semantics_2 = semantics_2.T.reshape(1, -1)
    else:
        semantics_1 = clean_values(model1, x_train)  # shape: (n_nodes, n_samples)
        semantics_2 = clean_values(model2, x_train)

    # Align nodes using optimal assignment
    cost_matrix = cdist(semantics_1, semantics_2, metric='euclidean')
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Align semantics according to optimal mapping
    aligned_1 = semantics_1[row_ind]
    aligned_2 = semantics_2[col_ind]

    # Compute semantic difference (SSD per node)
    ssd = get_ssd(aligned_1, aligned_2)  # shape: (n_nodes,)

    # Map SSD to crossover weights
    weights = get_weights(ssd, alpha=1e-4, beta=0.4, epsilon=0.0)  # shape: (n_nodes,)
    return weights


def _validate_int_param(param_name, value, min_val, max_val):
    """
    Validates and converts a parameter to an integer if necessary.
    Raises errors or issues warnings as needed.
    """
    if not isinstance(value, int):
        print(f"Warning: {param_name} {value} is not an integer. Rounding to {int(value)}.")
        value = int(value)
    assert min_val <= value <= max_val, (
        f"{param_name} must be between {min_val} and {max_val}. Current value: {value}."
    )
    return value


def _split(m, points: [int]):
    previous_point = 0
    parts = []
    points = np.append(points, len(m))
    for point in points:
        parts.append(m[previous_point:point])
        previous_point = point

    return parts


def pairwise_minkowski_distance(xa, xb, p=2):
    """
    Computes pairwise Minkowski distance between two collections of 1D arrays.
    :param xa: First 2D array or 1D array
    :param xb: Second 2D array or 1D array
    :param p: Minkowski distance parameter
    :return: Pairwise Minkowski distances as a 2D array
    """
    xa = np.atleast_2d(xa)  # Ensure XA is 2D
    xb = np.atleast_2d(xb)  # Ensure XB is 2D
    distance = np.sum(np.abs(xa[:, None] - xb) ** p, axis=-1) ** (1 / p)
    return distance[0, 0]


def get_score(p_to_t, p_to_m, m_to_t):
    """

    @param p_to_t: distance from parent to target
    @param p_to_m: distance from potential mate to parent
    @param m_to_t: distance from potential mate to target
    """
    if p_to_m == 0:
        return np.inf
    return (m_to_t / p_to_m) * (1 + np.abs(m_to_t - p_to_t))


def clean_values(model, x_train, include_output=False):
    """
    Gets the values matrix for each input.

    @param model: CGP model
    @param x_train: NumPy array of x inputs
    @param include_output: Whether output values are included in the set
    @return: Matrix of values with rows as node numbers and columns as inputs
    """
    # Get boolean masks for function and output nodes
    function_mask = model.model[:, model.model_keys['NodeType']] == node_to_int('Function')
    output_mask = model.model[:, model.model_keys['NodeType']] == node_to_int('Output')
    # Compute matrix size
    model_length = np.sum(function_mask)
    if include_output:
        model_length += np.sum(output_mask)

    number_of_inputs = x_train.shape[0]
    values_matrix = np.zeros((model_length, number_of_inputs))

    # Iterate over inputs
    for i, input_value in enumerate(x_train):
        model(input_value)  # Run the model with the input

        # Extract function node values
        function_values = model.model[function_mask, model.model_keys['Value']]

        if include_output:
            # Extract output node values
            output_values = model.model[output_mask, model.model_keys['Value']]
            values = np.concatenate([function_values, output_values])
        else:
            values = function_values

        values_matrix[:, i] = values

    # Ensure all values are finite, replacing NaNs/Infs with 0
    values_matrix[~np.isfinite(values_matrix)] = 0

    return values_matrix


"""
Uy, N. Q., Hien, N. T., Hoai, N. X., & OΓÇÖNeill, M. (2010). Improving the generalisation ability of genetic
 programming with semantic similarity based crossover. In Genetic Programming: 13th European Conference, 
 EuroGP 2010, Istanbul, Turkey, April 7-9, 2010. Proceedings 13 (pp. 184-195). Springer Berlin Heidelberg.
"""


# instead of picking a point at random, we get the SSD of all points and then later can pick N according to operator
# , alpha=1e-4, beta=0.4
def get_ssd(v_matrix_1, v_matrix_2):
    assert v_matrix_1.shape == v_matrix_2.shape, "Internal Semantic Matrices must be the same size"
    return np.sum(np.abs(v_matrix_1 - v_matrix_2), axis=1) / v_matrix_1.shape[1]


def get_weights(ssd_matrix: np.ndarray, alpha: float = 1e-4, beta: float = 0.4, epsilon: float = 0.0):
    # Set weights to zero where conditions are met
    weight_matrix = np.where((ssd_matrix < alpha) | (ssd_matrix >= beta), epsilon, ssd_matrix)
    total = np.sum(weight_matrix)

    # Normalize weights or use uniform distribution if the sum is zero
    return weight_matrix / total if total != 0 else np.full_like(ssd_matrix, 1 / ssd_matrix.size)



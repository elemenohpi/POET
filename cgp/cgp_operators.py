import numpy as np


def add(x, y):
    return x + y


def sub(x, y):
    with np.errstate(invalid='ignore'):
        return x - y


def mul(x, y):
    return x * y


def div(x, y):
    # Compute denominator safely
    denom = np.sqrt(1 + y ** 2)

    # Replace infinities or NaNs in denominator with a small constant
    denom = np.where(np.isfinite(denom), denom, 1e-8)

    # Avoid division by (near) zero
    denom = np.where(denom < 1e-8, 1e-8, denom)

    # Compute result
    result = x / denom

    # Replace any NaNs or infs in result with 0.0
    result = np.where(np.isfinite(result), result, 0.0)

    return result

def op_and(x, y):
    return x & y

def op_or(x, y):
    return x | y

def op_not(x, y):
    return not x

def op_xor(x, y):
    return x ^ y

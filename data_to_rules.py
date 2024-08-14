#written with chatgpt
import pandas as pd
import numpy as np
import sys

def split_sequence(sequence, max_size, min_size):
    """Splits a sequence into substrings of a size between min_size and max_size."""
    substrings = []
    start = 0
    while start < len(sequence):
        size = min(max_size, len(sequence) - start)
        if size < min_size:
            break
        substrings.append(sequence[start:start + size])
        start += size
    return substrings

def process_sequences(sequences, max_size, min_size):
    """Processes a list of sequences and returns a list of substrings."""
    rules = []
    for sequence in sequences:
        rules.extend(split_sequence(sequence, max_size, min_size))
    return rules

def initialize_arrays(length):
    """Initializes weights and status arrays."""
    weights = np.round(np.random.random(length)*10, 2)
    status = np.ones(length)
    return weights, status

def create_dataframe(patterns, weights, status):
    """Creates a DataFrame from patterns, weights, and status arrays."""
    return pd.DataFrame({
        'pattern': patterns,
        'weight': weights,
        'status': status
    })

def main():
    dataset = pd.read_csv(sys.argv[1])
    max_size = int(sys.argv[2])
    min_size = int(sys.argv[3])
    output_csv = sys.argv[4]+".csv"

    sequences = dataset['sequence']
    rules = process_sequences(sequences, max_size, min_size)

    weights, status = initialize_arrays(len(rules))
    df = create_dataframe(rules, weights, status)

    df.to_csv(output_csv, index=False)
    print(f"DataFrame saved to {output_csv}")

if __name__ == "__main__":
    main()

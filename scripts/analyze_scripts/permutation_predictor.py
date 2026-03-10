"""Generate novel high-fitness protein sequences using a POET model.

Creates 1000 random permutations of the dataset's amino-acid composition,
scores them with the model, then runs an evolutionary loop for n generations
keeping the top 15 each round.  Mutations are position swaps so every
candidate is always a valid permutation.

Usage:
  python scripts/analyze_scripts/permutation_predictor.py <model_csv> [options]

Examples:
  python scripts/analyze_scripts/permutation_predictor.py output/Mar-02-2026-upar_epoch_0/models/cleaned/model_40.csv -c configs/exp_charclass_vargaps.ini

  # 200 generations, save results
  python scripts/analyze_scripts/permutation_predictor.py output/.../model_40.csv -c configs/exp_charclass_vargaps.ini -n 200 -o top15.csv
"""

import argparse
import os
import sys
import random
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import eletility
import fitness as F
import individual as I


def load_model_and_alignment(model_path, config, data_path=None):
    """Load a model and compute the linear alignment from training data.

    Returns (Individual, Fitness, slope, intercept).
    """
    if data_path is not None:
        config = dict(config)
        config["learn_data"] = data_path

    ind = I.Individual(config)
    ind.makeFromFile(model_path)

    objF = F.Fitness(config)

    # Get raw predictions on training data to derive alignment
    objF.resetIndividual(ind)
    raw_preds = []
    for seq, actual in zip(objF.sequences, objF.fitness_values):
        _, pred = objF.eval(seq, actual, ind, True)
        raw_preds.append(pred)

    align = np.polyfit(raw_preds, objF.fitness_values, 1)
    return ind, objF, float(align[0]), float(align[1])


def predict(seq, ind, objF, slope, intercept):
    """Score a single sequence and return the aligned predicted fitness."""
    objF.resetIndividual(ind)
    _, raw = objF.eval(seq, 0.0, ind, True)
    return raw * slope + intercept


def generate_permutation(amino_acids):
    """Return a random permutation of the given amino-acid list as a string."""
    perm = list(amino_acids)
    random.shuffle(perm)
    return "".join(perm)


def mutate(seq, num_swaps=1):
    """Mutate a sequence by swapping random positions."""
    chars = list(seq)
    for _ in range(num_swaps):
        i, j = random.sample(range(len(chars)), 2)
        chars[i], chars[j] = chars[j], chars[i]
    return "".join(chars)


def main():
    parser = argparse.ArgumentParser(
        description="Predict novel high-fitness protein permutations."
    )
    parser.add_argument("model", help="Path to the model CSV file.")
    parser.add_argument(
        "-c", "--config", default="config.ini",
        help="Path to POET config file (default: config.ini)."
    )
    parser.add_argument(
        "-d", "--data", default=None,
        help="Override dataset CSV (for alignment and exclusion)."
    )
    parser.add_argument(
        "-n", "--generations", type=int, default=100,
        help="Number of evolutionary generations (default: 100)."
    )
    parser.add_argument(
        "-p", "--population", type=int, default=1000,
        help="Population size per generation (default: 1000)."
    )
    parser.add_argument(
        "-k", "--keep", type=int, default=15,
        help="Number of top sequences to keep each generation (default: 15)."
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=None,
        help="Random seed for reproducibility."
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Save final top sequences to this CSV file."
    )
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    # Load config
    configparser = eletility.ConfigParser()
    config = configparser.read(args.config)

    # Load model and compute alignment
    print("Loading model and computing alignment...")
    ind, objF, slope, intercept = load_model_and_alignment(
        args.model, config, args.data
    )

    # Determine amino-acid composition from the dataset
    dataset_seqs = set(objF.sequences)
    ref_seq = objF.sequences[0]
    amino_acids = sorted(ref_seq)
    print(f"Amino-acid composition: {''.join(amino_acids)} (from '{ref_seq}')")
    print(f"Dataset sequences excluded: {len(dataset_seqs)}")
    print(f"Population: {args.population}  |  Keep: {args.keep}  |  "
          f"Generations: {args.generations}")
    print()

    # Verify all dataset sequences share the same composition
    ref_counter = Counter(ref_seq)
    for s in objF.sequences:
        if Counter(s) != ref_counter:
            print(f"WARNING: '{s}' has different composition than '{ref_seq}'. "
                  "Falling back to generating from '{ref_seq}' composition.")
            break

    # Generate initial population (excluding dataset sequences)
    print("Generating initial population...")
    population = set()
    attempts = 0
    max_attempts = args.population * 20
    while len(population) < args.population and attempts < max_attempts:
        seq = generate_permutation(amino_acids)
        if seq not in dataset_seqs:
            population.add(seq)
        attempts += 1
    population = list(population)
    print(f"  Created {len(population)} unique sequences")

    # Score initial population
    print("Scoring initial population...")
    scored = []
    for seq in population:
        fit = predict(seq, ind, objF, slope, intercept)
        scored.append((seq, fit))
    scored.sort(key=lambda x: x[1], reverse=True)

    # Print initial top
    print(f"\n--- Initial top {args.keep} ---")
    for rank, (seq, fit) in enumerate(scored[:args.keep], 1):
        print(f"  {rank:>3}.  {seq}  {fit:.6f}")

    # Evolutionary loop
    for gen in range(1, args.generations + 1):
        # Keep the elite
        elite = scored[:args.keep]
        elite_seqs = {s for s, _ in elite}

        # Generate children by mutating the elite
        children_per_parent = args.population // args.keep
        new_pop = set()
        for parent_seq, _ in elite:
            new_pop.add(parent_seq)  # keep parent
            child_attempts = 0
            while len([s for s in new_pop if s not in elite_seqs]) < children_per_parent * (list(elite_seqs).index(parent_seq) + 1 if parent_seq in elite_seqs else 1):
                child = mutate(parent_seq, num_swaps=random.randint(1, 3))
                if child not in dataset_seqs:
                    new_pop.add(child)
                child_attempts += 1
                if child_attempts > children_per_parent * 5:
                    break

        # Fill up to population size from additional mutations if needed
        fill_attempts = 0
        while len(new_pop) < args.population and fill_attempts < args.population * 5:
            parent_seq = random.choice(elite)[0]
            child = mutate(parent_seq, num_swaps=random.randint(1, 3))
            if child not in dataset_seqs:
                new_pop.add(child)
            fill_attempts += 1

        # Score
        scored = []
        for seq in new_pop:
            fit = predict(seq, ind, objF, slope, intercept)
            scored.append((seq, fit))
        scored.sort(key=lambda x: x[1], reverse=True)

        if gen % 10 == 0 or gen == 1 or gen == args.generations:
            best = scored[0]
            worst_elite = scored[min(args.keep - 1, len(scored) - 1)]
            print(f"  Gen {gen:>4}  |  best: {best[1]:.6f}  |  "
                  f"top-{args.keep} worst: {worst_elite[1]:.6f}  |  "
                  f"pool: {len(scored)}")

    # Final results
    top = scored[:args.keep]
    print(f"\n{'='*60}")
    print(f"  Top {args.keep} predicted sequences after {args.generations} generations")
    print(f"{'='*60}")
    for rank, (seq, fit) in enumerate(top, 1):
        print(f"  {rank:>3}.  {seq}  predicted_fitness={fit:.6f}")

    # Save to CSV
    if args.output:
        df = pd.DataFrame({
            "rank": list(range(1, len(top) + 1)),
            "sequence": [s for s, _ in top],
            "predicted_fitness": [round(f, 6) for _, f in top],
        })
        df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")


if __name__ == "__main__":
    main()

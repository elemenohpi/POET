"""Rank any-to-any chimera candidates with a trained token-mode POET model."""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import chimera_codec as CC
import eletility
import fitness as F
import individual as I
import sequence_codec as SC


def _load_model_and_alignment(model_path, config):
    individual = I.Individual(config)
    individual.makeFromFile(model_path)
    fitness = F.Fitness(config)

    raw_predictions = []
    for seq, actual in zip(fitness.sequences, fitness.fitness_values):
        fitness.resetIndividual(individual)
        _, prediction = fitness.eval(seq, actual, individual, True)
        raw_predictions.append(prediction)

    if len(set(raw_predictions)) <= 1:
        slope, intercept = 1.0, 0.0
    else:
        slope, intercept = np.polyfit(raw_predictions, fitness.fitness_values, 1)
    return individual, fitness, float(slope), float(intercept)


def rank_candidates(
    model_path,
    config,
    backbone,
    max_swaps,
    exact_swaps=None,
    donor_species=None,
    reuse_donors=True,
    top=100,
):
    individual, fitness, slope, intercept = _load_model_and_alignment(
        model_path, config
    )

    scored = []
    swap_counts = [exact_swaps] if exact_swaps else range(1, max_swaps + 1)
    for swap_count in swap_counts:
        for tokens, swaps in CC.generate_candidates(
            backbone,
            swap_count,
            donor_species=donor_species,
            reuse_donors=reuse_donors,
        ):
            sequence = SC.join_tokens(tokens, config)
            fitness.resetIndividual(individual)
            _, raw_score = fitness.eval(sequence, 0.0, individual, True)
            predicted = raw_score * slope + intercept
            scored.append(
                {
                    "backbone": backbone,
                    "swap_count": swap_count,
                    "sequence": sequence,
                    "predicted_fitness": predicted,
                    "raw_score": raw_score,
                    "swaps": ";".join(
                        "S{:02d}<={}_{:02d}".format(*swap) for swap in swaps
                    ),
                }
            )

    scored.sort(key=lambda row: row["predicted_fitness"], reverse=True)
    return scored[:top], len(scored), slope, intercept


def main():
    parser = argparse.ArgumentParser(
        description="Rank any-to-any chimera candidates with a POET model."
    )
    parser.add_argument("model", help="Trained POET model CSV.")
    parser.add_argument("-c", "--config", default="configs/chimera_multihance.ini")
    parser.add_argument("--backbone", choices=["B2", "B3"], default="B3")
    parser.add_argument("--max-swaps", type=int, default=1)
    parser.add_argument(
        "--exact-swaps",
        type=int,
        default=None,
        help="Rank only candidates with this exact swap count.",
    )
    parser.add_argument(
        "--donor-species",
        choices=["B2", "B3"],
        default=None,
        help="Defaults to the opposite species from the backbone.",
    )
    parser.add_argument(
        "--no-reuse-donors",
        action="store_true",
        help="Do not reuse donor regions in multi-swap candidates.",
    )
    parser.add_argument("--top", type=int, default=100)
    parser.add_argument("-o", "--output", default=None)
    args = parser.parse_args()

    config = eletility.ConfigParser().read(args.config)
    top_rows, total, slope, intercept = rank_candidates(
        args.model,
        config,
        args.backbone,
        args.max_swaps,
        args.exact_swaps,
        args.donor_species,
        not args.no_reuse_donors,
        args.top,
    )

    df = pd.DataFrame(top_rows)
    print(
        "Scored {} candidates. Alignment: predicted = raw * {:.6f} + {:.6f}".format(
            total, slope, intercept
        )
    )
    print(df.to_string(index=False))
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
        print("Saved {}".format(args.output))


if __name__ == "__main__":
    main()

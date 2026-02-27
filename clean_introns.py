"""
Clean intron (inactive) nodes from a POET model CSV.

Removes rows where status != 1, deduplicates, and rounds weights to 2 decimals.

Usage:
    python clean_introns.py <model_directory> <model_number>
    e.g. python clean_introns.py output 0
    → reads output/model_0.csv → writes output/cleaned/model_0.csv
"""

import sys
import pandas as pd
from pathlib import Path


def clean_model(model_path: str, output_path: str | None = None):
    """Clean a single model CSV file by removing inactive rules.

    Parameters
    ----------
    model_path : str
        Path to the model CSV.
    output_path : str, optional
        Where to save cleaned model.  Defaults to ``<dir>/cleaned/<filename>``.
    """
    model = pd.read_csv(model_path, index_col=0)
    cleaned = model[model["status"] == 1].drop_duplicates()
    cleaned = cleaned.round(2)

    if output_path is None:
        parent = Path(model_path).parent
        cleaned_dir = parent / "cleaned"
        cleaned_dir.mkdir(parents=True, exist_ok=True)
        output_path = cleaned_dir / Path(model_path).name

    cleaned.to_csv(output_path)
    print(f"Cleaned model saved to {output_path}  ({len(cleaned)} active rules)")
    return cleaned


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python clean_introns.py <model_directory> <model_number>")
        sys.exit(1)

    model_dir = sys.argv[1]
    model_num = sys.argv[2]
    filepath = f"{model_dir}/model_{model_num}.csv"
    print(f"Cleaning: {filepath}")
    clean_model(filepath)

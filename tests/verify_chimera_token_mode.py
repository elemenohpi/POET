"""Verify token-mode POET on the chimera MultiHance dataset."""

import os
import random as rand
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "chimera"))

import chimera_codec as CC
import eletility
import optimizer
import pop as population
import sequence_codec as SC
from build_chimera_dataset import build_dataset


def _validate_training_data(path, config):
    df = pd.read_csv(path)
    if df.empty:
        raise AssertionError("chimera training data is empty")
    for sequence in df["sequence"]:
        tokens = SC.split_tokens(sequence, config)
        if len(tokens) != CC.DEFAULT_REGION_COUNT:
            raise AssertionError("sequence does not have 12 slot tokens: " + sequence)
        for idx, token in enumerate(tokens, start=1):
            info = CC.parse_token(token)
            if info["target_slot"] != idx:
                raise AssertionError("token slot/order mismatch: " + sequence)


def main(config_path="configs/chimera_multihance.ini", gens=10, seed=333):
    config = eletility.ConfigParser().read(config_path)

    build_dataset(
        Path("chimera_init_data/ChimeraICPdata.xlsx"),
        Path(config["learn_data"]),
        Path(config["alphabet_data"]),
        Path("chimera_init_data/chimera_token_manifest.csv"),
        Path("data/chimera_multihance_raw_replicates.csv"),
        Path("chimera_init_data/chimera_skipped_groups.csv"),
    )
    _validate_training_data(config["learn_data"], config)

    config["runs"] = str(gens)
    config["seed"] = str(seed)
    config["workers"] = "1"
    config["output_evo"] = "output/chimera_verify/evo.csv"
    config["output_model"] = "output/chimera_verify/model.csv"
    Path("output/chimera_verify").mkdir(parents=True, exist_ok=True)

    rand.seed(seed)
    pop = population.Population(config)
    opt = optimizer.Optimizer(config, pop)
    opt.optimize()

    model = pd.read_csv(config["output_model"])
    evo_lines = [
        line.strip()
        for line in Path(config["output_evo"]).read_text().splitlines()
        if line.strip()
    ]
    if model.empty:
        raise AssertionError("verification model is empty")
    if len(evo_lines) - 1 < gens:
        raise AssertionError("evolution log has fewer rows than requested gens")
    final_best = evo_lines[-1].split(",")[1]

    print()
    print("=" * 60)
    print("Chimera token-mode verification")
    print("  Generations : {}".format(gens))
    print("  Training rows: {}".format(len(pd.read_csv(config["learn_data"]))))
    print("  Model rules : {}".format(len(model)))
    print("  Final best  : {}".format(final_best))
    print("=" * 60)
    print("PASS - token-mode chimera run completed.")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/chimera_multihance.ini"
    g = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    s = int(sys.argv[3]) if len(sys.argv) > 3 else 333
    main(cfg, g, s)

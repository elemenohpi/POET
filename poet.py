import time
import random as rand

import archivist
import pop as population
import optimizer
from cli import build_parser, manage_input


def main():
    print("\n\n######################################################\n")
    print("POET V2.0b \n")
    print("######################################################\n")
    print("Configuring the application...\n")

    parser = build_parser()
    args = parser.parse_args()

    config = manage_input(args)

    rand.seed(int(config["seed"]))

    arch = archivist.Archivist(config)
    arch.setup()

    pop = population.Population(config)
    opt = optimizer.Optimizer(config, pop)

    start = time.time()
    opt.optimize()
    elapsed = time.time() - start

    print(
        "\nPOET evolution completed in {:.2f} seconds ({:.2f} minutes).".format(
            elapsed, elapsed / 60
        )
    )

    # ── Optional CGP post-processing layer ──
    use_cgp = config.get("use_cgp", "False").strip().lower() == "true"
    if use_cgp:
        try:
            from cgp_runner import run_cgp_on_model

            cgp_start = time.time()
            run_cgp_on_model(config)
            cgp_elapsed = time.time() - cgp_start
            print(
                "CGP post-processing completed in {:.2f} seconds ({:.2f} minutes).".format(
                    cgp_elapsed, cgp_elapsed / 60
                )
            )
        except ImportError as e:
            print(f"\nWARNING: CGP post-processing skipped — {e}")
            print(
                "To enable CGP, place the cgp/ package from "
                "https://github.com/MarkKocherovsky/cgp_crossover (branch With-SMAC) "
                "in the POET project root."
            )

    total = time.time() - start
    print(
        "\nTotal experiment completed in {:.2f} seconds ({:.2f} minutes).".format(
            total, total / 60
        )
    )


if __name__ == "__main__":
    main()

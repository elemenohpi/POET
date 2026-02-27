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
        "\nExperiment completed in {:.2f} seconds ({:.2f} minutes).".format(
            elapsed, elapsed / 60
        )
    )


if __name__ == "__main__":
    main()

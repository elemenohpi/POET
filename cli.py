import argparse

import eletility
from commands import measure_dataset_against_models, compare_models
from hpcc import hpcc


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Finds a model to predict fitness value of any given protein sequence. "
            "Fitness can be manually defined to any protein characteristic but our "
            "main goal is to predict CEST ability of proteins"
        )
    )
    parser.add_argument(
        "-config", help="Takes the config file to configure the application"
    )
    parser.add_argument("-o", help="Output file name. Include extension")
    parser.add_argument("-r", help="Number of GP iterations to find a model.")
    parser.add_argument(
        "-hpcc",
        help="Runs the experiment on the HPCC servers through slurm jobs",
        action="store_true",
    )
    parser.add_argument("-seed", help="The random seed")
    parser.add_argument(
        "-mo", help="Specifies a path to the model output. Include extension"
    )
    parser.add_argument(
        "-predict",
        help="Predict potential protein sequences. Use jointly with -seqsize and -iter",
        action="store_true",
    )
    parser.add_argument(
        "-f", help="Gets path to a model as its input and returns the fitness of it"
    )
    parser.add_argument(
        "-md",
        help="Expects a model directory. Returns a table of predictions, actual values and RMSE",
    )
    parser.add_argument(
        "-c", nargs="*", help="Compares the fitness of all given models"
    )
    parser.add_argument(
        "-al", nargs="*", help="Computes the average length of all given models"
    )
    parser.add_argument("-learn", help="Path to the learn data (format: csv)")
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=0,
        help="Number of parallel worker processes for fitness evaluation. "
        "0 = auto-detect CPU count (default), 1 = sequential.",
    )
    parser.add_argument(
        "--experiment",
        nargs="+",
        metavar="CONFIG",
        help="Compare two or more config files across multiple seeds. "
        "Use with --replicates and --gens.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=10,
        help="Number of replicates (seeds) for experiment mode (default: 10).",
    )
    parser.add_argument(
        "--gens",
        type=int,
        default=100,
        help="Number of generations per replicate in experiment mode (default: 100).",
    )
    return parser


def manage_input(args):
    configparser = eletility.ConfigParser()
    config = configparser.read(args.config if args.config else "config.ini")

    if args.seed:
        config["seed"] = int(args.seed)
    if args.r:
        config["runs"] = int(args.r)
    if args.o:
        config["output_evo"] = args.o
    if args.mo:
        config["output_model"] = args.mo
    if args.learn:
        config["learn_data"] = args.learn

    config["workers"] = str(args.workers)

    if args.experiment:
        from experiment import run_experiment

        configs = [configparser.read(p) for p in args.experiment]
        run_experiment(configs, args.replicates, args.gens, args.workers)
        exit()

    if args.f:
        raise NotImplementedError("-f is not yet tested")

    if args.md:
        measure_dataset_against_models(config, args.md)
        exit()

    if args.c:
        compare_models(config, args.c)
        exit()

    if args.hpcc:
        hpcc()
        exit()

    if args.al:
        raise NotImplementedError("-al is not yet tested")

    if args.predict:
        print(
            "Predicting proteins:\nPopulation pool is set to be 1000000\n================================"
        )
        count = int(input("Enter prediction count: "))
        seq_size = int(input("Enter protein sequences size: "))
        iterations = int(
            input(
                "Enter number of evolutionary iterations (Larger values results in more confident and "
                "yet similar predictions.\n Lower values makes room for novelty but the prediction might"
                " not be as accurate): "
            )
        )
        model = input("Enter the path to the predictor model(s) directory: ")
        config_path = input(
            "Enter the path to the config file used to run the experiment: "
        )
        config = configparser.read(config_path)
        import predictor as P

        pred = P.Predictor(count, seq_size, iterations, config, model)
        pred.predict()
        exit()

    return config

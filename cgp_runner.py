import argparse
import os
import eletility
import pandas as pd
import numpy as np
import random as rand

from pathlib import Path
from cgp import CartesianGP
from cgp.cgp_operators import add, sub, mul, div

def add_arguments_to_parser(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--config",
        default="config.ini",
        help="Path to config.ini file (default: config.ini)",
    )
    return parser


def manage_input(args) -> dict:
    configparser = eletility.ConfigParser()
    # args.config is a PATH string here
    config = configparser.read(args.config)
    return config


def rules_to_constants(model:pd.DataFrame, peptides:pd.DataFrame) -> tuple[list, list]:
    """
        The model has to have, at minimum, a "pattern" and "weight" column
        The peptides have to have a "Sequence" and a "Score" column
    """
    input_data = []
    target_output = []
    for i, peptide in enumerate(peptides.itertuples()):
        target_output.append(float(peptide.Score))
        sequence = peptide.Sequence
        datum = []
        for j, motif in enumerate(model.itertuples()):
            if motif.pattern in sequence:
                datum.append(motif.weight)
            else:
                datum.append(0.0)
        input_data.append(datum)
    return input_data, target_output

def set_up_cgp(x, y, seed):
    trial_number = 0
    max_generations = 1000
    model_size = 32
    xover_type = None
    xover_rate = 0.0
    max_parents = 1
    max_children = 4
    mutation_type = "full"
    mutation_rate = 1.0
    selection_type = "paretoelite"
    fitness_function = "correlation"
    test_problem_key = "poet"
    n_points = 1
    tournament_size = 4
    n_elites = 1
    step_size = 20
    asex = True

    model_parameters = {
        'max_size': model_size,
        'inputs': x.shape[-1],
        'outputs': 1,
        'arity': 2,
        'constants': np.array([])
    }
    function_bank = {'add': add, 'sub': sub, 'mul': mul, 'div': div}

    if asex or max_parents < max_children:
        mutation_breeding = True
    else:
        mutation_breeding = False
    CHECKPOINT_PATH = os.path.join(os.environ.get("SCRATCH", "/tmp"), "ckpt")
    CHECKPOINT_FILE = f"{CHECKPOINT_PATH}/{test_problem_key}_trial_{trial_number}_ckpt.pkl"

    evolution_module = CartesianGP(
            parents=max_parents,
            children=max_children,
            max_generations=max_generations,
            mutation=mutation_type,
            selection=selection_type,
            xover=xover_type,
            fixed_length=True,
            fitness_function=fitness_function,
            model_parameters=model_parameters,
            n_points=n_points,
            n_elites=n_elites,
            tournament_size=tournament_size,
            function_bank=function_bank,
            mutation_breeding=mutation_breeding,
            checkpoint_filename=CHECKPOINT_FILE,
            one_dimensional_xover=False,
            seed=seed,
            tuning=False
        )
    return evolution_module
    

def main():
    parser = argparse.ArgumentParser()
    add_arguments_to_parser(parser)
    args = parser.parse_args()

    config = manage_input(args)

    rand.seed(config["seed"])
    model_path = config["prediction_model"]

    model = pd.read_csv(model_path)

    data_path = config["learn_data"]
    unseen_data_path = config["unseen_data"]

    data = pd.read_csv(data_path)
    unseen_data = pd.read_csv(data_path)

    input_data, target_fitness = rules_to_constants(model, data)
    input_data = np.array(input_data)
    target_fitness = np.array(target_fitness)

    input_test_data, test_fitness = rules_to_constants(model, unseen_data)
    input_test_data = np.array(input_test_data)
    test_fitness = np.array(test_fitness)

    target_fitness = np.array(target_fitness)
    evolver = set_up_cgp(input_data, target_fitness, int(config["seed"]))

    run_path = f'output'
    Path(run_path).mkdir(parents=True, exist_ok=True)
    
    best_model, best_test_model = evolver.fit(input_data, input_test_data, target_fitness, test_fitness, step_size=10)
    
    evolver.save_metrics(run_path)
    print('---')
    best_test_model.print_model()
    df = pd.DataFrame(best_test_model.model)
    df.to_csv(f'{run_path}/best_model.csv', index=True)
    print(f'Complexity: {best_test_model.count_active_nodes()}')

if __name__ == "__main__":
    main()


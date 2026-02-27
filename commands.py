import os

import individual as I
import fitness as F


def measure_dataset_against_models(config, path):
    all_files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    csv_files = [f for f in all_files if f.endswith(".csv")]

    models = []
    for file in csv_files:
        model = I.Individual(config)
        model.makeFromFile(os.path.join(path, file))
        models.append(model)

    fitness = F.Fitness(config)
    seq_fitness_tuples, individuals_evaluations = fitness.model_vs_dataset(
        config, models
    )

    header = "sequence,fitness," + ",".join(
        f.split(".")[0] + " Prediction" for f in csv_files
    )
    rows = [header]
    for seq_index, seq_fitness in enumerate(seq_fitness_tuples):
        cols = [str(seq_fitness[0]), str(seq_fitness[1])]
        cols += [
            str(round(individual_evals[seq_index], 2))
            for individual_evals in individuals_evaluations
        ]
        rows.append(",".join(cols))

    with open("table.csv", "w") as text_file:
        text_file.write("\n".join(rows) + "\n")


def compare_models(config, paths):
    paths = paths[0]
    files = [f for f in os.listdir(paths) if os.path.isfile(os.path.join(paths, f))]

    model = I.Individual(config)
    f = F.Fitness(config)
    avg = 0
    best = 100000
    best_model = ""
    for file in files:
        if file.split(".")[-1] != "csv":
            continue
        file = os.path.join(paths, file)
        model.makeFromFile(file)
        fitness, test = f.measureTotal(model)
        avg += fitness
        print(
            "Pro-Predictor: Fitness ({}) of {}: {} test: {}".format(
                config["fitness_alg"], file, fitness, test
            )
        )
        if fitness < best:
            best_model = file
            best = fitness
    avg /= len(files)
    print(
        "Pro-Predictor: Best model: {} with {}: {} Average {}: {}".format(
            best_model, config["fitness_alg"], best, config["fitness_alg"], avg
        )
    )

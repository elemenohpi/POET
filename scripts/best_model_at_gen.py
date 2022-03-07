import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt


def find_best(path, n):
    path = path[0]
    files = [f for f in listdir(path) if isfile(join(path, f))]
    dataframes = []
    for file in files:
        # print(file)

        filetokens = file.split(".")
        if filetokens[len(filetokens) - 1] != "csv":
            continue
        file = path + file
        df = pd.read_csv(file, index_col=0)
        dataframes.append(df)

    values = []

    for df in dataframes:
        values.append(df["best fitness"].iloc[n])

    max_value = max(values)
    max_index = values.index(max_value)

    training_fitness = dataframes[max_index]["best fitness"].iloc[n]
    test_fitness = dataframes[max_index][" best test"].iloc[n]
    rule_count = dataframes[max_index][" best rule count"].iloc[n]
    unused_rule_count = dataframes[max_index][" best unused rule count"].iloc[n]
    used_rule_count = rule_count - unused_rule_count

    print("Best Model Performance:\nTraining Fitness: {}\nTest Fitness: {}\nRule Count: {}\n"
          "Unused Rule Count: {}\nUsed Rule Count: {}".format(training_fitness, test_fitness, rule_count,
                                                              unused_rule_count, used_rule_count))


def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="returns the best model performance at gen N given set of evo files")
    parser.add_argument("-n", help="N")

    args = parser.parse_args()

    if args.f is not None and args.n is not None:
        find_best(args.f, int(args.n))
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

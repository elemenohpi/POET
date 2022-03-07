import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np


def plotIt(path, title):
    q25s = []
    q75s = []
    medians = []
    gens = list(range(0, 5000))
    path = path[0]
    files = [f for f in listdir(path) if isfile(join(path, f))]
    dataframes = []
    for file in files:
        filetokens = file.split(".")
        if filetokens[len(filetokens) - 1] != "csv":
            continue
        file = path + file
        df = pd.read_csv(file, index_col=0)
        dataframes.append(df)

    for i in range(5000):
        fittest_at_time = []
        for df in dataframes:
            fittest_at_time.append(df["best fitness"][i])

        fittest_df = pd.DataFrame(fittest_at_time)
        quantiles = fittest_df.quantile([0.25, 0.75])
        medians.append(fittest_df.median().values[0])
        q25s.append(quantiles[0].iloc[0])
        q75s.append(quantiles[0].iloc[1])

    plt.style.use('fast')

    # plot
    fig, ax = plt.subplots()

    plt.xlabel("Generations")
    plt.ylabel("Fitness (RMSE)")
    plt.title(title)

    ax.fill_between(gens, q25s, q75s, alpha=.5, linewidth=0)
    ax.plot(gens, medians, linewidth=2.0)

    # plt.show()
    plt.savefig("./gen_evo/fig.png", dpi=300)

def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="provide directory to plot median and quartiles. only considers files "
                                              "that start with 'new' because of the output of the duplication script")
    parser.add_argument("-t", help="plot title")

    args = parser.parse_args()

    if args.f is not None and args.t is not None:
        plotIt(args.f, args.t)
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

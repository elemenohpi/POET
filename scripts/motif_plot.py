import argparse
import eletility
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

F = eletility.Files()


def create_motif_frequency_map(path):
    files = [f for f in listdir(path) if isfile(join(path, f))]

    motif_frequency_map = {}

    for file in files:
        file = path + file

        with open(file, encoding="utf-8-sig") as my_file:
            lines = my_file.readlines()

        model_motifs = []

        for line in lines:
            line = line.strip()
            tokens = line.split(",")
            if not tokens[0].isnumeric():
                continue

            if tokens[3] == "0":
                continue

            motif = tokens[1]

            if motif not in model_motifs:
                model_motifs.append(tokens[1])

        map_keys = motif_frequency_map.keys()

        for motif in model_motifs:
            if motif in map_keys:
                motif_frequency_map[motif] += 1
            else:
                motif_frequency_map[motif] = 1

    sorted_map = sort_motif_frequency_map(motif_frequency_map)

    motif_list = []
    frequency_list = []

    for index, key in enumerate(reversed(sorted_map.keys())):
        if sorted_map[key] > 5:
            motif_list.append(key)
            frequency_list.append(sorted_map[key])

    plot(motif_list, frequency_list)


def plot(motifs, freq):
    plt.style.use('fast')

    # plot
    fig, ax = plt.subplots()

    plt.xlabel("Discovered Motifs")
    plt.ylabel("Frequency of Appearance in POET Models")
    plt.title("Most Frequently Discovered Motifs Among the Best Models of 50 Repeats of Epoch 8")
    plt.xticks(rotation=90)

    ax.bar(motifs, freq, width=1, edgecolor="white", linewidth=0.7)

    plt.show()


def sort_motif_frequency_map(map):
    sorted_map = {k: v for k, v in sorted(map.items(), key=lambda item: item[1])}
    return sorted_map


def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="prints out motif frequency")

    args = parser.parse_args()

    if args.f is not None:
        create_motif_frequency_map(args.f[0])
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

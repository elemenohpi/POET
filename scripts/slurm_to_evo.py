import argparse
import eletility

from os import listdir
from os.path import isfile, join

F = eletility.Files()


def fix(path):
    path = path[0]
    files = [f for f in listdir(path) if isfile(join(path, f))]
    file_counter = 1
    for file in files:
        if file[:5] != "slurm":
            continue
        file = path + file
        name = file.split("\\")
        new_name = str(file_counter) + "_evo.csv"
        file_counter += 1
        directory = name[:len(name) - 1]
        path = ""
        for part in directory:
            path += part + "\\"

        new_name = "gen_evo/" + new_name

        with open(file, encoding="utf-8-sig") as myfile:
            lines = myfile.readlines()

        F.truncate(new_name)
        F.writeLine(new_name, "best fitness, best test, best rule count, best unused rule count, average fitness, "
                              "average test, average rule count, average unused rule count")

        newlines = []
        prev_gen = 100000000
        for line in reversed(lines):
            tokens = line.strip().split(" ")
            gen = tokens[0].split(":")[0]
            if not gen.isnumeric():
                continue
            if int(gen) >= prev_gen:
                continue
            prev_gen = int(gen)
            constructed_line = gen + "," + tokens[2] + "," + tokens[4] + "," + tokens[6] + "," + tokens[8] + "," \
                               + tokens[11] + "," + tokens[13] + "," + tokens[15] + "," + tokens[17]
            newlines.append(constructed_line)
            if gen == "0":
                break

        for line in reversed(newlines):
            F.writeLine(new_name, line)

def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="provide directory to slurm files to remove duplicates and order "
                                              "problems and create evo files in the gen_evo folder. run from main "
                                              "folder")

    args = parser.parse_args()

    if args.f is not None:
        fix(args.f)
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

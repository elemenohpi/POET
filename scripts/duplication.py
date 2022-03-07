import argparse
import eletility

from os import listdir
from os.path import isfile, join

F = eletility.Files()


def fix(path):
    path = path[0]
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        file = path + file
        name = file.split("\\")
        new_name = "new_" + name[len(name) - 1]
        directory = name[:len(name) - 1]
        path = ""
        for part in directory:
            path += part + "\\"
        new_name = path + new_name

        with open(file, encoding="utf-8-sig") as myfile:
            lines = myfile.readlines()

        new_lines = [""]
        index = 0
        for line in lines:
            if line.strip() != new_lines[index]:
                # print("'", line.strip(), "'", "\n", "'",new_lines[index].strip(), "'", "\n\n")
                new_lines.append(line.strip())
                index += 1
        new_lines = new_lines[1:]
        for line in new_lines:
            F.writeLine(new_name, line)

def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="provide directory to remove duplicates")

    args = parser.parse_args()

    if args.f is not None:
        fix(args.f)
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

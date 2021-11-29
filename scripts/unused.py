import argparse
import eletility
from os import listdir
from os.path import isfile, join

F = eletility.Files()

motifs = {}

def add(motif):
    pass

def fix(path):
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        new_name = path + "processed_" + file

        with open(path+file, encoding="utf-8-sig") as my_file:
            lines = my_file.readlines()

        new_lines = []
        for line in lines:
            data = line.split(",")
            if data[3].strip() != "0" and data[1] != "pattern":
                # new_lines.append(line.strip())
                try:
                    motifs[data[1]] += 1
                except:
                    motifs[data[1]] = 1

    frequency = {}
    
    keys = motifs.keys()
    for key in keys:
        count = motifs[key]
        try:
            frequency[count] += 1
        except:
            frequency[count] = 1
    
    # for key in motifs.keys():
    #     print(key + "," + repr(motifs[key]))
    
    
    
    print(frequency)
    exit()

#        for line in new_lines:
#            F.writeLine(new_name, line)


def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="fixes the duplication problem in all given files")

    args = parser.parse_args()

    if args.f is not None:
        fix(args.f[0])
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

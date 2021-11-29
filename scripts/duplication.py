import argparse
import eletility

F = eletility.Files()


def fix(files):
    for file in files:
        name = file.split("\\")
        new_name = "new_" + name[len(name) - 1]
        directory = name[:len(name) - 1]
        path = ""
        for part in directory:
            path += part + "\\"
        new_name = path + new_name

        lines = []
        with open(file, encoding="utf-8-sig") as myfile:
            lines = myfile.readlines()
        new_lines = [""]
        index = 0
        for line in lines:
            if line.strip() is not new_lines[index].strip():
                new_lines.append(line.strip())
                index += 1
        new_lines = new_lines[1:]
        for line in new_lines:
            F.writeLine(new_name, line)


def main():
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('-f', nargs='*', help="fixes the duplication problem in all given files")

    args = parser.parse_args()

    if args.f is not None:
        fix(args.f)
    else:
        raise "arguments not provided"


if __name__ == "__main__":
    main()

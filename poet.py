import argparse
# from os import listdir
# from os.path import isfile, join

# import pandas as pd
import pop as population
import optimizer as optimizer
import predictor as P
import individual as I
# import math
import fitness as F
import os.path
import random as rand
from subprocess import call
import archivist as archivist
# import pandas as pd
import eletility
from datetime import date





def add_arguments_to_parser(parser):
    parser.add_argument("-config", help="Takes the config file to configure the application")
    parser.add_argument('-o', help="Output file name. Include extension")
    parser.add_argument('-r', help='Number of GP iterations to find a model.')
    parser.add_argument('-hpcc', help="Runs the experiment on the HPCC servers through slurm jobs", action="store_true")
    parser.add_argument('-seed', help='The random seed')
    parser.add_argument("-mo", help="Specifies a path to the model output. Include extension")
    parser.add_argument('-predict',
                        help='Number of potential protein sequences you want the program to predict. Must be '
                             'jointly used with -seqsize and -iter', action="store_true")
    parser.add_argument('-f',
                        help="Gets path to a model as it's input and returns the fitness of it")
    parser.add_argument('-md',
                        help="Expects a model to be given. Returns a table of predictions, actual values and RMSE for that model")
    parser.add_argument('-c', nargs='*', help="Compares the fitness of all given models")
    parser.add_argument('-al', nargs='*', help="Computes the average length of all given models")
    #
    parser.add_argument('-learn', help='Path to the learn data (format: csv)')
    # parser.add_argument('-translation',
    #                     help='Path to the translation table (format: csv, default: )')
    # parser.add_argument('-pop',
    #                     help='Path to the initial population files. If not specified, this application uses random initial population as default (format: csv)')
    # parser.add_argument('-model',
    #                     help='Path to a generated model to determine a given protein\'s fitness. You will need to use -seq after this command')
    # parser.add_argument('-seq',
    #                     help='A sequence to be tested using an already specified model. This command only runs if it\'s used jointly with the -model command')
    # parser.add_argument('-seqsize', help='Size of the protein sequences for prediction')
    # parser.add_argument('-iter', help='Number of iterations to predict/find potential proteins')
    #
    # # parser.add_argument('-archive', nargs='*',
    # #                     help='Setups the default output directories if necessary and archives existing files/results')

    return parser


def manage_input(args):
    configparser = eletility.ConfigParser()
    if args.config:
        config = configparser.read(args.config)
    else:
        config = configparser.read("config.ini")

    if args.seed:
        config["seed"] = int(args.seed)

    # Runs = Number of evolutionary generations
    if args.r:
        # number of evolutionary generations
        config["runs"] = int(args.r)

    if args.o:
        config["output_evo"] = args.o

    if args.mo:
        config["output_model"] = args.mo

    if args.learn:
        config["learn_data"] = args.learn

    # get the fitness of a given model
    if args.f:
        raise "Testing needed"
        modelFitness(args.f)

    if args.md:
        path = args.md
        measure_dataset_against_models(config, path)
        exit()

    # compare a bunch of models
    if args.c:
        # ToDo:: double check the correlations vs. logs. Slightly different values are found in the results. I suspect it's due to either how the models are being saved or how the correlation function is returning the r value. could be a non-issue and caused by the cross validation.
        compare_models(config, args.c)
        exit()

    # run on hpcc
    if args.hpcc:
        hpcc()
        exit(1)

    if args.al:
        raise "Testing needed"
        averageLength(args.al)
        exit()

    if args.predict:
        print("Predicting proteins:\n================================".format(args.predict))
        count = int(input("Enter prediction count: "))
        seq_size = int(input("Enter protein sequences size: "))
        iterations = int(input("Enter number of evolutionary iterations (Larger values results in more confident and "
                               "yet similar predictions.\n Lower values makes room for novelty but the prediction might"
                               " not be as accurate): "))
        model = input("Enter the path to the predictor model(s): ")
        predictor = P.Predictor(count, seq_size + 1, iterations, config, model)
        predictor.predict()
        exit()

    return config


def main():
    print("\n\n######################################################\n")
    print("POET V2.0b \n")
    print("######################################################\n")

    print("Configuring the application...\n")

    # Argument descriptions
    parser = argparse.ArgumentParser(
        description='Finds a model to predict fitness value of any given protein sequence. Fitness can be manually defined to any protein characteristic but our main goal is to predict CEST ability of proteins')

    parser = add_arguments_to_parser(parser)
    args = parser.parse_args()

    config = manage_input(args)
    arch = archivist.Archivist(config)

    # Setting up the random seed
    rand.seed(config["seed"])

    # Setups the architecture of the project ToDo:: Might be a good idea to improve this or remove it altogether
    arch.setup()

    pop = population.Population(config)
    # for i in pop.pop: # display all individuals
    #     print(i.print())
    
    opt = optimizer.Optimizer(config, pop)
    opt.optimize() # Optimize the population


    # elif args.pop == None and args.seq == None and args.model == None:

    # elif args.seq != None and args.model != None:
    #     proteinSeq = args.seq
    #     modelPath = args.model
    #     model = I.Individual()
    #     model.makeFromFile(modelPath)
    #     fObj = F.Fitness()
    #     # prediction = fObj.predict(proteinSeq, model)
    #     error, prediction = fObj.eval(proteinSeq, 12.7, model, True)
    #     print(prediction)
    # else:
    #     raise ("Invalid arguements.")
    #     pass
    # return


def measure_dataset_against_models(config, path):
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    models = []
    for file in files:
        if ".csv" not in file:
            continue
        model = I.Individual(config)
        model.makeFromFile(os.path.join(path, file))
        models.append(model)
    f = F.Fitness(config)
    seq_fitness_tuples, individuals_evaluations = f.model_vs_dataset(config, models)
    csv = "sequence,fitness,"
    for file in files:
        if ".csv" not in file:
            continue
        filename = file.split(".")[0]
        csv += filename + " Prediction,"
    csv = csv[:-1] + "\n"
    for seq_index, seq_fitness in enumerate(seq_fitness_tuples):
        row = ""
        row += str(seq_fitness[0]) + "," + str(seq_fitness[1]) + ","
        for individual_evals in individuals_evaluations:
            row += str(round(individual_evals[seq_index], 2)) + ","
        row = row[:-1] + "\n"
        csv += row
    text_file = open("table.csv", "w")
    n = text_file.write(csv)
    text_file.close()
    pass



def modelFitness(path):
    model = I.Individual()
    model.makeFromFile(path)
    f = F.Fitness()
    fitness, test = f.measureTotal(model)
    print("Pro-Predictor: Fitness (RMSE) of {}: {}".format(path, fitness))
    exit(1)


def compare_models(config, paths):
    paths = paths[0]
    files = [f for f in os.listdir(paths) if os.path.isfile(os.path.join(paths, f))]

    model = I.Individual(config)
    avg = 0
    best = 100000
    bestModel = ''
    for file in files:
        filetokens = file.split(".")
        if filetokens[len(filetokens) - 1] != "csv":
            continue
        file = os.path.join(paths, file)
        # file = paths + file
        model.makeFromFile(file)
        f = F.Fitness(config)
        fitness, test = f.measureTotal(model)
        avg += fitness
        print("Pro-Predictor: Fitness ({}) of {}: {} test: {}".format(config["fitness_alg"], file, fitness, test))
        if fitness < best:
            bestModel = file
            best = fitness
    avg /= len(paths)
    print("Pro-Predictor: Best model: {} with {}: {} Average {}: {}".format(bestModel, config["fitness_alg"], best,
                                                                            config["fitness_alg"], avg))
    exit(1)


def averageLength(paths):
    model = I.Individual()
    for index, path in enumerate(paths):
        model.makeFromFile(path)
        f = F.Fitness()
        fitness = f.measureTotal(model)
        sumL = 0
        counter = 0
        for rule in model.rules:
            if rule.status:
                try:
                    if math.isnan(rule.pattern):
                        continue
                except:
                    pass
                sumL += len(rule.pattern)
                counter += 1
        if counter == 0:
            print("Pro-Predictor: {}: Old Model. Old Models are Not Supported".format(path))
            continue
        avg = sumL / counter
        print("Pro-Predictor: {}: Fitness (RMSE): {} Average Rule Length: {}".format(path, fitness, avg))
    exit(1)


def hpcc():
    today = date.today()
    day = today.strftime("%b-%d-%Y")

    print("POET is running on the HPCC experiment mode. Please enter the following information\n"
          "======================================")
    title = input("Experiment Title: ")
    hours = int(input("Requested Hours: "))
    reps = int(input("Number of Experiment Repeats: "))
    runs = int(input("Evolutionary Generations: "))
    seed = int(input("Starting Random Seed: "))
    config = input("Experiment configuration file: ")
    confirmation_text = "#======================================\n#Experiment Title: {}\n#Hours: {}\n" \
                        "#Repeats: {}\n#Generations: {}\n#Seed: {}\n#Config: {}".format(title, hours, reps, runs, seed,
                                                                                        config)

    print(confirmation_text)

    confirm = input("To confirm the above settings, enter YES: ")
    if confirm.lower() != "yes":
        print("Aborting!")
        exit()

    title = day + "-" + title

    if os.path.exists("output/{}".format(title)):
        raise "An experiment with the same title exists"
    else:
        directory = os.path.join("output", title)
        subs_directory = os.path.join(directory, "subs")
        slurms_directory = os.path.join(directory, "slurms")
        logs_directory = os.path.join(directory, "logs")
        models_directory = os.path.join(directory, "models")
        errors_directory = os.path.join(directory, "errors")
        os.makedirs(directory)
        os.makedirs(subs_directory)
        os.makedirs(logs_directory)
        os.makedirs(models_directory)
        os.makedirs(errors_directory)
        os.makedirs(slurms_directory)

    content = confirmation_text
    with open(config, "r") as cfile:
        content += "\n" + cfile.read()

    file_handler = eletility.Files()
    file_handler.writeTruncate(os.path.join(directory, "config.ini"), content)

    for i in range(reps):
        filename = os.path.join(subs_directory, "{}_{}.sb".format(title, i))
        file = open(filename, "w")
        file.write("#!/bin/bash --login\n")
        file.write("\n########## SBATCH Lines for Resource Request ##########\n\n")
        file.write(
            "#SBATCH --time={}:02:00             # limit of wall clock time - how long the job will run (same as -t)\n".format(
                hours))
        file.write(
            "#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)\n")
        file.write(
            "#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)\n")
        file.write("#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)\n")
        file.write(
            "#SBATCH --mem-per-cpu=8G            # memory required per allocated CPU (or core) - amount of memory (in bytes)\n")
        file.write(
            "#SBATCH --job-name {}_{}      # you can give your job a name for easier identification (same as -J)\n".format(
                title, i))
        file.write(
            "#SBATCH --error={}/{}_{}.err      # you can give your job a name for easier identification (same as -J)\n".format(
                errors_directory, title, i))
        file.write(
            "#SBATCH --output={}/{}_{}.txt      # you can give your job a name for easier identification (same as -J)\n".format(
                slurms_directory, title, i))

        # SBATCH --error=%j.err
        file.write("\n########## Command Lines to Run ##########\n\n")
        file.write("module pandas\n")
        # file.write("module load GCC/6.4.0-2.28 OpenMPI  ### load necessary modules, e.g\n")
        file.write("cd ~/POET\n")
        log_file_path = logs_directory + "/evo_{}.csv".format(i)
        model_file_path = models_directory + "/model_{}.csv".format(i)
        file.write(
            "srun -n 1 python poet.py -config {} -r {} -o {} -mo {} -seed {}\n"
            "".format(config, runs, log_file_path, model_file_path, i + seed))
        file.write("cd batch\n")
        file.write("scontrol show job $SLURM_JOB_ID     ### write job information to output file")
        file.close()
        call(["sbatch", os.path.join(subs_directory, "{}_{}.sb".format(title, i))])


if __name__ == "__main__":
    main()

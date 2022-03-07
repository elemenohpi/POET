# POET: Protein Optimization Enhancing Tool
# Authors: 
#			Iliya Alavy - Department of Computer Science and Engineering - Michigan State University
# 		    Alexander Bricco - Department of Bioengineering -  Michigan State University
#			All Rights Reserved @ Michigan State University

import argparse
from os import listdir
from os.path import isfile, join

import pandas as pd
import pop as Population
import settings
import optimizer as Optimizer
import predictor as P
import individual as I
import math
import fitness as F
import os
import random as R
from subprocess import call
import archivist as Archivist
import pandas as pd


def main():
    print("\n\n######################################################\n")
    print("POET V2.0b \n")
    print("######################################################\n")

    print("Configuring the application...\n")

    # Arguement descriptions
    parser = argparse.ArgumentParser(
        description='Finds a model to predict fitness value of any given protein sequence. Fitness can be manually defined to any protein characteristic but our main goal is to predict CEST ability of proteins')

    parser.add_argument('-learn', default=settings.default_learn,
                        help='Path to the learn data (format: csv, default: ' + settings.default_learn + ')')
    parser.add_argument('-translation', default=settings.TT,
                        help='Path to the translation table (format: csv, default: ' + settings.TT + ')')
    parser.add_argument('-pop',
                        help='Path to the initial population files. If not specified, this application uses random initial population as default (format: csv)')
    parser.add_argument('-model',
                        help='Path to a generated model to determine a given protein\'s fitness. You will need to use -seq after this command')
    parser.add_argument('-seq',
                        help='A sequence to be tested using an already specified model. This command only runs if it\'s used jointly with the -model command')
    parser.add_argument('-predict',
                        help='Number of potential protein sequences you want the program to predict. Must be jointly used with -seqsize and -iter')
    parser.add_argument('-seqsize', help='Size of the protein sequences for prediction')
    parser.add_argument('-iter', help='Number of iterations to predict/find potential proteins')
    parser.add_argument('-hpcc', help="Number of replications you need to queue on hpcc. Uses the default config file")
    parser.add_argument('-o', help="Output file name")
    parser.add_argument('-r', default=settings.runs, help='Number of GP iterations to find a model.')
    parser.add_argument('-seed', help='The random seed')
    parser.add_argument('-f',
                        help="Gets path to a model as it's input and returns the fitness of it")  # ToDo:: Code this part.
    parser.add_argument('-c', nargs='*', help="Compares the fitness of all given models")
    parser.add_argument('-al', nargs='*', help="Computes the average length of all given models")
    parser.add_argument('-archive', nargs='*',
                        help='Setups the default output directories if necessary and archives existing files/results')
    parser.add_argument('-md',
                        help="Expects a model to be given. Returns a table of predictions, actual values and RMSE for that model")

    args = parser.parse_args()

    # Read the tables
    settings.learn_df = pd.read_csv(args.learn)
    settings.TT = pd.read_csv(args.translation)

    arch = Archivist.Archivist()

    if args.archive != None:
        arch.setup(True)
        print("Archiving completed Successfully!")
        exit(1)

    # Setting up the random seed
    if args.seed != None:
        settings.seed = int(args.seed)
    R.seed(settings.seed)

    # total number of runs
    if args.r != None:
        settings.runs = int(args.r)

    # output file name
    if args.o != None:
        settings.output_file_name = args.o

    # Should be after initializing the output file name ^^^
    arch.setup()

    # run on hpcc
    if args.hpcc != None:
        hpcc(int(args.hpcc), 24, 10000, 3333)
        exit(1)

    # get the fitness of a given model
    if args.f is not None:
        modelFitness(args.f)

    if args.md is not None:
        path = args.md
        files = [f for f in listdir(path) if isfile(join(path, f))]
        models = []
        for file in files:
            if ".csv" not in file:
                continue
            model = I.Individual()
            model.makeFromFile(join(path, file))
            models.append(model)
        f = F.Fitness()
        seq_fitness_tuples, individuals_evaluations = f.model_vs_dataset(models)
        csv = "sequence,fitness,"
        for file in files:
            if ".csv" not in file:
                continue
            filename = file.split(".")[0]
            csv += filename + " Abs Error" + "," + filename + " Prediction,"
        csv = csv[:-1] + "\n"
        for seq_index, seq_fitness in enumerate(seq_fitness_tuples):
            row = ""
            row += str(seq_fitness[0]) + "," + str(seq_fitness[1]) + ","
            for individual_evals in individuals_evaluations:
                row += str(round(individual_evals[seq_index][0], 2)) + "," + str(round(individual_evals[seq_index][1], 2)) + ","
            row = row[:-1] + "\n"
            csv += row
        text_file = open("table.csv", "w")
        n = text_file.write(csv)
        text_file.close()


        exit()

    # compare a bunch of models
    if args.c != None:
        compareModels(args.c)

    # calculate the average length of the models
    if args.al != None:
        averageLength(args.al)

    # Make a dictionary out of the translation table
    for i, row in settings.TT.iterrows():
        settings.dic[row[0]] = row[1]

    # Translate the data using the translation table in case we are not using a numeric translation table
    if not is_numeric():
        for i, row in settings.learn_df.iterrows():
            translatedSeq = ""
            origSeq = row[0]
            for j in range(len(origSeq)):
                try:
                    translatedSeq += settings.dic[origSeq[j]]
                except IndexError:
                    print("Could not find all the required data in the Translation Table. Exiting...")
                    exit(0)
            settings.learn_df.iloc[i, 0] = translatedSeq

    # ToDo:: I never checked the prediction part
    if args.predict != None:
        if args.seqsize == None or args.iter == None:
            raise ("Prediction mode needs -seqsize and -iter in order to work")
        else:
            predictor = P.Predictor(int(args.predict), int(args.seqsize) + 1, int(args.iter))
            predictor.predict()
            return
    elif args.pop == None and args.seq == None and args.model == None:
        pop = Population.Population()
        opt = Optimizer.Optimizer(pop)
        opt.optimize()
    elif args.seq != None and args.model != None:
        proteinSeq = args.seq
        modelPath = args.model
        model = I.Individual()
        model.makeFromFile(modelPath)
        fObj = F.Fitness()
        prediction = fObj.predict(proteinSeq, model)
        # error, prediction = fObj.eval(proteinSeq, 12.7, model, True)
        print(prediction)
    else:
        raise ("Invalid arguements.")
        pass
    return


def is_numeric():
    codes = settings.TT['code']
    for i in range(codes.size):
        try:
            float(codes[i])
        except ValueError:
            return False
    return True


def modelFitness(path):
    model = I.Individual()
    model.makeFromFile(path)
    f = F.Fitness()
    fitness, test = f.measureTotal(model)
    print("Pro-Predictor: Fitness (RMSE) of {}: {}".format(path, fitness))
    exit(1)


def compareModels(paths):
    paths = paths[0]
    files = [f for f in listdir(paths) if isfile(join(paths, f))]

    model = I.Individual()
    avg = 0
    best = 100000
    bestModel = ''
    for file in files:
        filetokens = file.split(".")
        if filetokens[len(filetokens) - 1] != "csv":
            continue
        file = paths + file
        model.makeFromFile(file)
        f = F.Fitness()
        fitness, test = f.measureTotal(model)
        avg += fitness
        print("Pro-Predictor: Fitness (RMSE) of {}: {} test: {}".format(file, fitness, test))
        if fitness < best:
            bestModel = file
            best = fitness
    avg /= len(paths)
    print("Pro-Predictor: Best model: {} with RMSE: {} Average RMSE: {}".format(bestModel, best, avg))
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


def hpcc(reps, hours, runs, startingSeed):
    for i in range(reps):
        filename = "subs/{}.sb".format(i)
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
            "#SBATCH --job-name POET_rep_{}      # you can give your job a name for easier identification (same as -J)\n".format(
                i))
        file.write(
            "#SBATCH --error=POET_rep_{}.err      # you can give your job a name for easier identification (same as -J)\n".format(
                i))

        # SBATCH --error=%j.err
        file.write("\n########## Command Lines to Run ##########\n\n")
        file.write("module pandas\n")
        # file.write("module load GCC/6.4.0-2.28 OpenMPI  ### load necessary modules, e.g\n")
        file.write("cd ~/POET1\n")
        file.write("srun -n 1 python poet.py -r {} -o {} -seed {}\n".format(runs, i, i + startingSeed))
        file.write("cd batch\n")
        file.write("scontrol show job $SLURM_JOB_ID     ### write job information to output file")
        file.close()
        call(["sbatch", "subs/{}.sb".format(i)])


main()

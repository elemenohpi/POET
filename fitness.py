# Authors: Iliya Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import random as R
import individual as I
import math
from threading import Thread
from scipy import stats
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

from tqdm import tqdm
import re
import statistics

import matplotlib.pyplot as plt



class Fitness:

    def __init__(self, config):
        self.config = config

        # Open training file
        try:
            self.learn_file = pd.read_csv(f"data/Training/{config['learn_data']}.csv", sep=';') # data/learnall.csv by default
            print(f'[INFO] Training data from {config["learn_data"]} are correctly loaded')
        except AttributeError:
            print("[ERROR] 'learn_data' in config.ini is not set")
            exit()

        self.mode = int(config["pattern_mode"])  # 0 is the summation mode and 1 is the multiplication mode
                                                 # of the motif weights
        self.k = 0  # for 10-fold cross validation implementation


    def model_vs_dataset(self, config, individuals):
        seq_fitness_tuples = []
        individuals_evaluations = []
        values = []

        for i, row in self.learn.iterrows():
            # row[0] is seq and row[1] is fitness
            seq_fitness_tuples.append((row[0], row[1]))
            values.append(row[1])

        for individual in individuals:
            predictions = []
            responses = []
            for seq_fit_tuple in seq_fitness_tuples:
                error, prediction = self.eval(seq_fit_tuple[0], seq_fit_tuple[1], individual, True)
                predictions.append(prediction)
            if config["fitness_alg"] == "correlation":
                align = np.polyfit(predictions, values, 1)
                for prediction in predictions:
                    response = prediction * align[0] + align[1]
                    responses.append(response)
                individuals_evaluations.append(responses)
            elif config["fitness_alg"] == "RMSE":
                individuals_evaluations.append(predictions)
        # print(individuals_evaluations)
        return seq_fitness_tuples, individuals_evaluations


    def create_dico_predictions(self, rule, sequence):
        try:
            regex_rule = re.compile(rule.pattern) # la compilation des regex est redondante
        except :
            print('[ERROR] Compilation error: fitness.py 93', rule.pattern, rule.tree_shape)
            exit()

        # Search match in 'sequence'. Non-overlapping matches of the expression
        object_find = re.finditer(regex_rule, sequence)

        for match in object_find: # search motif len >=2 et <=6
            motif = match.group()


            if len(motif) >= 2:
                if motif not in rule.db_motif.keys():
                    rule.db_motif[motif] = 1
                else:
                    rule.db_motif[motif] += 1


    def motif_searcher(self, sequence, individual, dbmotif, training_fitness):
        # On evalue chaque regles de l'individu
        for rule in individual.rules: # List of rules of the indiivdual
            self.create_dico_predictions(rule, sequence)









    def eval_RMSE(self, sequence, actualFitness, individual):
        measuredFitness = 0.0  # Fitness of individual i

        for rule in individual.rules:

            try:
                regex_rule = re.compile(rule.pattern) # la compilation des regex est redondante
            except :
                print('[ERROR] Compilation error: fitness.py 93', rule.pattern, rule.tree_shape)
                exit()

            # Search match in 'sequence'. Non-overlapping matches of the expression
            object_find = re.finditer(regex_rule, sequence)

            counter = 0

            for match in object_find: # search motif len >=2 et <=6
                motif = match.group()
                counter+=1 # il y a match

            if counter != 0: # comme il y a match on change le status
                if rule.status == 0:
                    rule.status = 1

                measuredFitness += rule.weight # on ajoute le poids a la valeur de fitness
        
        # Calculate the error after the estimation
        error = abs(measuredFitness - float(actualFitness))
        return error, measuredFitness


    def eval_motif_searcher(self, sequence, actualFitness, individual):
        '''
        function to evaluate each regex on the different sequence from the dataset

        sequence: a training peptide sequence from the training file
        actualFitness:
        individual: the target individual for the evaluation

        return: 
        '''

        measuredFitness = 0.0  # Fitness of individual i

        # On evalue chaque regles de l'individu
        for rule in individual.rules: # List of rules of the indiivdual
            
            try:
                regex_rule = re.compile(rule.pattern) # la compilation des regex est redondante
            except :
                print('[ERROR] Compilation error: fitness.py 93', rule.pattern, rule.tree_shape)
                exit()

            # Search match in 'sequence'. Non-overlapping matches of the expression
            object_find = re.finditer(regex_rule, sequence)

            full = 0        # check if complete sequence is found
            list_motif = [] # limit the redondance of motif
            nbr_match_motif = 0

            for match in object_find: # search motif len >=2 et <=6
                motif = match.group()

                # Find a good motif
                if len(motif) > 1 and len(motif) <= 11:
                    full += len(motif)

                    # remise a 0 si les jokers sont repetes .{X}
                    if '.{' in str(regex_rule) or '.+' in str(regex_rule) or '(..){' in str(regex_rule):
                        rule.score -= 10
                    else:
                        # motif not already found in the sequence
                        if motif not in list_motif:
                            nbr_match_motif+=1
                            list_motif.append(motif)
                            rule.score += (len(motif)*len(motif))
                            # print(sequence, match.span(),'---', motif, '+',len(motif)*len(motif), 'points' )

                            if rule.status == 0:
                                rule.status = 1
                                individual.usedRulesCount += 1

                # Size 1 - Singleton
                elif len(match.group()) == 1:
                    # print('taille 1',match,'---', match.span(),'---', match.group(), '-1 point')
                    rule.score -= 10

            # match with all the sequence
            if full == len(sequence):
                rule.score-=len(sequence)

            # Aucun match n'a ete trouvé (Penalty) - les individus avec + de rules seront penalise
            if nbr_match_motif == 0:
                rule.score -= 1 

            # print(f'SCORE de {rule.pattern} sur la sequence:',sequence, '=', rule.score)
            measuredFitness += rule.score

        return measuredFitness


    def eval(self, sequence, actualFitness, individual, returnPrediction=False):
        pos = 0
        measuredFitness = 0.0  # Fitness of individual i

        while pos < len(sequence):
            # Check every rule
            for rule in individual.rules:
                # Find the length of the rule
                try:
                    length_pattern = len(rule.pattern)
                    # print(rule.pattern)

                except:
                    print('Empty rule')
                    # happens in the case of empty rules.
                    continue

                # Continue if the rule length is more than the unchecked sequence length    
                if length_pattern > (len(sequence) - pos) or length_pattern == 0:
                    continue
                

                reverse_pattern = rule.pattern[::-1]
                # Normal or Reverse
                # print( rule.pattern, "vs", sequence[pos: (pos + length_pattern)])
                if ((rule.pattern == sequence[pos: (pos + length_pattern)]) or reverse_pattern == sequence[
                                                                                          pos: (pos + length_pattern)]):
                    # rule is found. Update its status and the usedRulesCount to avoid further computation

                    if rule.status == 0:
                        rule.status = 1
                        individual.usedRulesCount += 1

                    # Check the multiplication/summation mode. We always use the summation mode tho. 
                    if self.mode == 0:
                        measuredFitness += rule.weight
                    elif self.mode == 1:
                        measuredFitness *= rule.weight
                    else:
                        raise "wrong mode"

                    # If the rule is found, we don't wanna check any other "shorter" rules on the very same position. We wanna update the positioning and look for other rules. If no rule was found, the position updates anyway after the end of this loop so the difference is only the break command which happens only if the rule is found. This defenitely takes more processing power and considers overlapping rules but this is the way to go to consider all the possible cases. Also, note that this raises the issue of exponential slow-down when the number of rules in a table increase. Therefore, I assume we need some sort of bloat removal when the speed of the tool has decreased drastically.
                    # break
                    # print(rule.pattern, rule.weight)

            # Update the position 
            pos += 1
        
        # Calculate the error after the estimation
        error = abs(measuredFitness - actualFitness)
        # exit()
        if returnPrediction is True:
            return error, measuredFitness
        return error


    def eval_old(self, sequence, actualFitness, individual, returnPrediction=False):
        # This is the starting position of the sequence that we're looking at each time.
        pos = 0
        measuredFitness = 0.0  # Fitness of individual i

        # Iterate through the sequence
        while pos < len(sequence):
            # print(pos, len(sequence))

            # Check every rule
            for rule in individual.rules:
                # Find the length of the rule
                try:
                    length_pattern = len(rule.pattern)
                    # print(rule.pattern)

                except:
                    print('Empty rule')
                    # happens in the case of empty rules.
                    continue

                # Continue if the rule length is more than the unchecked sequence length	
                if length_pattern > (len(sequence) - pos) or length_pattern == 0:
                    continue
                

                reverse_pattern = rule.pattern[::-1]
                # Normal or Reverse
                # print( rule.pattern, "vs", sequence[pos: (pos + length_pattern)])
                if ((rule.pattern == sequence[pos: (pos + length_pattern)]) or reverse_pattern == sequence[
                                                                                          pos: (pos + length_pattern)]):
                    # rule is found. Update its status and the usedRulesCount to avoid further computation

                    if rule.status == 0:
                        rule.status = 1
                        individual.usedRulesCount += 1

                    # Check the multiplication/summation mode. We always use the summation mode tho. 
                    if self.mode == 0:
                        measuredFitness += rule.weight
                    elif self.mode == 1:
                        measuredFitness *= rule.weight
                    else:
                        raise "wrong mode"

                    # If the rule is found, we don't wanna check any other "shorter" rules on the very same position. We wanna update the positioning and look for other rules. If no rule was found, the position updates anyway after the end of this loop so the difference is only the break command which happens only if the rule is found. This defenitely takes more processing power and considers overlapping rules but this is the way to go to consider all the possible cases. Also, note that this raises the issue of exponential slow-down when the number of rules in a table increase. Therefore, I assume we need some sort of bloat removal when the speed of the tool has decreased drastically.
                    # break
                    # print(rule.pattern, rule.weight)

            # Update the position 
            pos += 1

        # Calculate the error after the estimation
        error = abs(measuredFitness - actualFitness)
        # exit()
        if returnPrediction is True:
            return error, measuredFitness
        return error

    def resetIndividual(self, individual):
        '''
        Reset values (usedRulesCount and rules.status) of an individual
        '''
        individual.usedRulesCount = 0
        individual.fitness = 0 # set the fitness

        # Make sure that every rule status is 0
        for rule in individual.rules:
            rule.status = 0
            rule.score = 0
            rule.db_motif = {}

    def measure_dataset(self, individual):
        train_error = 0.0
        for i, row in self.learn.iterrows():
            sequence = row[0]
            actual_fitness = row[1]
            error = self.eval(sequence, actual_fitness, individual)
            train_error += error ** 2

        trainSize = len(self.learn.index)
        MSE_train = train_error / trainSize
        RMSE_train = math.sqrt(MSE_train)

        return RMSE_train, 0

    def measureTotal(self, individual, dbmotif):
        """
        Mesure the fitness value for 1 individual
        """

        self.resetIndividual(individual) # Reset rules.status/rule.score/fitness of an individual
        
        # print('je vais evaluer:')
        # individual.print()
        # print(f'avec une fitness: {individual.fitness}')
        
        if self.config["fitness_alg"] == "RMSE":

            chunks = [0] * 10
            print('Chunks',chunks)
            for i, row in self.learn.iterrows():
                sequence = row[0]
                actual_fitness = row[1]
                j = int(i / (len(self.learn.index) / 10))
                if j == 10:
                    j = 9
                if j == 11:
                    print(i, int(len(self.learn.index)), int(i / (len(self.learn.index) / 10)))
                    exit()
                error = self.eval(sequence, actual_fitness, individual)
                chunks[j] += error ** 2

            test = 0
            train = 0
            for i in range(10):
                test += chunks[i]
                for j in range(10):
                    if i == j:
                        continue
                    train += chunks[j]

            testSize = int(len(self.learn.index) / 10)
            trainSize = len(self.learn.index) - testSize

            train = train / trainSize / 10
            test = test / testSize / 10

            RMSE_train = math.sqrt(train)
            RMSE_test = math.sqrt(test)
            return RMSE_train, RMSE_test
        elif self.config["fitness_alg"] == "correlation":
            predictions = []
            CEST_measurements = []
            for i, row in self.learn.iterrows():
                sequence = row[0]
                actual_fitness = row[1]
                CEST_measurements.append(actual_fitness)
                _, prediction = self.eval(sequence, actual_fitness, individual, True)
                predictions.append(prediction)
            pearsonr = stats.pearsonr(predictions, CEST_measurements)
            if math.isnan(pearsonr[0]):
                return 1, 0
            return 1 - pearsonr[0] ** 2, 0
        elif self.config["fitness_alg"] == "regex_motif":

            for i, row in self.learn_file.iterrows(): # training data
                # Training sequences and fitness values
                sequence, training_fitness = row[0], row[1] 
                
                # on evalue chaque rule de l'individu pour avoir 1 score par rule 
                # evaluate the individual
                score = self.eval_motif_searcher(sequence, training_fitness, individual)
            
            # change the fitness of the individual
            individual.fitness = score
        elif self.config["fitness_alg"] == "regex":
            predictions = []
            CEST_measurements = []

            for i, row in self.learn_file.iterrows():
                sequence, actual_fitness = row[0], row[1]
                print(sequence, actual_fitness)

                CEST_measurements.append(actual_fitness)

                _, prediction = self.eval(sequence, actual_fitness, individual, True)

                predictions.append(prediction)
            # pearsonr = stats.pearsonr(predictions, CEST_measurements)
            # if math.isnan(pearsonr[0]):
            #     return 1, 0
            # return 1 - pearsonr[0] ** 2, 0

        elif self.config["fitness_alg"] == "motif":
            for i, row in self.learn_file.iterrows(): # training data
                if ',' in row[0]:
                    row = row[0].strip().split(',')
                elif ';' in row:
                    row = row[0].strip().split(';')
                # Training sequences and fitness values
                sequence = row[0]
                for j, rule in enumerate(individual.rules): # List of rules of the indiivdual
                    # Cree la db_motif prediction pour chaque regle/regex sur chaque seq d'entrainement
                    self.create_dico_predictions(rule, sequence)
                



            for rule in individual.rules: # List of rules of the indiivdual
                # print(rule.db_motif)
                # print("REGEX:", rule.pattern)

                for motif, occ in rule.db_motif.items():


                        # !!! PRENDRE EN COMPTE LA LONGUEUR DU MOTIF,
                        # PLUS IL EST LONG MIEUX C EST


                    if motif in dbmotif.keys(): # la dbmotif du training set = tous les motifs du training set
                        if dbmotif[motif][2] == 1: # cest > 12.5
                            diff = abs(dbmotif[motif][0] - occ)
                            # print("diff", diff)
                            indice = round(dbmotif[motif][1],2)
                            # print("indice", indice)
                            x = diff / (indice/10)

                            # # Alternative
                            y = occ / dbmotif[motif][0]
                            y = 1 - abs(math.log10(y))
                            y = y * indice
                            
                            rule.score += y



                        else: # bad cest
                            diff = abs(dbmotif[motif][0] - occ)
                            # print("diff", diff)
                            indice = round(dbmotif[motif][1],2)
                            # print("indice", indice)
                            x = diff / (indice/10)

                            # Alternative
                            y = occ / dbmotif[motif][0]
                            y = 1 - abs(math.log10(y))
                            y = y * indice

                            # if x <= 0.0:
                                # rule.score -= indice
                            # else:
                                # rule.score -= (1 - math.log(x))
                            rule.score -= y*3

                        



                        # print("x", x)

                        # print(rule.pattern,motif,occ, dbmotif[motif][0], round(dbmotif[motif][1],2))
                        # print(x)


                            # print(diff, indice, x, (1-math.log(x)) )
                        # print(abs(occ-dbmotif[motif][0]))
                        # rule.score += abs(occ-dbmotif[motif][0]) / round(dbmotif[motif][1],2)
                    # else:
                    #     rule.score -= 0.5

                # La somme des score des regex = fitness total de l'individu
                individual.fitness += rule.score
        elif self.config["fitness_alg"] == "RMSE2":
            predictions = []
            CEST_measurements = []

            for i, row in self.learn_file.iterrows():
                if ',' in row[0]:
                    row = row[0].strip().split(',')
                elif ';' in row:
                    row = row[0].strip().split(';')

                sequence = row[0]
                actual_fitness = float(row[1])

                # Add all CEST value
                CEST_measurements.append(actual_fitness)

                # evaluate individual on the sequence
                error, measuredFitness = self.eval_RMSE(sequence, actual_fitness, individual)

                predictions.append(measuredFitness)

            pearsonr = stats.pearsonr(predictions, CEST_measurements)
            if math.isnan(pearsonr[0]):
                individual.fitness = 1.0

            # return 1 - pearsonr[0] ** 2, 0

            individual.fitness = round(1 - pearsonr[0] ** 2,3)





    def predict(self, sequence, individual):
        raise "Make use of eval() instead of predict"
        fitness = 0.0
        pos = 0

        # Iterate through the sequence
        while pos < len(sequence):
            # Check every rule
            found = False
            for i, rule in enumerate(individual.rules):
                if rule.status == 0:
                    continue

                try:
                    length = len(rule.pattern)
                except:
                    continue
                # Continue if the rule length is more than the unchecked sequence length
                if length > (len(sequence) - pos) or length == 0:
                    continue
                reverse_pattern = rule.pattern[::-1]
                if (pos + length < len(sequence) and (
                        rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                             pos: (pos + length)]):
                    # print("pos is at " + str(pos) + " sequence is: " + sequence + " found: " + sequence[pos : (pos + length)] + " weight: " + str(rule.weight))
                    if self.mode == 0:
                        fitness += rule.weight
                    elif self.mode == 1:
                        fitness *= rule.weight
                    break
            pos += 1
        return fitness


def plot_fitness():
    CEST_value = [50,40,30,20,12.5,10,5,1,0] #9
    difference = [10,9,8,7,6,5,4,3,2,1] #11 
   
    for diff in difference:
        print('--- DIFF = ',diff, '---')
        l = []
        for cest in CEST_value:

            x = 1-math.log10(diff)
            print(diff, x, math.log10(diff)) # log neperien
            x = x * cest
            print(cest, x)
            l.append(x)

        plt.plot(CEST_value, l, label=str(diff))
    plt.legend()
    plt.xlabel('CEST value')
    plt.ylabel('Score')
    plt.grid()
    plt.show()

# plot_fitness()



                            
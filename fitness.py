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
import math

import threading
from multiprocessing import Lock, Process

import os


import motif_bank




class Fitness:

    def __init__(self, config):
        self.config = config

        # Open training file
        try:
            self.learn_file = pd.read_csv(f"data/Training/{config['learn_data']}.csv", sep=';')
            print(f'[INFO] Training data from {config["learn_data"]} are correctly loaded')
        except AttributeError:
            print("[ERROR] 'learn_data' in config.ini is not set")
            exit()

        self.k = 0  # for 10-fold cross validation implementation
        self.min_size_motif = int(config['min_motif_size_db'])
        self.max_size_motif = int(config['max_motif_size_db'])
        self.penalty = 5
        self.max_rules = int(config['maximum_rule_count'])




        self.training_db_motif = {} # used only with the fitness function 'motif'
        self.db_strat = config["db_motif"]
        self.alphabet = config["learn_data"]
        self.MAX_motif_size = int(config["max_motif_size_db"])
        self.MIN_motif_size = int(config["min_motif_size_db"])

        self.build_db_motif()





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


    def create_dico_predictions(self, list_rules, sequence):
        '''
        This function add to the rule dictionnary, all motifs identified
        by the regex in the sequence

        rule: a regex
        sequence: the sequence to analyze with the regex 
        '''

        for rule in list_rules:
            try:
                regex_rule = re.compile(rule.pattern)
            except :
                print('[ERROR] Compilation error: fitness.py 93', rule.pattern, rule.tree_shape)
                exit()
    
            # Search match in 'sequence'. Non-overlapping matches of the expression    
            for match in re.finditer(regex_rule, sequence): # search motif len >=2 et <=6
                motif = match.group()
    
                if self.min_size_motif <= len(motif) <= self.max_size_motif:
                    rule.nb_match+=1
                    rule.db_motif[motif] = rule.db_motif.get(motif, 0) + 1




        


    # def motif_searcher(self, sequence, individual, dbmotif, training_fitness):
    #     # On evalue chaque regles de l'individu
    #     for rule in individual.rules: # List of rules of the indiivdual
    #         self.create_dico_predictions(rule, sequence)




    # Build the training database with all motifs

    def build_db_motif(self):
        for i in range (self.MIN_motif_size, self.MAX_motif_size+1):
            # Strategies
            if self.db_strat == 'std':  # Best CEST value
                self.training_db_motif = motif_bank.standard_motif_bank(i, self.alphabet)
            if self.db_strat == 'avg':  # Average CEST value 
                self.training_db_motif = motif_bank.average_motif_bank(i, self.alphabet)
            if self.db_strat == 'occ':  # Number of motif in each class define the classe
                self.training_db_motif = motif_bank.occ_motif_bank(i, self.alphabet)



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
        actualFitness: the fitness/value of the sequence
        individual: the target individual for the evaluation

        return: measuredFitness
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
                            individual.nb_match +=1
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
            rule.nb_match = 0

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



    def logistic_function(self, x, k, x0, L):
        return 1/(L / (np.exp(-k*(x-x0))))


    # Anciennement measureTotal()
    def evaluate_individual(self, individual):
        """
        Mesure the fitness value for 1 individual
        """
        dbmotif = self.training_db_motif





        
        # Reset rules.status/rule.score/fitness of an individual
        self.resetIndividual(individual) 
        
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
        # elif self.config["fitness_alg"] == "correlation":
                # predictions = []
                # CEST_measurements = []
                # for i, row in self.learn.iterrows():
                #     sequence = row[0]
                #     actual_fitness = row[1]
                #     CEST_measurements.append(actual_fitness)
                #     _, prediction = self.eval(sequence, actual_fitness, individual, True)
                #     predictions.append(prediction)
                # pearsonr = stats.pearsonr(predictions, CEST_measurements)
                # if math.isnan(pearsonr[0]):
                #     return 1, 0
                # return 1 - pearsonr[0] ** 2, 0
        # elif self.config["fitness_alg"] == "regex_motif":
    
                # for i, row in self.learn_file.iterrows(): # training data
                #     # Training sequences and fitness values
                #     sequence, training_fitness = row[0], row[1] 
                    
                #     # on evalue chaque rule de l'individu pour avoir 1 score par rule 
                #     # evaluate the individual
                #     score = self.eval_motif_searcher(sequence, training_fitness, individual)
                
                # # change the fitness of the individual
                # individual.fitness = score
        # elif self.config["fitness_alg"] == "regex":
                # predictions = []
                # CEST_measurements = []
    
                # for i, row in self.learn_file.iterrows():
                #     sequence, actual_fitness = row[0], row[1]
                #     print(sequence, actual_fitness)
    
                #     CEST_measurements.append(actual_fitness)
    
                #     _, prediction = self.eval(sequence, actual_fitness, individual, True)
    
                #     predictions.append(prediction)
                # # pearsonr = stats.pearsonr(predictions, CEST_measurements)
                # # if math.isnan(pearsonr[0]):
                # #     return 1, 0
                # # return 1 - pearsonr[0] ** 2, 0
        # elif self.config["fitness_alg"] == "RMSE2":
            # predictions = []
            # CEST_measurements = []

            # for i, row in self.learn_file.iterrows():
            #     if ',' in row[0]:
            #         row = row[0].strip().split(',')
            #     elif ';' in row:
            #         row = row[0].strip().split(';')

            #     sequence = row[0]
            #     actual_fitness = float(row[1])

            #     # Add all CEST value
            #     CEST_measurements.append(actual_fitness)

            #     # evaluate individual on the sequence
            #     error, measuredFitness = self.eval_RMSE(sequence, actual_fitness, individual)

            #     predictions.append(measuredFitness)

            # pearsonr = stats.pearsonr(predictions, CEST_measurements)
            # if math.isnan(pearsonr[0]):
            #     individual.fitness = 1.0

            # # return 1 - pearsonr[0] ** 2, 0

            # individual.fitness = round(1 - pearsonr[0] ** 2,3)
        elif self.config["fitness_alg"] == "motif":
            fitness = 0
            # Main fitness function based on motifs

            # Create the db motif for each rule/regex of each individual on each training sequence
            for i, row in self.learn_file.iterrows(): # training data
                delimiter = ',' if ',' in row[0] else ';'
                row = row[0].strip().split(delimiter)
                sequence = row[0] # Training sequences

                # Create the motif db for each rule on each sequence
                self.create_dico_predictions(individual.rules, sequence)

            nbr_regex = len(individual.rules)
            # print('Nombre de regle', x)
            for rule in individual.rules: 

                rule.score = 0

                for motif, occ in rule.db_motif.items():
                    if motif in dbmotif.keys(): # la dbmotif du training set = tous les motifs du training set
    
                        motif_occ = dbmotif[motif][0]
                        weight = round(dbmotif[motif][1],2)
                        motif_class = dbmotif[motif][2]
                        len_motif = len(motif)

                        y = self.fitness_function(occ, motif_occ, weight, len_motif)

                        # Class of the motif (1=good, 0=bad)
                        if motif_class == 1: # cest >= 12.5                            
                            rule.score += y /10
                        else: # # cest > 12.5
                            rule.score -= (y * self.penalty) /10


                # std_dev = 5 # bien definir si on change le nombre max de regle
                # mean = int(self.max_rules/2)

                # print('Score de la regle', rule.score)

                # normal distribution
                # pds = (1/(std_dev * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((nbr_regex - mean) / std_dev) ** 2)

                # Logistic regression
                # print('NBR regex:', len(individual.rules))
                x = len(individual.rules) # nbr_regex_max
                k = 0.05

                pds = self.logistic_function(x, k, 0, 1)
                # print('PDS',pds)

                

                # pds = pds*10
                # print('Poids:', pds)

                # print('Rule',len(individual.rules), rule.score, pds, (rule.score*pds))
                # print('Score/pds = ', rule.score/pds, rule.score)

                # utile ???


                if rule.score == 0:
                    rule.score = -10 * pds
                    # score = -1000
                else:
                    rule.score = rule.score * pds

                # individual.fitness += rule.score
                fitness += rule.score

            individual.fitness = fitness
            # print('Fitness:',fitness)
            # individual.print()

        return individual





        


    def fitness_function(self, occ_motif, occ_db, weight, len_motif):
        y = occ_motif / occ_db
        y = 1 - abs(math.log10(y))
        y = (y * weight) * len_motif # *len_motif permet de prendre en compte la taille des motifs

        return round(y,2)


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



                            
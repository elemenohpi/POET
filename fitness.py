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



class Fitness:

    def __init__(self, config):
        self.config = config

        # Open training file
        try:
            self.learn_file = pd.read_csv(config["learn_data"]) # data/learnall.csv by default
            print(f'[INFO] Training data from {config["learn_data"]} are correctly loaded')
        except AttributeError:
            print("[ERROR] settings.learn_df is not set")
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




    def regex_eval(self, sequence, actualFitness, individual):
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
            
            nbr_match_motif = 0 # plus le nombre de match est grand, moins c'est bien
            nbr_singleton = 0 # Match d'une seule lettre

            singleton = False

            try:
                regex_rule = re.compile(rule.pattern) # la compilation des regex est redondante
            except :
                print('[ERROR] fitness.py 93', rule.pattern, rule.tree_shape)
                exit()

            # print('Evaluation de la regex',i, '-', rule.pattern)

            # cherche les matches dans la sequence cible - correspondances non chevauchantes de l'expression
            object_find = re.finditer(regex_rule, sequence)
            full = 0

            for match in object_find: # search motif len >=2 et <=6

                motif = match.group()

                if len(motif) > 1 and len(motif) <= 6:
                    # print(match.span(),'---', motif, '+',len(motif), 'points' )
                    nbr_match_motif+=1
                    full += len(match.group())
                    rule.score += len(match.group())

                    if rule.status == 0:
                        rule.status = 1
                        individual.usedRulesCount += 1

                # Size 1
                elif len(match.group()) == 1:
                    nbr_singleton+=1
                    # print('taille 1',match,'---', match.span(),'---', match.group(), '-1 point')
                    singleton = True

            if full == len(sequence):
                rule.score-=len(sequence)

            # Aucun match n'a ete trouvé (Penalty) - les individus avec + de rules seront penalise
            # if nbr_match_motif == 0 and nbr_singleton==0:
            #     rule.score -= 1 

            if singleton:
                # print('Singleton -1')
                rule.score -= 1

            # print(f'SCORE de {rule.pattern} sur la sequence:',sequence, '=', rule.score)
            measuredFitness += rule.score

        return measuredFitness


    def eval(self, sequence, actualFitness, individual, returnPrediction=False):
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

    def measureTotal(self, individual):
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

        elif self.config["fitness_alg"] == "regex":
            CEST_measurements = []

            for i, row in self.learn_file.iterrows(): # training data
                # Training sequences and fitness values
                sequence, training_fitness = row[0], row[1] 

                # Add training (True) fitness values (repetitive !)
                CEST_measurements.append(training_fitness)  
                
                # print('sequence', i, '=', sequence, 'fitness=', training_fitness)
                # on evalue chaque rule de l'individu pour avoir 1 score par rule 
                # et donc un score totale ou moyen final
                score = self.regex_eval(sequence, training_fitness, individual)
            
            individual.fitness = score

        



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

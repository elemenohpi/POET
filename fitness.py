# Authors: Iliya Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import random as R
import individual as I
import math
from threading import Thread
from scipy import stats
import numpy as np

import warnings
warnings.filterwarnings('ignore')

class Fitness:
    def __init__(self):
        try:
            self.learn = settings.learn_df
        except AttributeError:
            print("WARNING: settings.learn_df is not set")

        self.mode = settings.pattern_mode  # 0 is the summation mode and 1 is the multiplication mode
        self.k = 0  # for 10-fold cross validation implementation
        pass

    def model_vs_dataset(self, individuals):
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
            if settings.fitness_alg == "correlation":
                align = np.polyfit(predictions, values, 1)
                for prediction in predictions:
                    response = prediction * align[0] + align[1]
                    responses.append(response)
                individuals_evaluations.append(responses)
            elif settings.fitness_alg == "RMSE":
                individuals_evaluations.append(predictions)
        # print(individuals_evaluations)
        return seq_fitness_tuples, individuals_evaluations

    def eval(self, sequence, actualFitness, individual, returnPrediction=False):
        # This is the starting position of the sequence that we're looking at each time. 
        pos = 0
        measuredFitness = 0.0  # Fitness of individual i

        # Iterate through the sequence
        while pos < len(sequence):
            # Check every rule
            for rule in individual.rules:
                # Find the length of the rule
                try:
                    length = len(rule.pattern)
                except:
                    # happens in the case of empty rules.
                    continue

                # Continue if the rule length is more than the unchecked sequence length	
                if length > (len(sequence) - pos) or length == 0:
                    continue

                reverse_pattern = rule.pattern[::-1]
                # Normal or Reverse
                if ((rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                          pos: (pos + length)]):
                    # rule is found. Update its status and the usedRulesCount to avoid further computation
                    if rule.status == 0:
                        rule.status = 1
                        individual.usedRulesCount += 1

                    # Check the multiplication/summation mode. We always use the summation mode tho. 
                    if self.mode == 0:
                        measuredFitness += rule.weight
                    elif self.mode == 1:
                        measuredFitness *= rule.weight

                    # If the rule is found, we don't wanna check any other "shorter" rules on the very same position. We wanna update the positioning and look for other rules. If no rule was found, the position updates anyway after the end of this loop so the difference is only the break command which happens only if the rule is found. This defenitely takes more processing power and considers overlapping rules but this is the way to go to consider all the possible cases. Also, note that this raises the issue of exponential slow-down when the number of rules in a table increase. Therefore, I assume we need some sort of bloat removal when the speed of the tool has decreased drastically.
                    # break
                    # print(rule.pattern, rule.weight)

            # Update the position 
            pos += 1

        # Calculate the error after the estimation
        error = abs(measuredFitness - actualFitness)
        if returnPrediction is True:
            return error, measuredFitness
        return error

    def resetIndividual(self, individual):
        individual.usedRulesCount = 0

        # Make sure that every rule status is 0
        for rule in individual.rules:
            rule.status = 0

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
        self.resetIndividual(individual)
        if settings.fitness_alg == "RMSE":
            # ToDo:: this is not using threading. Change the setting name
            if settings.enable_threading == "yes":
                chunks = [0] * 10
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
            else:
                raise "uncharted territory"
        elif settings.fitness_alg == "correlation":
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
        else:
            print(settings.fitness_alg)
            raise "Uncharted territory"





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

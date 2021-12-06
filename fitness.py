# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import random as R
import settings
import individual as I
import math
from threading import Thread

class Fitness:
    def __init__(self):
        self.learn = settings.learn_df
        self.mode = settings.pattern_mode  # 0 is the summation mode and 1 is the multiplication mode
        self.k = 0  # for 10-fold cross validation implementation
        pass

    def eval(self, sequence, actualFitness, individual):
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
                if length < (len(sequence) - pos) or length == 0:
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
                    break

            # Update the position 
            pos += 1

        # Calculate the error after the estimation
        error = (measuredFitness - actualFitness) ** 2
        return error

    def resetIndividual(self, individual):
        individual.usedRulesCount = 0

        # Make sure that every rule status is 0
        for rule in individual.rules:
            rule.status = 0

    def thread_measure(self, individual, j, results):
        test_error = 0.0
        train_error = 0.0
        for i, row in self.learn.iterrows():
            sequence = row[0]
            actual_fitness = row[1]
            if int(self.learn.size / 10) * j <= i < int(self.learn.size / 10) * (j + 1):
                test_error += self.eval(sequence, actual_fitness, individual)
                continue
            train_error += self.eval(sequence, actual_fitness, individual)
        results[j] = (train_error, test_error)

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


        # if settings.enable_threading == "yes":
        #     threads = [None] * 10
        #     results = [None] * 10
        #
        #     for i in range(10):
        #         # threads[i] = Thread(target=self.thread_measure, args=(individual))
        #         threads[i] = Thread(target=self.thread_measure, args=(individual, i, results))
        #         threads[i].start()
        #         pass
        #
        #     for i in range(10):
        #         threads[i].join()
        #
        #     MSE_train = 0
        #     MSE_test = 0
        #
        #     for i in range(10):
        #         MSE_train += results[i][0]
        #         MSE_test += results[i][1]
        #
        #     testSize = int(int(self.learn.size) / 10)
        #     trainSize = self.learn.size - int(int(self.learn.size) / 10)
        #     MSE_train = MSE_train / trainSize / 10
        #     MSE_test = MSE_test / testSize / 10
        #     RMSE_train = math.sqrt(MSE_train)
        #     RMSE_test = math.sqrt(MSE_test)
        #
        #     return RMSE_train, RMSE_test
        # elif settings.enable_threading != "no":
        #     raise "incorrect value for enable_threading. correct options: no / yes / maybe"
        #
        # test_error = 0.0
        # train_error = 0.0
        # for j in range(10):
        #     # Iterating through the training set
        #     for i, row in self.learn.iterrows():
        #         sequence = row[0]
        #         actual_fitness = row[1]
        #         if int(self.learn.size / 10) * j <= i < int(self.learn.size / 10) * (j + 1):
        #             test_error += self.eval(sequence, actual_fitness, individual)
        #             continue
        #
        #         train_error += self.eval(sequence, actual_fitness, individual)
        #
        # testSize = int(int(self.learn.size) / 10)
        # trainSize = self.learn.size - int(int(self.learn.size) / 10)
        # MSE_train = train_error / trainSize / 10
        # MSE_test = test_error / testSize / 10
        # RMSE_train = math.sqrt(MSE_train)
        # RMSE_test = math.sqrt(MSE_test)
        #
        # return RMSE_train, RMSE_test

    def measure(self, individual):
        raise "update this function according to measureTotal"
        if self.k == 10:
            self.k = 0
        fitness = 0.0
        sum = 0.0
        for i, row in self.learn.iterrows():
            if i >= int(self.learn.size / 10) * self.k and i < int(self.learn.size / 10) * (self.k + 1):
                continue
            sequence = row[0]
            actualFitness = row[1]
            pos = 0
            measuredFitness = 0.0  # Fitness of individual i
            if self.mode == 1:
                measuredFitness = 1.0

            # Iterate through the sequence
            while pos < len(sequence):
                # Check every rule
                found = False
                for rule in individual.rules:
                    length = len(rule.pattern)
                    reverse_pattern = rule.pattern[::-1]
                    # Normal or Reverse
                    if (pos + length < len(sequence) and (
                            rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                                 pos: (pos + length)]):
                        # rule is found
                        individual.usedRules[rule.pattern] = True
                        # print("pos is at " + str(pos) + " sequence is: " + sequence + " found: " + sequence[pos : (pos + length)] + " weight: " + str(rule.weight))
                        if self.mode == 0:
                            measuredFitness += rule.weight
                        elif self.mode == 1:
                            measuredFitness *= rule.weight
                        pos = pos + length
                        found = True
                        break
                if not found:
                    pos += 1
                if found and length == 0:
                    pos += 1
            # print("pos: "+str(pos)+" len: "+str(len(sequence)))
            error = (measuredFitness - actualFitness) ** 2
            # ToDo:: Checkout limiting measuredFitness to possitive values

            fitness += error
        # print ("actualFitness: " + str(actualFitness) +" measuredFitness: "+str(measuredFitness)+" error: "+str(error)+" accuracy: "+ str(accuracy))
        # print("Training:: {}- Actual: {}, Estimate: {} Squared Error: {}".format(i, actualFitness, measuredFitness, error))
        MSE = fitness / int(self.learn.size - int(self.learn.size / 10))
        RMSE = math.sqrt(MSE)
        # print("ENDDDD:: MSE: {}, RMSE: {}".format(MSE, RMSE))
        return RMSE

    def validate(self, individual):
        raise ("update this function according to measureTotal")
        if self.k == 10:
            self.k = 0
        fitness = 0.0
        sum = 0.0
        for i, row in self.learn.iterrows():
            if not (i >= int(self.learn.size / 10) * self.k and i < int(self.learn.size / 10) * (self.k + 1)):
                continue
            sequence = row[0]
            actualFitness = row[1]
            pos = 0
            measuredFitness = 0.0  # Fitness of individual i
            if self.mode == 1:
                measuredFitness = 1.0

            # Iterate through the sequence
            while pos < len(sequence):
                # Check every rule
                found = False
                for rule in individual.rules:
                    length = len(rule.pattern)
                    reverse_pattern = rule.pattern[::-1]
                    # Normal or Reverse
                    if (pos + length < len(sequence) and (
                            rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                                 pos: (pos + length)]):
                        # rule is found
                        individual.usedRules[rule.pattern] = True
                        # print("pos is at " + str(pos) + " sequence is: " + sequence + " found: " + sequence[pos : (pos + length)] + " weight: " + str(rule.weight))
                        if self.mode == 0:
                            measuredFitness += rule.weight
                        elif self.mode == 1:
                            measuredFitness *= rule.weight
                        pos = pos + length
                        found = True
                        break
                if not found:
                    pos += 1
                if found and length == 0:
                    pos += 1
            # print("pos: "+str(pos)+" len: "+str(len(sequence)))
            error = (measuredFitness - actualFitness) ** 2
            fitness += error
        # print ("actualFitness: " + str(actualFitness) +" measuredFitness: "+str(measuredFitness)+" error: "+str(error)+" accuracy: "+ str(accuracy))
        # print("Validation:: {}- Actual: {}, Estimate: {} Squared Error: {}".format(i, actualFitness, measuredFitness, error))
        MSE = fitness / int(self.learn.size / 10)
        # print(self.learn.size / 10)
        RMSE = math.sqrt(MSE)
        # print("ENDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD:: MSE: {}, RMSE: {}".format(MSE, RMSE))
        return RMSE

    def measure_unseen(self, individual):
        raise ("update this function according to measureTotal")
        if self.unseen.size <= 0:
            return 0
        fitness = 0.0
        sum = 0.0
        for i, row in self.unseen.iterrows():
            sequence = row[0]
            actualFitness = row[1]
            pos = 0
            measuredFitness = 0.0  # Fitness of individual i
            if self.mode == 1:
                measuredFitness = 1.0

            # Iterate through the sequence
            while pos < len(sequence):
                # Check every rule
                found = False
                for rule in individual.rules:
                    length = len(rule.pattern)
                    # Not tested #ToDo::
                    reverse_pattern = rule.pattern[::-1]
                    if (pos + length < len(sequence) and (
                            rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                                 pos: (pos + length)]):
                        # print("pos is at " + str(pos) + " sequence is: " + sequence + " found: " + sequence[pos : (pos + length)] + " weight: " + str(rule.weight))
                        if self.mode == 0:
                            measuredFitness += rule.weight
                        elif self.mode == 1:
                            measuredFitness *= rule.weight
                        pos = pos + length
                        found = True
                        break
                if not found:
                    pos += 1
                if found and length == 0:
                    pos += 1
            error = (measuredFitness - actualFitness) ** 2
            fitness += error
        # print ("actualFitness: " + str(actualFitness) +" measuredFitness: "+str(measuredFitness)+" error: "+str(error)+" accuracy: "+ str(accuracy))
        MSE = fitness / self.learn.size
        RMSE = math.sqrt(MSE)
        return RMSE

    def predict(self, sequence, individual):
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
                reverse_pattern = rule.pattern[::-1]
                if (pos + length < len(sequence) and (
                        rule.pattern == sequence[pos: (pos + length)]) or reverse_pattern == sequence[
                                                                                             pos: (pos + length)]):
                    # print("pos is at " + str(pos) + " sequence is: " + sequence + " found: " + sequence[pos : (pos + length)] + " weight: " + str(rule.weight))
                    if self.mode == 0:
                        fitness += rule.weight
                    elif self.mode == 1:
                        fitness *= rule.weight
                    pos = pos + length
                    found = True
                    break
            if not found:
                pos += 1
            if found and length == 0:
                pos += 1
        return fitness

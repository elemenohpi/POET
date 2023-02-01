import os
import random

import fitness as F
import random as R
import individual as I
import copy
import pandas as pd

class Sequence:
    def __init__(self, pattern, fitness):
        self.pattern = pattern
        self.fitness = fitness
        pass


# all I know is that iter is a lie... #ToDo::
class Predictor:

    def __init__(self, count, size, iterations, config, model_path):
        self.config = config
        self.model_path = model_path
        self.count = count
        self.size = size
        self.iterations = iterations
        self.pop = []
        self.output = []
        self.popSize = 10000
        codes = pd.read_csv("data/translation/amino_to_amino.csv")
        codes = codes["code"].tolist()
        self.codes = codes
        self.hydroDic = {
            "I": -0.31,
            "L": -0.56,
            "F": -1.13,
            "V": 0.07,
            "M": -0.23,
            "P": 0.45,
            "W": -1.85,
            "J": 0.17,
            "T": 0.14,
            "E": 2.02,
            "Q": 0.58,
            "C": -0.24,
            "Y": -0.94,
            "A": 0.17,
            "S": 0.13,
            "N": 0.42,
            "D": 1.23,
            "R": 0.81,
            "G": 0.01,
            "H": 0.96,
            "K": 0.99
        }
        pass

    def predict(self):
        random.seed(int(self.config["seed"]))
        objF = F.Fitness(self.config)
        files = [f for f in os.listdir(self.model_path) if os.path.isfile(os.path.join(self.model_path, f))]
        ensemble = []
        for model in files:
            if model.split(".")[-1] != "csv":
                continue
            individual = I.Individual(self.config)
            individual.makeFromFile(os.path.join(self.model_path, model))
            # individual.remove_unexpressed()
            ensemble.append(individual)

        self.populate(ensemble)

        # self.sort()

        learn_data = pd.read_csv(self.config["learn_data"])
        data_set_rules = learn_data["sequence"].tolist()

        # for i in range(self.count):
        #     print(repr(i + 1) + ":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness, 2))

        for iter_index in range(1, self.iterations+1):
            if iter_index % 1 == 0:
                print("Iteration", iter_index, "complete")
            # self.sort()
            # for i in range(self.count):
            # 	print(repr(i+1)+":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness,2))

            improvement = False

            # mutate
            for myi, seq in enumerate(self.pop):
                randomsite = R.randint(0, self.size - 1)
                # print("randomsite", randomsite)
                rand = R.randint(0, len(self.codes) - 1)
                newseq = copy.deepcopy(seq)
                # print("seq", seq.pattern)
                # print(newseq.pattern)
                # exit()
                newamino = self.codes[rand]
                # print("newamino", newamino)
                # try:
                # newseq.pattern[randomsite] = newamino
                chars = list(newseq.pattern)
                chars[randomsite] = newamino
                newseq.pattern = "".join(chars)
                # except Exception(e):
                #     print(e)
                #     print("newseq pattern", newseq.pattern)
                #     print("newamino", newamino)
                #     print("randomsite", randomsite)
                #     exit()
                # print(seq.pattern)
                # print(newseq.pattern)
                if newseq.pattern in data_set_rules:
                    fitness = 0
                    print("\n already in the dataset \n")
                elif self.hydrophobic(newseq.pattern):
                    fitness = 0
                else:
                    fitness = 0
                    for model in ensemble:
                        _, tempF = objF.eval(newseq.pattern, 0, model, True)
                        fitness += tempF
                    fitness /= len(ensemble)

                # print("tempF", tempF)
                if fitness > seq.fitness:
                    # print("old seq", seq.pattern)
                    # print("new seq", newseq.pattern)
                    # print("fitness", fitness)
                    # print("old fitness", seq.fitness)
                    newseq.fitness = fitness
                    # print(seq.pattern + " -> " + newseq.pattern + " | fitness: " + repr(seq.fitness) + " -> " + repr(newseq.fitness) )
                    self.pop[myi] = newseq
                    # print(self.pop[myi].fitness)
                    # exit()
                    # improvement = True

            # if improvement:
            #     print(seq.fitness)
            #     seq = newseq
            #     print(seq.fitness)
            #     print(self.pop[myi].fitness, "should be updated")
            #     exit()
                # for i in range(self.count):
                #     print(repr(i + 1) + ":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness, 2))
        self.sort()
        for i in range(10):
            print(self.pop[i].pattern, self.pop[i].fitness)

    def sort_rules(self, rules):
        sortedR = []
        for rule in rules:
            placed = False
            for index, r in enumerate(sortedR):
                if rule.weight < r.weight:
                    sortedR.insert(index, rule)
                    placed = True
                    break
            if not placed:
                sortedR.append(rule)
        sortedR = reversed(sortedR)
        return sortedR

    # self.size
    # self.count

    def populate(self, ensemble):
        print("Populating the sequence pool...\n")
        # Read the best model and turn it to an individual

        objF = F.Fitness(self.config)
        for i in range(self.popSize):
            pattern = ""
            if self.size <= 0:
                self.size = R.randint(8, 12)
            for j in range(self.size):
                rand = R.randint(0, len(self.codes) - 1)
                pattern += self.codes[rand]
            if self.hydrophobic(pattern):
                fitness = 0
            else:
                fitness = 0
                for model in ensemble:
                    _, tempF = objF.eval(pattern, 0, model, True)
                    fitness += tempF
                fitness /= len(ensemble)

            sequence = Sequence(pattern, fitness)
            self.pop.append(sequence)
        # for seq in self.pop:
        #     if len(seq.pattern) != self.size:
        #         print(len(seq.pattern), seq.pattern)
        # exit()

    def hydrophobic(self, pattern):
        temp = 0
        for char in pattern:
            temp += self.hydroDic[char]
        if temp >= 0:
            return False
        else:
            return True

    def sort(self):
        print("Sorting the sequences...\n")
        n = len(self.pop)
        # Traverse through all array elements
        for i in range(n):
            # Last i elements are already in place
            for j in range(0, n - i - 1):
                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if self.pop[j].fitness < self.pop[j + 1].fitness:
                    self.pop[j], self.pop[j + 1] = self.pop[j + 1], self.pop[j]

# Restrict the amino acid search
# impirical research

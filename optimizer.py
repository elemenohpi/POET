import fitness as F
import pandas as pd
import archivist as Archivist
import copy
import random as R
import individual as I
import rule as Rule


class Optimizer:
    def __init__(self, config, population):
        self.config = config
        self.P = population
        self.runs = int(config["runs"])
        self.tournamentSize = int(config["tournament_size"])
        self.logInterval = int(config["pop_log_interval"])
        self.crossRate = float(config["crossover_unused_selection_chance"])
        self.ruleSize = int(config["maximum_rule_size"])
        self.ruleCount = int(config["maximum_rule_count"])
        self.minWeight = float(config["rule_weight_min"])
        self.maxWeight = float(config["rule_weight_max"])
        self.output_evo = config["output_evo"]
        self.output_model = config["output_model"]
        self.mAR = float(config["mut_add_rule"])
        self.mRR = float(config["mut_remove_rule"])
        self.mCW = float(config["mut_change_weight"])
        self.mATP = float(config["mut_add_to_pattern"])
        self.mRFP = float(config["mut_remove_from_pattern"])
        self.mCWmin = 0
        self.mCWmax = 1
        codes = pd.read_csv("data/translation/amino_to_amino.csv")
        self.codes = codes["code"].tolist()

    def optimize(self):
        fitness = F.Fitness(self.config)
        arch = Archivist.Archivist(self.config)

        num_workers = int(self.config.get("workers", 0))
        fitness.start_workers(num_workers)

        log_string = (
            "best fitness, best test, best rule count, best unused rule count, average fitness, "
            "average test, average rule count, average unused rule count "
        )
        arch.saveEvo(log_string)

        try:
            self._run_generations(fitness, arch)
        finally:
            fitness.stop_workers()

    def _run_generations(self, fitness, arch):
        for i in range(self.runs):
            # Evaluate all individuals (parallel or sequential)
            fitness.measureBatch(self.P.pop)

            avgFitness = 0.0
            avgTest = 0.0
            avgRuleCount = 0.0
            avgUsedRulesCount = 0.0
            bestIndividual = self.P.pop[0]

            for j in self.P.pop:
                avgFitness += j.fitness
                avgTest += j.test
                avgRuleCount += len(j.rules)
                avgUsedRulesCount += j.usedRulesCount
                if j.fitness < bestIndividual.fitness:
                    bestIndividual = j

            avgFitness = round(avgFitness / float(len(self.P.pop)), 3)
            avgTest = round(avgTest / float(len(self.P.pop)), 3)
            avgRuleCount = round(avgRuleCount / float(len(self.P.pop)), 0)
            avgUsedRulesCount = round(avgUsedRulesCount / float(len(self.P.pop)), 0)
            bestFitness = round(bestIndividual.fitness, 3)
            bestTest = round(bestIndividual.test, 3)
            bestRuleCount = len(bestIndividual.rules)
            bestUsedRulesCount = bestIndividual.usedRulesCount

            # Log the outcome before doing the changes to the population / generating a new population
            print_string = "{}: -b {} -bt {} -rc {} -urc {} ||| -a {} -at {} -arc {} -aurc {}".format(
                i,
                bestFitness,
                bestTest,
                bestRuleCount,
                bestUsedRulesCount,
                avgFitness,
                avgTest,
                avgRuleCount,
                avgUsedRulesCount,
            )

            log_string = "{},{},{},{},{},{},{},{},{}".format(
                i,
                bestFitness,
                bestTest,
                bestRuleCount,
                bestUsedRulesCount,
                avgFitness,
                avgTest,
                avgRuleCount,
                avgUsedRulesCount,
            )
            # Print the evolutionary log
            print(print_string, flush=True)

            # Log the evolution
            arch.saveEvo(log_string)

            # Elitism
            new_pop = [bestIndividual]

            # Save the best model
            data = []
            for rule in bestIndividual.rules:
                data.append([rule.pattern, rule.weight, rule.status])
            df = pd.DataFrame(data, columns=["pattern", "weight", "status"])
            arch.saveModel(df)

            # Select Parents (Tournament Selection) and Crossover
            for k in range(len(self.P.pop) - 1):
                tournament = []
                offspring = I.Individual(self.config)

                for j in range(self.tournamentSize):
                    tournament.append(self.P.pop[R.randint(0, len(self.P.pop) - 1)])

                tournament = self.sort_tournament(tournament)

                # We got two best parents
                parentA = copy.deepcopy(tournament[0])
                parentB = copy.deepcopy(tournament[1])

                # Do the crossover magic - Cluster crossover
                # Efficiency thing. Find the greater rule length
                lenA = len(parentA.rules)
                lenB = len(parentB.rules)
                maxLen = max(lenA, lenB)

                # we keep track of the rules we want to add to the offspring
                rules = []
                for j in range(maxLen):
                    if j < lenA:
                        ruleA = parentA.rules[j]
                        if ruleA.status == 1:
                            rules.append(ruleA)
                        elif (
                            R.random() < self.crossRate
                        ):  # we give unused rules some chance to get selected
                            rules.append(ruleA)
                    if j < lenB:
                        ruleB = parentB.rules[j]
                        if ruleB.status == 1:
                            rules.append(ruleB)
                        elif (
                            R.random() < self.crossRate
                        ):  # we give unused rules some chance to get selected
                            rules.append(ruleB)

                offspring.rules = rules
                offspring.bubbleSort()

                # Resize the offspring so it doesn't exceed the maximum allowed count
                while len(offspring.rules) > self.ruleCount:
                    countGreens = 0
                    for index in range(len(offspring.rules) - 1, -1, -1):
                        if countGreens >= index:
                            del offspring.rules[index]
                            break
                        else:
                            if offspring.rules[index].status == 0:
                                del offspring.rules[index]
                                break
                            else:
                                countGreens += 1

                new_pop.append(offspring)

            new_pop[0] = copy.deepcopy(new_pop[0])
            self.P.pop = new_pop

            # Mutations

            # We keep a copy of the elite
            elite = copy.deepcopy(self.P.pop[0])

            for indv in self.P.pop:
                needs_sort = False
                # On Model
                if R.random() <= self.mAR:
                    # add rule
                    self.mut_add_rule(indv)
                    needs_sort = True

                if R.random() <= self.mRR:
                    # remove rule
                    self.mut_remove_rule(indv)

                # On Rule
                for rule in indv.rules[:]:
                    if R.random() <= self.mCW:
                        # change weight
                        self.mut_change_weight(rule)
                    if R.random() <= self.mATP:
                        # add to pattern
                        self.mut_add_to_pattern(rule)
                        needs_sort = True
                    if R.random() <= self.mRFP:
                        # remove from pattern
                        self.mut_remove_from_pattern(rule)
                        needs_sort = True
                        if rule.pattern == "":
                            indv.rules.remove(rule)
                if needs_sort:
                    indv.bubbleSort()

            zeroFitness, testData = fitness.measureTotal(self.P.pop[0])

            # Check if elite got worse
            if elite.fitness < zeroFitness:
                self.P.pop[0] = elite

    def sort_tournament(self, t):
        t.sort(key=lambda x: x.fitness)

        if self.config["diversity_selection"] == "False":
            return t[:2]
        if self.config["diversity_selection"] != "True":
            raise ValueError(
                "Unknown value for diversity_selection: '{}'".format(
                    self.config["diversity_selection"]
                )
            )
        patterns = {rule.pattern for rule in t[0].rules}

        highest_closeness = 0
        highest_fitness = 0
        for i in range(1, len(t)):
            t[i].closeness = 0
            t[i].relative_fitness = 0
            for rule in t[i].rules:
                if rule.pattern in patterns:
                    t[i].closeness += 1
            if t[i].closeness > highest_closeness:
                highest_closeness = t[i].closeness
            if t[i].fitness > highest_fitness:
                highest_fitness = t[i].fitness

        if highest_closeness == 0:
            return t[:2]

        w = float(self.config["diversity_weight"])
        for i in range(1, len(t)):
            t[i].relative_fitness = (
                t[i].fitness / highest_fitness + w * t[i].closeness / highest_closeness
            )
        t2 = t[1:]
        t2.sort(key=lambda x: x.relative_fitness)
        return [t[0], t2[0]]

    # Add a random rule mutation
    def mut_add_rule(self, individual):
        if len(individual.rules) >= self.ruleCount:
            return
        pattern = ""
        weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
        # Add these many rules
        for i in range(R.randint(1, self.ruleSize)):
            # Rule size is calculated randomly, and now we need to select a random combination of codes with a
            # specified size
            randomchar = self.codes[R.randint(0, (len(self.codes) - 1))]
            pattern += randomchar
        rule = Rule.Rule(pattern, weight, 0)
        individual.rules.append(rule)

    # Remove rule mutation
    def mut_remove_rule(self, individual):
        if len(individual.rules) == 0:
            return
        tempRand = R.randint(0, len(individual.rules) - 1)
        del individual.rules[tempRand]

    # Add to weight mutation
    def mut_change_weight(self, rule):
        optRand = R.randint(0, 1)
        weightRand = round(R.uniform(self.mCWmin, self.mCWmax), 2)
        if optRand == 0:
            # Addition
            rule.weight += weightRand
        elif optRand == 1:
            # Substraction
            rule.weight -= weightRand

    # Alter patterns mutation (add letter)
    def mut_add_to_pattern(self, rule):
        if len(rule.pattern) >= self.ruleSize:
            return
        pattern = rule.pattern
        randomchar = self.codes[R.randint(0, len(self.codes) - 1)]
        if len(pattern) == 0:
            pattern = randomchar
        else:
            insPos = R.randint(0, len(pattern))
            pattern = pattern[0:insPos] + randomchar + pattern[insPos : (len(pattern))]
        rule.pattern = pattern

    # Alter patterns mutation (remove letter)
    def mut_remove_from_pattern(self, rule):
        if len(rule.pattern) == 0:
            return
        if len(rule.pattern) == 1:
            rule.pattern = ""
            return
        pattern = rule.pattern
        insPos = R.randint(0, len(pattern) - 1)
        pattern = pattern[0:insPos] + pattern[insPos + 1 : (len(pattern))]
        rule.pattern = pattern

    def removeExtra(self, indv):
        seen = {}
        to_remove = set()
        for i, rule in enumerate(indv.rules):
            if rule.pattern in seen:
                prev_idx = seen[rule.pattern]
                if rule.status == 0:
                    to_remove.add(i)
                elif indv.rules[prev_idx].status == 0:
                    to_remove.add(prev_idx)
                    seen[rule.pattern] = i
                else:
                    to_remove.add(i)
            else:
                seen[rule.pattern] = i
        indv.rules = [r for i, r in enumerate(indv.rules) if i not in to_remove]

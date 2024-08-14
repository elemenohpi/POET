# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import pop as Population
import fitness as F
import pandas as pd
import archivist as Archivist
import copy
import random as R
import individual as I
import rule as Rule
import time
from datetime import datetime


class Optimizer:
	def __init__(self, config, population,):
		self.config = config
		self.P = population
		self.runs = int(config["runs"])
		self.tournamentSize = int(config["tournament_size"])
		self.logInterval = int(config["pop_log_interval"])
		self.crossRate = float(config["crossover_unused_selection_chance"])
		self.ruleSize = int(config["maximum_rule_size"])
		self.minSize = int(config["minimum_rule_size"])
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
		pass

	def optimize(self):
		# We need an instance of the Fitness and Archivist classes for later usage
		fitness = F.Fitness(self.config)
		arch = Archivist.Archivist(self.config)

		# ToDo:: raise("Add dynamic mutation rates to escape from premature convergence.\nAdd Remove duplicate rules
		#  to make space for the new rules without changing the fitness")

		# growthLog = 0
		# mutRates = [self.mAR, self.mRR, self.mCW, self.mATP, self.mRFP, self.mCWmin, self.mCWmax]

		log_string = "best fitness, best test, best rule count, best unused rule count, average fitness, " \
					 "average test, average rule count, average unused rule count "
		arch.saveEvo(log_string)

		# For all the generations
		for i in range(self.runs):
			# To log the time
			start_time = int(round(time.time() * 1000))

			# Measure fitness
			avgFitness = 0.0
			avgTest = 0.0
			avgRuleCount = 0.0
			avgUsedRulesCount = 0.0

			# To calculate the best fitness
			# Temporarily set the first member as the best member

			self.P.pop[0].fitness, self.P.pop[0].test = fitness.measureTotal(self.P.pop[0])
			bestIndividual = self.P.pop[0]

			now = datetime.now()
			current_time = now.strftime("%H:%M:%S")
			# print("Current Time =", current_time)
			for j in self.P.pop:
				j.fitness, j.test = fitness.measureTotal(j)
				avgFitness += j.fitness
				avgTest += j.test
				avgRuleCount += len(j.rules)
				avgUsedRulesCount += j.usedRulesCount
				if j.fitness < bestIndividual.fitness:
					bestIndividual = j
			#
			# now = datetime.now()
			#
			# current_time = now.strftime("%H:%M:%S")

			# To calculate the average fitness and the average
			avgFitness = round(avgFitness / float(len(self.P.pop)), 3)
			avgTest = round(avgTest / float(len(self.P.pop)), 3)
			avgRuleCount = round(avgRuleCount / float(len(self.P.pop)), 0)
			avgUsedRulesCount = round(avgUsedRulesCount / float(len(self.P.pop)), 0)
			bestFitness = round(bestIndividual.fitness, 3)
			bestTest = round(bestIndividual.test, 3)
			bestRuleCount = len(bestIndividual.rules)
			bestUsedRulesCount = bestIndividual.usedRulesCount

			# Log the outcome before doing the changes to the population / generating a new population
			print_string = "{}: -b {} -bt {} -rc {} -urc {} ||| -a {} -at {} -arc {} -aurc {}".format(i, bestFitness,
																									  bestTest,
																									  bestRuleCount,
																									  bestUsedRulesCount,
																									  avgFitness,
																									  avgTest,
																									  avgRuleCount,
																									  avgUsedRulesCount)

			log_string = "{},{},{},{},{},{},{},{},{}".format(i, bestFitness,
															 bestTest,
															 bestRuleCount,
															 bestUsedRulesCount,
															 avgFitness, avgTest,
															 avgRuleCount,
															 avgUsedRulesCount)
			# Print the evolutionary log
			print(print_string, flush=True)

			# Log the evolution
			arch.saveEvo(log_string)

			# Create a copy of the population
			newPop = copy.deepcopy(self.P)
			newPop.pop.clear()

			# Elitism
			newPop.pop.append(bestIndividual)

			# Save the best model
			data = []
			for rule in bestIndividual.rules:
				data.append([rule.pattern, rule.weight, rule.status])
			df = pd.DataFrame(data, columns=['pattern', 'weight', 'status'])
			arch.saveModel(df)

			# Select Parents (Tournament Selection) and Crossover
			for k in range(len(self.P.pop) - 1):
				tournament = []
				offspring = I.Individual(self.config)

				for j in range(self.tournamentSize):
					# Randomly append as much individuals to the tournament as we need
					tournament.append(self.P.pop[R.randint(0, len(self.P.pop) - 1)])

				tournament = self.sort_tournament(tournament)

				# We got two best parents
				parentA = copy.deepcopy(tournament[0])
				parentB = copy.deepcopy(tournament[1])

				# Do the crossover magic - Cluster crossover
				# Efficiency thing. Find the greater rule length
				lenA = len(parentA.rules)
				lenB = len(parentB.rules)
				maxLen = lenA
				if lenA < lenB:
					maxLen = lenB

				# we keep track of the rules we want to add to the offspring
				rules = []
				for j in range(maxLen):
					if j < lenA:
						ruleA = parentA.rules[j]
						if ruleA.status == 1:
							rules.append(ruleA)
						elif R.random() < self.crossRate:  # we give unused rules some chance to get selected
							rules.append(ruleA)
					if j < lenB:
						ruleB = parentB.rules[j]
						if ruleB.status == 1:
							rules.append(ruleB)
						elif R.random() < self.crossRate:  # we give unused rules some chance to get selected
							rules.append(ruleB)

				offspring.rules = rules
				offspring.bubbleSort()

				# Resize the offspring so it doesn't exceed the maximum allowed count
				while len(offspring.rules) > self.ruleCount:
					countGreens = 0
					for index in range(len(offspring.rules) - 1, -1, -1):
						if (countGreens >= index):
							del (offspring.rules[index])
							break
						else:
							if (offspring.rules[index].status == 0):
								del (offspring.rules[index])
								break
							else:
								countGreens += 1

				newPop.pop.append(offspring)

			self.P.pop = copy.deepcopy(newPop.pop)

			# Mutations

			# We keep a copy of the elite
			elite = copy.deepcopy(self.P.pop[0])

			for indv in self.P.pop:
				# On Model
				if R.random() <= self.mAR:
					# add rule
					self.mut_add_rule(indv)

				if R.random() <= self.mRR:
					# remove rule
					self.mut_remove_rule(indv)

				# On Rule
				for rule in indv.rules:
					if R.random() <= self.mCW:
						# change weight
						self.mut_change_weight(rule)
					if R.random() <= self.mATP:
						# add to pattern
						self.mut_add_to_pattern(rule)
					if R.random() <= self.mRFP:
						# remove from pattern
						self.mut_remove_from_pattern(rule)
						if len(rule.pattern) <= self.minSize:
							indv.rules.remove(rule)
				indv.bubbleSort()

			zeroFitness, testData = fitness.measureTotal(self.P.pop[0])

			# Check if elite got worse
			if elite.fitness < zeroFitness:
				self.P.pop[0] = elite

			# self.P.pop[0].print()

			# for indv in self.P.pop:
			# 	# print(len(indv.rules))
			# 	self.removeExtra(indv)
			# 	# print("after: {}".format(len(indv.rules)))

			# self.P.pop[0].print()

			# To Log the time
			end_time = int(round(time.time() * 1000)) - start_time
		# print(end_time)
		pass

	def sort_tournament(self, t):
		n = len(t)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n - i - 1):
				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				if t[j].fitness > t[j + 1].fitness:
					t[j], t[j + 1] = t[j + 1], t[j]

		if self.config["diversity_selection"] == "False":
			return t
		if self.config["diversity_selection"] != "True":
			raise "unknown value for diversity_selection"
		patterns = [rule.pattern for rule in t[0].rules]

		highest_closeness = 0
		highest_fitness = 0
		for i in range(1, n):
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
			return t

		w = float(self.config["diversity_weight"])
		for i in range(1, n):
			t[i].relative_fitness = t[i].fitness / highest_fitness + w * t[i].closeness / highest_closeness
		t2 = t[1:n]
		t2.sort(key=lambda x: x.relative_fitness, reverse=False)
		return t[0], t2[0]

	# Add a random rule mutation
	def mut_add_rule(self, individual):
		if len(individual.rules) >= self.ruleCount:
			return
		pattern = ""
		weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
		# Add these many rules
		for i in range(R.randint(self.minSize, self.ruleSize)):
			# Rule size is calculated randomly, and now we need to select a random combination of codes with a
			# specified size
			randomchar = self.codes[R.randint(0, (len(self.codes) - 1))]
			pattern += randomchar
		rule = Rule.Rule(pattern, weight, 0)
		individual.rules.append(rule)

	# Remove rule mutation
	def mut_remove_rule(self, individual):
		# print("mrr")
		fitness = F.Fitness(self.config)
		if len(individual.rules) == 0:
			return
		tempRand = R.randint(0, len(individual.rules) - 1)
		del (individual.rules[tempRand])

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
		pass

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
			pattern = pattern[0:insPos] + randomchar + pattern[insPos:(len(pattern))]
		rule.pattern = pattern

	# Alter patterns mutation (remove letter)
	def mut_remove_from_pattern(self, rule):
		if len(rule.pattern) <= self.minSize:
			return
		pattern = rule.pattern
		insPos = R.randint(0, len(pattern) - 1)
		pattern = pattern[0:insPos] + pattern[insPos + 1:(len(pattern))]
		rule.pattern = pattern

	def removeExtra(self, indv):
		# removeList = []
		for j, rule in enumerate(indv.rules):
			for k in range(j + 1, len(indv.rules)):
				if indv.rules[k].pattern == indv.rules[j].pattern:
					# removeList.append(indv.rules[j])
					if indv.rules[j].status == 0:
						indv.rules.remove(indv.rules[j])
						j -= 1
					elif indv.rules[k].status == 0:
						indv.rules.remove(indv.rules[k])
						k -= 1
					break
				# for i in removeList:
# 	indv.rules.remove(i)

# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University
import rule as Rule
import random as R
import pandas as pd

import eletility

import regex



class Individual:
	# Constructor 
	def __init__(self, config):
		self.rules = [] # List of all rules/motifs/pattern
		self.usedRules = {}
		self.usedRulesCount = 0
		self.ruleSize = int(config["maximum_rule_size"]) # max number of characters in rules
		self.maxRuleCount = int(config["maximum_rule_count"])
		self.minWeight = float(config["rule_weight_min"])
		self.maxWeight = float(config["rule_weight_max"])
		self.fitness = 0
		self.test = 0
		self.extra = {}

		self.depth_tree = int(config["max_depth_tree"])
		self.min_braces = int(config["min_braces"])
		self.max_braces = int(config["max_braces"])



	def remove_unexpressed(self):
		to_be_removed = []
		for index, rule in enumerate(self.rules):
			if rule.status == 0:
				to_be_removed.append(rule)
		for rule in to_be_removed:
			self.rules.remove(rule)

	def makeFromFile(self, file):
		# print ("Creating an Individual from file...")
		self.rules = []
		self.usedRules = {}
		self.fitness = 0
		tempIndv = pd.read_csv(file)
		# print(tempIndv)
		tmpPatterns = tempIndv['pattern']
		tmpWeights = tempIndv['weight']
		try:
			tmpStatus = tempIndv['status']
		except:
			print("Pro-Predictor: model {} does not have status column, putting 0 for all.".format(file))
			tmpStatus = [0] * len(tmpPatterns)
		# print("size " + str(tmpPatterns.size))
		for i in range(0, len(tmpPatterns) - 1):
			# print (str(i) + " " + str(tmpPatterns[i]) + " " + str(tmpWeights[i]))
			rule = Rule.Rule(tmpPatterns[i], tmpWeights[i], tmpStatus[i])
			self.rules.append(rule)
		pass

	def init_pattern(self, init_method):
		# Load alphabet
		codes = pd.read_csv("data/translation/amino_to_amino.csv")
		codes = codes["code"].tolist()

		# Add random number of rules between 1 and maxRuleCount / 3
		for i in range(R.randint(1, int(self.maxRuleCount / 3))):
			# any value within the given interval [minWeight, maxWeight] (default 0-10) is equally likely to be drawn by uniform
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

			# Add these many rules according to the number of characters (ruleSize) to create the pattern
			pattern = ""
			for j in range(R.randint(1, self.ruleSize)):
				# Rule size is calculated randomly, and now we need to select a random combination of codes with a specified size
				randomchar = codes[R.randint(0, (len(codes) - 1))] # draw a random char
				pattern += randomchar

			# Create a new RE pattern
			if init_method == 'half':
				pattern_re = regex.indi_half(self.depth_tree, self.min_braces, self.max_braces)
			elif init_method == 'grow':
				pattern_re = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)
			elif init_method == 'full':
				pattern_re = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
			else:
				print('Error for initialisation method, choose between: grow, full or half')
				exit()

			rule = Rule.Rule(pattern_re, weight, 0)

			# Create a Rule instance
			# rule = Rule.Rule(pattern, weight, 0)
			self.rules.append(rule)

		self.bubbleSort() # Sort the array


	def init_formula(self):
		raise Exception("uncharted territories")
		pass

	def bubbleSort(self):
		# Sort the array Rules. The first element has the highest size and the last the smallest size
		n = len(self.rules)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n - i - 1):
				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				if len(self.rules[j].pattern) < len(self.rules[j + 1].pattern):
					self.rules[j], self.rules[j + 1] = self.rules[j + 1], self.rules[j]

	def print(self):
		for kh, rule in enumerate(self.rules):
			print(f"{kh}- {rule.pattern} - {rule.weight} - {rule.status}")


def main():
	configparser = eletility.ConfigParser()
	config = configparser.read("config.ini")

	new_indi = Individual(config)
	init_method = 'half'

	new_indi.init_pattern(init_method)
	
	new_indi.print()



if __name__ == '__main__':
	main()
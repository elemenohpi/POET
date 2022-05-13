# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University
import rule as Rule
import random as R
import pandas as pd


class Individual:
	# Constructor 
	def __init__(self, config):
		self.rules = []
		self.usedRules = {}
		self.usedRulesCount = 0
		self.ruleSize = int(config["maximum_rule_size"])
		self.maxRuleCount = int(config["maximum_rule_count"])
		self.minWeight = float(config["rule_weight_min"])
		self.maxWeight = float(config["rule_weight_max"])
		self.fitness = 0
		self.test = 0
		self.extra = {}

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

	def init_pattern(self):
		codes = pd.read_csv("data/translation/amino_to_amino.csv")
		codes = codes["code"].tolist()

		for i in range(R.randint(1, int(self.maxRuleCount / 3))):
			pattern = ""
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
			# Add these many rules
			for j in range(R.randint(1, self.ruleSize)):
				# Rule size is calculated randomly, and now we need to select a random combination of codes with a specified size
				randomchar = codes[R.randint(0, (len(codes) - 1))]
				pattern += randomchar
			rule = Rule.Rule(pattern, weight, 0)
			self.rules.append(rule)

		# for i in self.rules:
		# 	print(str(i.pattern) + " => " + str(i.weight))

		# print ("\n\n")
		self.bubbleSort()

	# for i in self.rules:
	# 	print(str(i.pattern) + " => " + str(i.weight))

	def init_formula(self):
		raise Exception("uncharted territories")
		pass

	def bubbleSort(self):
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
			print("{}- {} - {} - {}".format(kh, rule.pattern, rule.weight, rule.status))

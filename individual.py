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
		self.rules.sort(key=lambda x: len(x.pattern), reverse=True)

	def print(self):
		for kh, rule in enumerate(self.rules):
			print("{}- {} - {} - {}".format(kh, rule.pattern, rule.weight, rule.status))

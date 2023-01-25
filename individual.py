# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University
import rule as Rule
import random as R
import pandas as pd

import eletility

import regex
import re



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



	def check_pattern(self, pattern): # NS
		if pattern == None:
			return False
		else:
			# add all rules of the individual
			try:
				re.compile(pattern)
				return True
			except:
				return False


	def init_pattern(self, init_method):
		# Load alphabet
		codes = pd.read_csv("data/translation/amino_to_amino.csv")
		codes = codes["code"].tolist()

		# Add random number of rules between 1 and maxRuleCount / 3
		for i in range(R.randint(1, int(self.maxRuleCount / 3))):
			# any value within the given interval [minWeight, maxWeight] (default 0-10) is equally likely to be drawn by uniform
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

			# Add these many rules according to the number of characters (ruleSize) to create the pattern
			# pattern = ""
			# for j in range(R.randint(1, self.ruleSize)):
			# 	# Rule size is calculated randomly, and now we need to select a random combination of codes with a specified size
			# 	randomchar = codes[R.randint(0, (len(codes) - 1))] # draw a random char
			# 	pattern += randomchar




			# Create a new RE pattern
			if init_method == 'half':
				pattern_re, tree = regex.indi_half(self.depth_tree, self.min_braces, self.max_braces)
			
			elif init_method == 'grow':
				pattern_re, tree = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)
				# pattern_re, tree = regex.indi_grow(R.randint(2, self.depth_tree+1), self.min_braces, self.max_braces)

			elif init_method == 'full':
				pattern_re, tree = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
			else:
				print('Initialisation method Error, choose between: grow, full or half')
				exit()





			# rule = Rule.Rule(pattern, weight, 0) # Create a Rule instance
			# self.rules.append(rule)


			if self.check_pattern(pattern_re):
				rule = Rule.Rule(pattern_re, weight, 0, tree) # Create a Rule instance
				
				# checkpoint
				if rule.tree_shape[0] != 'cat' and rule.tree_shape[0] != '|':
					print('[ERROR] Root is incorrect:', rule.tree_shape)
					exit()

				self.rules.append(rule)

		self.bubbleSort() # Sort the array according to the size of motifs


	def init_formula(self):
		raise Exception("uncharted territories")
		pass

	def bubbleSort_old(self):
		# Sort the array Rules. The first element has the highest size and the last the smallest size
		n = len(self.rules)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n - i - 1):
				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				try:
					if len(self.rules[j].pattern) < len(self.rules[j + 1].pattern):
						self.rules[j], self.rules[j + 1] = self.rules[j + 1], self.rules[j]
				except TypeError:
					print('[ERROR] indi.py l143')

					print(self.rules[j].pattern)
					print(self.rules[j+1].pattern)

	# def get_leaf_nodes(self):

		# last_layer = 0
		# dict_layer = {}
		# # from last layer to the root (reverse direction )
		# for j, layer in enumerate(range(self.depth_tree-1, -1, -1)):
	
		# 	if j == 0:
		# 		last_layer = layer+1
	
		# 	# Count all nodes in the layer
		# 	for i in range (2**(self.depth_tree - (layer+1))):
		# 		node+=1

		# 		if (self.depth_tree-layer) not in dict_layer.keys():
		# 			dict_layer[self.depth_tree-layer] = []
		# 			dict_layer[self.depth_tree-layer].append(node-1)
		# 		else:
		# 			dict_layer[self.depth_tree-layer].append(node-1)

		# return dict_layer[last_layer]

	# def get_or_nodes(self, rule_tree):
	# 	counter = 0
	# 	for node in rule_tree:
	# 		if node == '|':
	# 			counter +=1
	# 	return counter

	def bubbleSort(self):
		self.compute_complexity()
		# Sort the array Rules. The first element has the highest size and the last the smallest size
		n = len(self.rules)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n - i - 1):
				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				try:
					if self.rules[j].complexity < self.rules[j+1].complexity:
						self.rules[j], self.rules[j + 1] = self.rules[j + 1], self.rules[j]
				except TypeError:
					print('[ERROR] indi.py l199')

					print(self.rules[j].pattern)
					print(self.rules[j+1].pattern)

	def compute_fitness(self):
		self.fitness = 0
		for rule in self.rules:
			self.fitness += rule.score

	def compute_complexity(self):
		'''
		on compte le nombre de feuille (sans les None) ainsi que le nombre de or |
		puis on soustrait ce nombre au nombre de feuille pour avoir la complexity
		'''

		# Create a list with all leaf nodes
		# leaves_list = self.get_leaf_nodes() 

		for rule in self.rules:

			nbr_leaf = 0
			nbr_or = 0

			for node in rule.tree_shape:
				if node in regex.LAST:
					nbr_leaf+=1
				if type(node) == list:
					nbr_leaf+=1
				if node == '|':
					nbr_or +=1

			rule.complexity = nbr_leaf - nbr_or

	def remove_double_rules(self):
		list_rule = []
		final = []

		for rule in self.rules:

			if rule.pattern not in list_rule:
				list_rule.append(rule.pattern)
				final.append(rule)
				
		self.rules = []

		for rule in final:
			self.rules.append(rule)


	def print(self):
		for kh, rule in enumerate(self.rules):
			print(f"{kh+1}- {rule.pattern} - {rule.weight} - {rule.status} - {rule.complexity} - {round(rule.score,2)}")
			# print(rule.tree_shape)
		print("Fitness = ",self.fitness)

# def main():
# 	configparser = eletility.ConfigParser()
# 	config = configparser.read("config.ini")

# 	new_indi = Individual(config)
# 	init_method = 'half'

# 	new_indi.init_pattern(init_method)
	
# 	new_indi.print()



# if __name__ == '__main__':
# 	main()
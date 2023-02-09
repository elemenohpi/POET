# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University
import rule as Rule
import random as R
import pandas as pd

import eletility

import regex
import re

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

import copy
import time

from fuzzywuzzy import fuzz



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
		self.max_identity = float(config['identity'])



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


	def init_regex_pattern(self, init_method):

		# Add random number of rules between 1 and maxRuleCount / 5
		for i in range(R.randint(1, int(self.maxRuleCount/4) )):

			# Create a ramdom weight (default 0-10)
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

			# Create a new regex pattern according to a tree construction method
			if init_method == 'half':
				pattern_re, tree = regex.indi_half(self.depth_tree, self.min_braces, self.max_braces)
			
			elif init_method == 'grow':
				pattern_re, tree = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)
				# pattern_re, tree = regex.indi_grow(R.randint(2, self.depth_tree+1), self.min_braces, self.max_braces)

			elif init_method == 'full':
				pattern_re, tree = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
			else:
				print('[ERROR] Initialization method Error, choose between: grow, full or half in config.ini file')
				exit()

			# checkpoint 1
			if self.check_pattern(pattern_re):
				rule = Rule.Rule(pattern_re, weight, 0, tree) # Create a Rule instance
				
				# checkpoint 2
				if rule.tree_shape[0] != 'cat' and rule.tree_shape[0] != '|':
					print('[ERROR] Root is incorrect:', rule.tree_shape)
					exit()

				# Add the new regex to the list of rules
				self.rules.append(rule)

		self.complexitySort() # Sort the array according to the complexity of regex


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

	def complexitySort(self):
		# self.compute_complexity()
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
					if self.rules[j].score < self.rules[j+1].score:
						self.rules[j], self.rules[j + 1] = self.rules[j + 1], self.rules[j]
				except TypeError:
					print('[ERROR] indi.py l176')

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




	def remove_double_rules_old(self):
		set_pattern = set()
		final = []

		# Remove duplicates
		for rule in self.rules:
			if rule.pattern not in set_pattern:
				set_pattern.add(rule.pattern)
				final.append(rule)

		self.rules.clear()  # Reset the list of rules
		self.rules = copy.copy(final) # replaces the list of rules by a list without duplicates

		afac=set()

		for i, rule1 in enumerate(self.rules):
			for j, rule2 in enumerate(self.rules):
				if i == j:
					pass
				else:

					# print('R1', rule1.pattern, rule1.weight, rule1.score)
					# print('R2',rule2.pattern, rule2.weight, rule2.score)

					if rule1.weight == rule2.weight:
					# try:
						ratio = fuzz.ratio(rule1.pattern, rule2.pattern)
		
						if ratio >= self.max_identity:
							# to_remove = j + i + 1 if self.rules[j + i + 1].score < rule.score else i
							# self.rules.remove(rule2.pattern)
							# self.rules.pop(j)

							sup_tuple = (i, j)
							sorted_ = tuple(sorted(sup_tuple))
							afac.add(sorted_)
			
		x = []
		for t in afac:
			# tt.remove(tt[t[1]])
			if self.rules[t[0]].score >= self.rules[t[1]].score:
				x.append(self.rules[t[1]])
			else:
				x.append(self.rules[t[0]])

		for i in x:
			self.rules.remove(i)









	def remove_double_rules_old_old(self):
		'''
		Remove duplicate regex and regex with more than 80% identity.
		'''

		set_pattern = []
		final = []

		nbr_ali = 0

		# Remove duplicates
		for rule in self.rules:
			if rule.pattern not in set_pattern:
				set_pattern.append(rule.pattern)
				final.append(rule)

		self.rules = [] # Reset the list of rules
		self.rules = copy.copy(final) # replaces the list of rules by a list without duplicates

		# supprime les regex.pattern qui sont trop similaire
		# on check le % d identite entre les regles qui ont le meme poids
		# car se sont probablement des regex trop proche du a une mutation
		for i, rule in enumerate(self.rules):
			for j, rule2 in enumerate(self.rules[i+1:]):
			# for j, rule2 in enumerate([r for r in self.rules[i+1:] if r.weight == rule.weight]):

				# if i == j:
				# 	pass
				# else:
				# try:
				# same weight ~ID
				if rule.weight == rule2.weight:
					p1 = Seq(rule.pattern)
					p2 = Seq(rule2.pattern)

					# alignement of the 2 patterns
					alignments = pairwise2.align.globalxx(p1, p2)
					best_alignment = alignments[0]
					identity = (best_alignment[2] / len(best_alignment[0])) * 100
					nbr_ali+=1

					# print('ALIGNEMENT prend', t3, 's')

					if identity >= self.max_identity: # defaut = 49.0
						# remove the rule with the lowest score
						to_remove = j + i + 1 if self.rules[j + i + 1].score < rule.score else i
						self.rules.pop(to_remove)
						
				else:
					pass


					# if self.rules[i].weight == self.rules[j].weight:
						# p1 = self.rules[i].pattern
						# p2 = self.rules[j].pattern
	
						# # alignement of the 2 patterns
						# alignments = pairwise2.align.globalxx(p1, p2)
						# best_alignment = alignments[0]
						# identity = (best_alignment[2] / len(best_alignment[0])) * 100
						
						# if identity >= self.max_identity: # defaut = 51.0
						# 	# remove the rule with the lowest score
						# 	self.rules.pop(j if self.rules[j].score < self.rules[i].score else i)

						# 	break
				# except:
				# 	pass


		# print('NBR ALIGNEMENT prend', nbr_ali)



			



	




	def print(self):
		for kh, rule in enumerate(self.rules):
			print(f"{kh+1}- {rule.pattern} - {rule.weight} - {rule.status} - {rule.complexity} - {round(rule.score,2)}")
			# print(rule.tree_shape)
		print("Fitness = ", round(self.fitness,3), '\n')

# def main():
# 	configparser = eletility.ConfigParser()
# 	config = configparser.read("config.ini")

# 	new_indi = Individual(config)
# 	init_method = 'half'

# 	new_indi.init_pattern(init_method)
	
# 	new_indi.print()



# if __name__ == '__main__':
# 	main()
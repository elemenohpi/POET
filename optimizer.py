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

from tqdm import tqdm
import regex
import re
import matplotlib.pyplot as plt
import numpy as np





class Optimizer:
	def __init__(self, config, population):
		self.config = config
		self.P = population
		self.runs = int(config["runs"]) # number of generations
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


		self.depth_tree = int(config["max_depth_tree"])
		self.min_braces = int(config["min_braces"])
		self.max_braces = int(config["max_braces"])
		self.maxnodes = (2**self.depth_tree)-1 
		self.dict_layer = {}
		self.last_layer = self.which_layer()
		self.mPP = float(config["mut_ponctual_point"])
		self.mAA = float(config["mut_add_aa"])
		self.min_new_aa = int(config["min_new_aa"])
		self.max_new_aa = int(config["max_new_aa"])




		# behaviour metrics
		self.mutation_count = 0
		self.mRR_count = 0

		self.list_mutation = []
		self.verbose = config["verbose"]
		self.pop_size = float(len(self.P.pop))


	def save_best_model(self, arch, bestIndividual):
		# Save the best model    --- Create function
		data = []
		for rule in bestIndividual.rules:
			data.append([rule.pattern, rule.weight, rule.status])
		df = pd.DataFrame(data, columns=['pattern', 'weight', 'status'])
		arch.saveModel(df)

	def optimize(self):
		# We need an instance of the Fitness and Archivist classes for later usage
		fitness = F.Fitness(self.config)
		arch = Archivist.Archivist(self.config)

		# ToDo:: raise("Add dynamic mutation rates to escape from premature convergence.
		# \nAdd Remove duplicate rules
		#  to make space for the new rules without changing the fitness")

		# growthLog = 0
		# mutRates = [self.mAR, self.mRR, self.mCW, self.mATP, self.mRFP, self.mCWmin, self.mCWmax]

		log_string = "best fitness, best test, best rule count, best unused rule count, average fitness, " \
					 "average test, average rule count, average unused rule count "
		arch.saveEvo(log_string)

		# loop in all generations
		for run in tqdm(range(self.runs)):
			self.mutation_count = 0 # to count the number of a type of mutation
			now = datetime.now()

			# start_time = int(round(time.time() * 1000))

			# Measure fitness
			avgFitness,avgRuleCount,avgUsedRulesCount  = 0.0, 0.0, 0.0
			# avgTest = 0.0 # p-value ??????????

			# Temporarily set the first member of the pop as the best member, to calculate the best fitness
			starter = self.P.pop[0]
			fitness.measureTotal(starter)
			bestIndividual = starter

			if self.verbose == 'True':
				print('--- Best individual ---')
				bestIndividual.print()




			###########################################
			###                                     ###
			###             EVALUATION              ###
			###                                     ###
			###########################################

			# Evaluation Step
			for i, individual in enumerate(self.P.pop):
				# Compute the fitness value of each individual
				fitness.measureTotal(individual)

				# Set the best individual
				if individual.fitness > bestIndividual.fitness:
					bestIndividual = individual
			
				avgFitness += individual.fitness
				avgRuleCount += len(individual.rules) # number of rules of the individual
				avgUsedRulesCount += individual.usedRulesCount # number of rules matched with dataset (status=1)

			# current_time = now.strftime("%H:%M:%S")

			# To calculate the average metrics
			avgFitness = round(avgFitness / self.pop_size, 3)
			avgRuleCount = round(avgRuleCount / self.pop_size, 0)
			avgUsedRulesCount = round(avgUsedRulesCount / self.pop_size, 0)
			
			# To calculate the best metrics
			bestFitness = round(bestIndividual.fitness, 3)
			bestRuleCount = len(bestIndividual.rules)
			bestUsedRulesCount = bestIndividual.usedRulesCount

			# Log the outcome before doing the changes to the population / generating a new population
			# print_string = f"{run}: -b {bestFitness} -bt {bestTest} -rc {bestRuleCount} -urc {bestUsedRulesCount} ||| -a {avgFitness} -at {avgTest} -arc {avgRuleCount} -aurc {avgUsedRulesCount}"
			print_string = f"Run {run}: -b {bestFitness} -rc {bestRuleCount} -urc {bestUsedRulesCount} ||| -a {avgFitness} -arc {avgRuleCount} -aurc {avgUsedRulesCount}"
			log_string = f"{run},{bestFitness},{bestRuleCount},{bestUsedRulesCount},{avgFitness},{avgRuleCount},{avgUsedRulesCount}"
			
			# Print the evolutionary log
			if self.verbose == 'True':
				print(print_string, flush=True)

			# Log the evolution
			arch.saveEvo(log_string)

			# Create a copy of the population = New pop
			newPop = copy.deepcopy(self.P)
			newPop.pop.clear() # Erase all individual of pop
			newPop.pop.append(bestIndividual) # Elitism

			# Save the best model    --- Create function
			self.save_best_model(arch, bestIndividual)
			# data = []
			# for rule in bestIndividual.rules:
			# 	data.append([rule.pattern, rule.weight, rule.status])
			# df = pd.DataFrame(data, columns=['pattern', 'weight', 'status'])
			# arch.saveModel(df)



			###########################################
			###                                     ###
			###              CROSSOVER              ###
			###                                     ###
			###########################################


			# Select Parents + Crossover
			for k in range(len(self.P.pop) - 1):

				# Create a new empty offspring
				offspring = I.Individual(self.config) 

				# Selection of 2 parents for crossover with Tournament method
				parentA, parentB = self.tournament()

				# print('---- Parent A ----')
				# parentA.print()
				# print('\n---- Parent B ----')
				# parentB.print()

				# Run crossover between 2 parents
				rules = self.crossover(parentA, parentB)
				offspring.rules = rules

				offspring.bubbleSort() # Sort rules according to their length

				# Resize the offspring so it doesn't exceed the maximum allowed count
				while len(offspring.rules) > self.ruleCount:
					del (offspring.rules[-1])

				# add new offspring
				newPop.pop.append(offspring)

			# Replace the population by the newPop containing offspring
			self.P.pop = copy.deepcopy(newPop.pop)

			# We keep a copy of the elite
			elite = copy.deepcopy(self.P.pop[0])


			###########################################
			###                                     ###
			###              MUTATIONS              ###
			###                                     ###
			###########################################

			for indv in self.P.pop:

				# -- On Model --

				# add a new rule/regex    - OK
				if R.random() <= self.mAR:
					self.MUT_add_rule(indv)

				# remove a rule/regex     - OK
				if R.random() <= self.mRR:
					self.MUT_remove_rule(indv)


				# # -- On Rule --
				for rule in indv.rules:

					# change weight of the rule 	- OK
					if R.random() <= self.mCW:
						self.MUT_change_weight(rule)

					# Replace to pattern - point mutation (modify a branch of the tree)     - OK
					if R.random() <= self.mATP:
						self.MUT_replace_subtree(rule)
					
					# Replace the value of one node (invert or replace) 	- ~~~
					if R.random() <= self.mPP:
						self.MUT_ponctual_point(rule)

					# remove from pattern      - dont work
					# if R.random() <= self.mRFP:
					# 	# pick a random node (avoid None node and root)
					# 	self.MUT_remove_from_pattern(rule)

					# Add X new AA in a leaf node       - OK
					if R.random() <= self.mAA:
						random_node = self.pick_a_node(rule.tree_shape)
						self.MUT_add_alphabet(rule, random_node)

			
				indv.bubbleSort()
				# print('---------')
				# indv.print()


			fitness.measureTotal(self.P.pop[0])
			self.list_mutation.append(self.mutation_count) # pour compter les mutations qui ont eu lieu

			# Check if elite got worse, and replace if
			if elite.fitness > self.P.pop[0].fitness:
				self.P.pop[0] = elite


			# print(len(self.P.pop))

			# for indv in self.P.pop:
			# 	# print(len(indv.rules))
			# 	self.removeExtra(indv)
			# 	# print("after: {}".format(len(indv.rules)))

			# self.P.pop[0].print()

			# To Log the time
			# end_time = int(round(time.time() * 1000)) - start_time
		# print(end_time)
		# plt.bar([1,2,3,4,5,6,7,8,9,10], self.list_mutation)
		# plt.show()
			

	# a crossover between 2 rules, not between 2 individuals # NS
	def internal_crossover(self):
		pass


	#  obtain 2 offsprings, not one # NS
	def OP_crossover_twice(self):
		pass


	# Obtain 1 offspring
	def crossover(self, parentA, parentB):
		'''
		Function to perform the crossover between parentA 
		and ParentB. Combine all the rules of ParentA and ParentB with status 1.
		Rules with status 0 have a probability of 0.2 to be selected

		ParentA/B: parents(individuals) used for the crossover

		return: an array containing all selected the rules for the offspring
		'''

		# Do the crossover magic - Cluster crossover
		# Efficiency thing. Find the greater rule length
		lenA = len(parentA.rules)
		lenB = len(parentB.rules)
		maxLen = lenA
		rules = []

		if lenA < lenB:
			maxLen = lenB

		# we keep track of the rules we want to add to the offspring
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

		return rules



	# [ATTENTION]
	# Supprime les elements au milieu, qu'importe leur taille 
	# prendre en compte les status 0 . PE parce aqu'on a pas de status 0
	# [ATTENTION]

	# # Resize the offspring so it doesn't exceed the maximum allowed count
	# # priorité à la suppression des règles courtes
	# while len(offspring.rules) > self.ruleCount:
	# 	countGreens = 0
	# 	for index in range(len(offspring.rules) - 1, -1, -1):
	# 		print(i)
	# 		if (countGreens >= index):
	# 			del (offspring.rules[index])
	# 			break
	# 		else:
	# 			if (offspring.rules[index].status == 0):
	# 				del (offspring.rules[index])
	# 				break
	# 			else:
	# 				countGreens += 1

	# self.P.pop = copy.deepcopy(newPop.pop)



	def tournament(self):
		'''
		The tournament operator. Select the 2 best individuals

		return : two best parents A and B
		'''
		tournament = []
		draw = []

		# Un individu ne peut jamais etre tire 2x
		while len(tournament) < self.tournamentSize:
			tmpIndi = R.randint(0, len(self.P.pop) - 1)
			
			# Avoid to draw the same individual twice in a tournament
			if tmpIndi not in draw:
				draw.append(tmpIndi)
				tournament.append(self.P.pop[tmpIndi])

		# Sort individuals according to their fitness value
		tournament = self.sort_tournament(tournament) 

		# We got two best parents of the tournaments
		parentA = copy.deepcopy(tournament[0])
		parentB = copy.deepcopy(tournament[1])

		return parentA, parentB

	def sort_tournament(self, tournament):
		'''
		Function to sort each individual in the tournament according to their fitness value
		tournament: Array with *self.tournamentSize* individuals

		return a sorting array, from the best to the worst individual
		'''
		n = len(tournament)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n - i - 1):

				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				if tournament[j].fitness > tournament[j + 1].fitness:
					tournament[j], tournament[j + 1] = tournament[j + 1], tournament[j]

		if self.config["diversity_selection"] == "False":
			return tournament
		

		if self.config["diversity_selection"] != "True":
			raise "[ERROR] Unknown value for diversity_selection: choose between True or False"
		

		# # other method to add diversity (not used with regex)

		# patterns = [rule.pattern for rule in tournament[0].rules]
		# highest_closeness = 0
		# highest_fitness = 0
		# w = float(self.config["diversity_weight"])


		# for i in range(1, n):
		# 	tournament[i].closeness = 0
		# 	tournament[i].relative_fitness = 0

		# 	# if an individual have a same pattern as the best individual of the tournament
		# 	for rule in tournament[i].rules:
		# 		if rule.pattern in patterns:
		# 			tournament[i].closeness += 1

		# 	if tournament[i].closeness > highest_closeness:
		# 		highest_closeness = tournament[i].closeness

		# 	if tournament[i].fitness > highest_fitness:
		# 		highest_fitness = tournament[i].fitness

		# # No individual is close to the best
		# if highest_closeness == 0:
		# 	return tournament

		# for i in range(1, n):
		# 	tournament[i].relative_fitness = tournament[i].fitness / highest_fitness + w * tournament[i].closeness / highest_closeness
		
		# t2 = tournament[1:n]
		# t2.sort(key=lambda x: x.relative_fitness, reverse=False)

		# return tournament[0], t2[0]


	# Add a random rule mutation
	def MUT_add_rule(self, individual):
		'''
		Mutation operator. Add a new rule to the list of rules of an individual

		individual: the target individual selected with a probability of *self.mAR*
		'''

		if len(individual.rules) >= self.ruleCount:
			return
		else:
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

			# Create a new RE pattern
			method = R.choice(['full','grow'])

			if method == 'full':
				pattern_re, tree = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
			elif method == 'grow':
				pattern_re, tree = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)

			rule = Rule.Rule(pattern_re, weight, 0, tree)
			individual.rules.append(rule)
	
	# Remove rule mutation
	def MUT_remove_rule(self, individual):
		'''
		Mutation operator. Remove a rule to the list of rules of an individual

		individual: the target individual selected with a probability of *self.mRR*
		'''

		if len(individual.rules) == 0:
			return
		else:
			tempRand = R.randint(0, len(individual.rules) - 1)
			del (individual.rules[tempRand])

	# Change weight mutation
	def MUT_change_weight(self, rule):
		'''
		Mutation operator. Modify the weight of a rule by addtion or substraction of a value

		rule: the target rule selected with a probability of *self.mCW*
		'''

		optRand = np.random.choice([0,1], p=[0.5,0.5])
		weightRand = round(R.uniform(self.mCWmin, self.mCWmax), 2)
		
		if optRand == 0: 			# Addition
			rule.weight += weightRand
		elif optRand == 1: 			# Substraction
			rule.weight -= weightRand

	# Replace patterns mutation (replace sub-tree)
	def MUT_replace_subtree(self, rule):
		# Draw a random node
		random_node = self.pick_a_node(rule.tree_shape)

		# Define the layer of the random node
		for key_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = key_layer

		# Leaf case, node in the last layer
		if layer == self.last_layer:
			# Change the leaf and modify tree shape and pattern of the rule
			# if the leaf is a list
			if type(rule.tree_shape[random_node]) == list:
				if ('^' in rule.tree_shape[random_node]): # [^] case
					rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint(1, len(regex.ALPHABET)-1))
					rule.tree_shape[random_node].insert(0,'^')
				else: # [] case
					rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint(1, len(regex.ALPHABET)-1))
			# if the leaf is a single element
			else:
				rule.tree_shape[random_node] = R.choice(regex.LAST)

		else:
			# Generate a new sub-tree
			tree_method = np.random.choice(['full', 'grow'], p=[0.5,0.5])

			if tree_method == 'full':
				new_regex, new_tree = regex.indi_full( (self.depth_tree - layer)+1, self.min_braces, self.max_braces) 
			elif tree_method == 'grow':
				new_regex, new_tree = regex.indi_grow( (self.depth_tree - layer)+1, self.min_braces, self.max_braces) 

			rule.tree_shape[random_node] = new_tree[0]
	
			child = [random_node]
			i=0

			# Add all nodes of the new subtree
			while (i*2)+1 != self.maxnodes: 
				if i in child:
					child.append((i*2)+1) # left
					child.append((i*2)+2) # right
				i+=1

			# Replace nodes
			for i, node_index in enumerate(child):
				try:
					rule.tree_shape[node_index] = new_tree[i]
				except IndexError:
					rule.tree_shape[node_index] = None

		# New regex pattern
		rule.pattern = regex.tree2regex(rule.tree_shape)

	# Change/invert one node
	def MUT_ponctual_point(self, rule):

		# pick a random node (avoid None node)
		random_node = self.pick_a_node(rule.tree_shape)
		
		if rule.tree_shape[random_node] == 'cat':
			rule.tree_shape[random_node] = '|'

		elif rule.tree_shape[random_node] == '|':
			rule.tree_shape[random_node] = 'cat'

		elif rule.tree_shape[random_node] == '+':
			rule.tree_shape[random_node] = '{' + str(R.randint(self.min_braces, self.max_braces)) + '}'
		
		elif '{' in rule.tree_shape[random_node]:
			rule.tree_shape[random_node] = '{' + str(R.randint(self.min_braces, self.max_braces)) + '}'
		
		elif type(rule.tree_shape[random_node]) == list:
			if ('^' in rule.tree_shape[random_node]):
				rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint(1, len(regex.ALPHABET)-1))
				rule.tree_shape[random_node].insert(0,'^')
			else:
				rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint(1, len(regex.ALPHABET)-1))

		elif rule.tree_shape[random_node] == '[]':
			try:
				rule.tree_shape[random_node] = '[^]'
				rule.tree_shape[(random_node*2) + 1].insert(0,'^')
			except AttributeError:
				print('[ERROR]')
				print(rule.tree_shape,rule.tree_shape[random_node], random_node)
				print(rule.pattern)

		elif rule.tree_shape[random_node] == '[^]':
			rule.tree_shape[random_node] = '[]'
			if ('^' in rule.tree_shape[(random_node*2) + 1]):
				rule.tree_shape[(random_node*2) + 1].remove('^')

		elif rule.tree_shape[random_node] in regex.LAST:
			rule.tree_shape[random_node] = R.choice(regex.LAST)

		rule.pattern = regex.tree2regex(rule.tree_shape)


	# Alter patterns mutation (remove letter)    ------> DONT WORK
	def MUT_remove_from_pattern(self, rule):
		copysave = copy.deepcopy(rule.tree_shape)
		# print(">>> AVANT",rule.tree_shape)
		# print(">>> AVANT",rule.pattern)
		random_node = self.pick_a_node(rule.tree_shape)
		# print('------------------------------------------------------')
		# print('Noeud supprime', random_node, rule.tree_shape[random_node])
		# self.deleteTree(random_node, rule)


		# Define the layer of the random node
		for key_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = key_layer

		# print('layer=',layer)
		# Define last layer, and leaf nodes
		# right_subtree =  self.create_subtree_list(1)
		# left_subtree = self.create_subtree_list(2)
		# right_leaf = 0
		# left_leaf = 0

		# for i in right_subtree:
		# 	if rule.tree_shape[i] == None:
		# 		pass
		# 	else:
		# 		right_leaf = i

		# for i in left_subtree:
		# 	if rule.tree_shape[i] == None:
		# 		pass
		# 	else:
		# 		left_leaf = i

		if layer == self.last_layer: # Never delete a leaf node
			pass
		# elif random_node == right_leaf: # Never delete a leaf node
		# 	pass
		# elif random_node == left_leaf: # Never delete a leaf node
		# 	pass
		else:
			self.mutation_count+=1

			# Define list of node in child subtree
			child = self.create_subtree_list(random_node)
			# print('j influence', child)

			# Replace all child nodes by None value
			for i, node_index in enumerate(child):
				rule.tree_shape[node_index] = None

			# Define parent and brother of the random node
			parent = self.I_am_ur_parent(random_node)
			brother = self.I_am_ur_brother(random_node)

			# print('parent',parent, rule.tree_shape[parent])
			# print('brother',brother, rule.tree_shape[brother])

			# Define list of node in brother and parent subtree
			brother_tree = self.create_subtree_list(brother)
			parent_tree = self.create_subtree_list(parent)
			# print('parent_tree',parent_tree )
			# print('brother_tree',brother_tree)

			# Replace parent nodes by brother nodes
			for i, j in enumerate(parent_tree):
				try:
					rule.tree_shape[parent_tree[i]] = rule.tree_shape[brother_tree[i]]
				except IndexError:
					rule.tree_shape[parent_tree[i]] = None

		
		# Checkpoint
		if rule.tree_shape[0] != 'cat' and rule.tree_shape[0] != '|':
			rule.tree_shape[1] = rule.tree_shape[0]
			rule.tree_shape[0] = 'cat'

		try:
			rule.pattern = regex.tree2regex(rule.tree_shape)
		except:
			# Error of compilation, tree shape is erroneous. remove changing
			rule.tree_shape = copysave
			rule.pattern = regex.tree2regex(rule.tree_shape)


		

	# Add some new AA in leaf node
	def MUT_add_alphabet(self, rule, random_node):
		# Define the layer of the random node
		for key_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = key_layer

		if layer == self.last_layer: # Modify only leaf node (no list)
			if type(rule.tree_shape[random_node]) != list:
				
				nbr_new_aa = R.randint(self.min_new_aa, self.max_new_aa) # Number of AA to add
				for i in range(nbr_new_aa):
					aa = R.choice(regex.LAST)
					rule.tree_shape[random_node] += aa

				rule.pattern = regex.tree2regex(rule.tree_shape)


		


	def create_subtree_list(self, random_node):
		# Create list of node in list_nodes subtree
		list_nodes = [random_node]
		i=0

		# Add all nodes of the new subtree
		while (i*2)+1 != self.maxnodes: 
			if i in list_nodes:
				list_nodes.append((i*2)+1) # left
				list_nodes.append((i*2)+2) # right
			i+=1

		return list_nodes

	def I_am_ur_parent(self, indexNode):
		if indexNode == 0:
			return 0
		else:
			if indexNode%2 == 0:
				return int((indexNode-2)/2)
			else:
				return int((indexNode-1)/2)

	def I_am_ur_brother(self, indexNode):
		if indexNode == 0:
			return 0
		else:
			if indexNode%2 == 0:
				return int(indexNode-1)
			else:
				return int(indexNode+1)

	def pick_a_node(self, tree_shape):
			'''
			Pick a random node (avoid None node)
			tree_shape: the tree in which a random node is drawn
	
			return a integer corresponding to the index of a random node
			'''
			if len(tree_shape) == 0:
				print('[ERROR] len = 0')
				exit()
			if len(tree_shape) == None:
				print('[ERROR] len = None')
				exit()
			if tree_shape[0] != 'cat' and tree_shape[0] != '|':
				print('[ERROR] Root error')
				print('>>>',tree_shape)
				exit()
			
	
			status = False
			while status == False:
				# print('je suis bloque ici')
				# print(tree_shape)

				random_node = R.randint(0, len(tree_shape) - 1)
	
				if tree_shape[random_node] == None: # don't modify None node
					status = False
				elif random_node == 0: # don't touch the root
					status = False
				else: 
					status = True
	
			return random_node

	# def remove_or(self, rule, node):
	# 	'''
	# 	manages the parent nodes, because in some cases there can be errors, 
	# 	for example with the or '|'

	# 	rule: the rarget rule
	# 	node: the child node
	# 	'''

	# 	parent = self.I_am_ur_parent(node)

	# 	if rule.tree_shape[node] == None:
	# 		pass

	# 	elif rule.tree_shape[parent] == '[]':
	# 		self.remove_or(rule, parent)
	# 		rule.tree_shape[parent] = None
	# 		rule.tree_shape[node] = None

	# 	elif rule.tree_shape[parent] == '[^]':
	# 		self.remove_or(rule, parent)
	# 		rule.tree_shape[parent] = None
	# 		rule.tree_shape[node] = None

	# 	elif rule.tree_shape[parent] == '+':
	# 		self.remove_or(rule, parent)
	# 		rule.tree_shape[parent] = None
	# 		rule.tree_shape[node] = None

	# 	elif '}' in rule.tree_shape[parent]:
	# 		self.remove_or(rule, parent)
	# 		rule.tree_shape[parent] = None
	# 		rule.tree_shape[node] = None

	# 	elif rule.tree_shape[parent] == '|':
	# 		rule.tree_shape[parent] = 'cat'
	# 		rule.tree_shape[node] = None
	# 	else:
	# 		pass

	def which_layer(self):
	
		node = 0
	
		# from last layer to the root (reverse direction )
		for j, layer in enumerate(range(self.depth_tree-1, -1, -1)):
	
			if j == 0:
				self.last_layer = layer+1
	
			# Count all nodes in the layer
			for i in range (2**(self.depth_tree - (layer+1))):
				node+=1

				if (self.depth_tree-layer) not in self.dict_layer.keys():
					self.dict_layer[self.depth_tree-layer] = []
					self.dict_layer[self.depth_tree-layer].append(node-1)
				else:
					self.dict_layer[self.depth_tree-layer].append(node-1)

		return self.last_layer

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

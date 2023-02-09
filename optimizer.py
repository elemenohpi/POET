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

import motif_bank
from statistics import mean
import evaluation as Evaluation

import seaborn as sns
import scipy as sp

import multiprocessing

from Bio import pairwise2
from Bio.Seq import Seq





class Optimizer:

	def __init__(self, config, population):
		self.config = config 									# Config parameters
		self.P = population 									# the initial population
		self.runs = int(config["runs"]) 						# number of generations
		self.tournamentSize = int(config["tournament_size"])	# Size of the tournament during selection step
		self.reduction_tournament_size = int(config["red_tournament_size"]) #  Size of the tournament during reduction step
		# self.logInterval = int(config["pop_log_interval"])
		# self.crossRate = float(config["crossover_unused_selection_chance"])
		# self.ruleSize = int(config["maximum_rule_size"])
		self.ruleCount = int(config["maximum_rule_count"])		# Max number of regex/rules
		self.minWeight = float(config["rule_weight_min"])		# a remplacer par un genre d'identifiant
		self.maxWeight = float(config["rule_weight_max"])		# a remplacer par un genre d'identifiant
		# self.output_evo = config["output_evo"]
		# self.output_model = config["output_model"]

		# Mutation rates
		self.mAR = float(config["mut_add_rule"]) 				# Add Rule
		self.mRR = float(config["mut_remove_rule"]) 			# Remove Rule
		self.mReR = float(config["mut_replace_rule"]) 			# Replace Rule
		self.mCW = float(config["mut_change_weight"]) 			# Change Weight       # INUTILE !!!!
		self.mReS = float(config["mut_replace_subtree"]) 		# Replace Subtree
		self.mRFP = float(config["mut_remove_from_pattern"]) 	# Remove subtree
		self.mRN = float(config["mut_point"]) 					# Replace Node
		self.mAA = float(config["mut_add_aa"]) 					# Add Alphabet

		self.mCWmin = 0				# AFAC
		self.mCWmax = 1				# AFAC
		# codes = pd.read_csv("data/translation/amino_to_amino.csv")
		# self.codes = codes["code"].tolist()

		# Regex Parameters
		self.seed = int(config["seed"])							# Seed value
		self.depth_tree = int(config["max_depth_tree"])			# Max depth of tree shape
		self.min_braces = int(config["min_braces"])				# Min value in {}
		self.max_braces = int(config["max_braces"])				# Max value in {}
		self.verbose = config["verbose"] 						# Display some information
		self.pop_size = int(config["population_size"])			# Size of the initial population
		self.init_method = str(config["init_pop_met"])			# Individual construction method
		self.alphabet = config["learn_data"]					# Name of the alphabet used for the encoded data
		self.MAX_motif_size = int(config["max_motif_size_db"])	# Max size of motif identified by regex
		self.MIN_motif_size = int(config["min_motif_size_db"])	# Min size of motif identified by regex
		self.db_strat = config["db_motif"]						# Construction method of the motif database
		self.insertion_mtd = config["insertion_pop"]			# Method for generating the n+1 population
		self.reduction_method = config['reduction_method']		# Algorithm used during the reduction step

		self.maxnodes = (2**self.depth_tree)-1 					# Max number of nodes in the tree 
		self.dict_layer = {}									# Dictionnary with key= #layer and value=nodes in the layer
		self.last_layer = self.which_layer()					# Num of the last layer of the tree 
		self.stop_fitness = []									# Checkpoint to stop evolution

		# Mutation counters
		self.f_add_rule = []									#List with number of mutation Add Rule per run
		self.f_remove_rule = []									#List with number of mutation Remove Rule per run
		self.f_replace_rule = []								#List with number of mutation Replace Rule per run
		self.f_replace_subtree = []								#List with number of mutation Replace Subtree per run
		self.f_add_aa = []										#List with number of mutation Add Alphabet element per run
		self.f_replace_node = []								#List with number of mutation Replace Node per run
		self.f_remove_subtree = []								#List with number of mutation Remove Subtree per run



		# self.training_db_motif = {} # used only with the fitness function 'motif' !!!! enlever et mis dans fitness




	# def build_db_motif(self):
	# 	for i in range (self.MIN_motif_size, self.MAX_motif_size+1):
	# 		# Strategies
	# 		if self.db_strat == 'std': 	# Best CEST value
	# 			self.training_db_motif = motif_bank.standard_motif_bank(i, self.alphabet)
	# 		if self.db_strat == 'avg':	# Average CEST value 
	# 			self.training_db_motif = motif_bank.average_motif_bank(i, self.alphabet)
	# 		if self.db_strat == 'occ':	# Number of motif in each class define the classe
	# 			self.training_db_motif = motif_bank.occ_motif_bank(i, self.alphabet)

	

	def save_best_model(self, archive_instance, bestIndividual, verbose):
		''' Save the best model (each regex/rule)

		archive_instance: instance of Archivist
		bestIndividual: best individual/model to save
		verbose: Boolean

		'''
		data = []

		for regex in bestIndividual.rules:
			data.append([regex.weight, regex.pattern, regex.score])
		
		df = pd.DataFrame(data, columns=['weight','pattern', 'score'])
		archive_instance.saveModel(df, verbose)

	def stop_condition(self, stop=5):
		'''
		Interrupt the evolution process if the fitness of the best individual doesn't change anymore during {stop} runs.
		stop:  Number of runs before the process stops

		return a boolean, if True the process stops
		'''

		if len(self.stop_fitness) < stop:
			return False

		for i in range(len(self.stop_fitness) - stop, len(self.stop_fitness)):
			if self.stop_fitness[i] != self.stop_fitness[i-1]:
				return False
		return True










	def remove_double_rules(self, individual):
		'''
		Remove duplicate regex and regex with more than 80% identity.
		'''

		set_pattern = []
		final = []

		nbr_ali = 0

		# Remove duplicates
		for rule in individual.rules:
			if rule.pattern not in set_pattern:
				set_pattern.append(rule.pattern)
				final.append(rule)

		individual.rules = [] # Reset the list of rules
		individual.rules = copy.copy(final) # replaces the list of rules by a list without duplicates

		# supprime les regex.pattern qui sont trop similaire
		# on check le % d identite entre les regles qui ont le meme poids
		# car se sont probablement des regex trop proche du a une mutation
		for i, rule in enumerate(individual.rules):
			for j, rule2 in enumerate(individual.rules[i+1:]):
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

					if identity >= 49.0: # defaut = 49.0
						# remove the rule with the lowest score
						to_remove = j + i + 1 if individual.rules[j + i + 1].score < rule.score else i
						individual.rules.pop(to_remove)
						break
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
		return individual













	def optimize(self):
		'''
		Main function of the evolutionary process
		'''


		###########################################
		###                                     ###
		###            CONFIGURATION            ###
		###                                     ###
		###########################################

		R.seed(self.seed) # activate the seed
		print('[INFO] Seed:', self.seed)

		# Initialisation of instance of the Fitness/Archivist/Evaluation classes for later usage
		fitness_instance = F.Fitness(self.config)
		archive_instance = Archivist.Archivist(self.config)
		evaluate_instance = Evaluation.Eval(self.config)
		individual_instance = I.Individual(self.config)

		# Metrics
		dev_fitness_count = [] 	# List: std. dev. of fitness values of all individuals /run
		dev_rule_count = []    	# List: std. dev. of number of rule/regex of all individuals /run
		population_size = []	# List: population size /run


		# save the evolution information in a log
		log_string = "Run;best_fitness;best_rule_count;average_fitness;" \
					 "average_rule_count;Population_size;Time;Eval_time;" \
					 "Xover_time;Mut_time;mAR;mRR;mReR;mCW;mReS;mRFP;mRN;mAA"
		archive_instance.saveEvo(log_string)



		###########################################
		###                                     ###
		###             EVOLUTION               ###
		###                                     ###
		###########################################

		if self.verbose == 'True':
			print('[INFO] Running the Evolution process')

		# 1st evaluation of each individual + initialisation of the best individual
		x = time.time()
		bestIndividual = self.first_evaluation(fitness_instance, self.P.pop[0], self.verbose)

		# print('BEST INDI a la 1ere EVAL')
		# bestIndividual.print()
		# print('---------------')

		y = time.time()
		z = round(y-x, 2)
		# print('------- la 1ere EVAL prend', z, 's')


		# Main evolution loops
		for RUN in tqdm(range(self.runs)):
			start_time = time.time()

			# Initialization - Metrics Avg. fitness/Rules
			avgFitness = [] 	# Average Fitness of all individuals
			avgRuleCount = [] 	# Average number of rule/regex

			if self.stop_condition(stop=20):
				print(f'[WARNING] the evolutionary process was interrupted at run {RUN} ' \
					  f'because the fitness was not increasing anymore')
				break

			elite = copy.deepcopy(bestIndividual) 	# We keep a copy of the elite before crossover and mutation
			
			# print('ELITE au debut du RUN')
			# elite.print()
			# print('---------------')

			# bestIndividual = max( self.P.pop, key=lambda x: x.fitness)
			# print('BEST INDI au debut du RUN')
			# bestIndividual.print()
			# print('---------------')


			# Sort the population according to the fitness value and get the best individual
			# bestIndividual = max( self.P.pop, key=lambda x: x.fitness)
			elite = copy.deepcopy(bestIndividual) 	# We keep a copy of the elite before crossover and mutation



			# Create a new empty population to store the futur offspring
			newPop = Population.Population(self.config, empty=True)


			###########################################
			###                                     ###
			###              CROSSOVER              ###
			###                                     ###
			###########################################	
			

			x = time.time()

			# crossover step
			p1 = 0
			p2 = 0
			p3 = 0
			p4 = 0

			for _ in range(len(self.P.pop)): 

				# self.P.pop[_].print()

				a1,a2,a3,a4 = self.run_crossover(I, newPop)
				p1+=a1
				p2+=a2
				p3+=a3
				p4+=a4

			# print('------- les OFFSPRING prennent', p1,p2,p3,p4, 's')


			y = time.time()
			crossover_time = round(y-x, 2)
			# print('------- le CROSSOVER prend', crossover_time, 's')













			# # Replace the population by the newPop containing only the offspring
			# if self.insertion_mtd == 'generational':
			# 	self.P.pop.clear() # Erase all individual of pop
			# 	self.P.pop = copy.deepcopy(newPop.pop)

			# # Add offspring in the previous generation to mix parent + offspring
			# elif self.insertion_mtd == 'mix':
			# 	self.P.pop += newPop.pop
			# else:
			# 	print('[ERROR] The population insertion method is wrong, '\
			# 		  'choose between generational or mix, in the config.ini file')
			# 	exit()


			###########################################
			###                                     ###
			###              MUTATIONS              ###
			###                                     ###
			###########################################

			dico_cnt_mut = {'cnt_add_rule':0,
							'cnt_remove_rule':0,
							'cnt_replace_rule':0,
							'cnt_replace_subtree':0,
							'cnt_add_aa':0,
							'cnt_replace_node':0,
							'cnt_remove_subtree':0}




			x = time.time()

			# for indv in self.P.pop: # les mutations sont effectuees sur les parents + offspring uniquement
			for indv in newPop.pop: # les mutations sont effectuees sur les offspring uniquement
				# R.seed(self.seed)

				self.run_mutations(indv, dico_cnt_mut) 	# actuellement chaque individu a une proba p de subir chaque mutation
														# on peut aussi imaginer qu'un individu a une proba p de subir qu'une seule mutation
				indv.remove_double_rules_old() # Duplicates are removed + the rules too similar
				indv.complexitySort() # Sort individual's rules according to their complexity

			y = time.time()
			mut = round(y-x, 2)
			# print('------- Etape 1 -  mutation + remove prend', mut, 's')

			# x = time.time()

			# pool = multiprocessing.Pool(process=100)

			# indis2 = pool.map(self.remove_double_rules, newPop.pop)

			# # newPop.pop.clear() # Erase all individual of pop
			# # newPop.pop = copy.deepcopy(indis2)

			# pool.close()
			# pool.join()

			# y = time.time()
			# mut = round(y-x, 2)
			# print('------- Etape 2 -  mutation prend', mut, 's')


			# x = time.time()

			# for indv in newPop.pop: # les mutations sont effectuees sur les offspring uniquement
			# 	indv.remove_double_rules_old() # Duplicates are removed + the rules too similar
			# 	indv.complexitySort() # Sort individual's rules according to their complexity





				# Verifier ici quon supprime bien les doubles parce quil y a un bug


	
			# y = time.time()
			# mut = round(y-x, 2)
			# print('------- Etape 3 -  mutation prend', mut, 's')
			# for indv in newPop.pop: # les mutations sont effectuees sur les offspring uniquement

			# 	# indv.remove_double_rules() # Duplicates are removed + the rules too similar
				
			# 	indv.complexitySort() # Sort individual's rules according to their complexity




				# a commenter si //
				# indv.fitness = fitness.evaluate_individual(indv)


			# recalculate the fitness value because operators can impact the individual //

			start_eval_time = time.time()
			pool2 = multiprocessing.Pool()

			# renvoie une liste de nouveaux individus
			tprPop = pool2.map(fitness_instance.evaluate_individual, newPop.pop)
			
			pool2.close()
			pool2.join()

			stop_eval_time = time.time()
			eval_time = round((stop_eval_time - start_eval_time), 2)
			# print(f'----- Evaluation step took {eval_time} sec.')





			# for j, i in enumerate(indis2):
			# 	i.print()
			# 	print(j, '+++')
			# for j, i in enumerate(newPop.pop):
			# 	i.print()
			# 	print(j, '---')

			# for j, i in enumerate(indis):
			# 	i.print()
			# 	print(j, '***')



			# bestIndividual = max(indis, key=lambda x: x.fitness)

			# print('bestIndividual apre Xover et mutation')
			# bestIndividual.print()
			# print('---------------')


			# on selectionne le nouveau meilleur individu
			bestIndividual = copy.deepcopy(max(tprPop, key=lambda x: x.fitness))
			




			# bestIndividual.remove_double_rules_old()
			# bestIndividual.compute_fitness()

			# print('Nouveau bestIndividual apre Xover mutation remove double et selection')
			# bestIndividual.print()
			# print('---------------')



			# print('Elite apre Xover et mutation')
			# elite.print()
			# print('---------------')


			if self.config["fitness_alg"] == "motif":
				if bestIndividual.fitness < elite.fitness:
					# bestIndividual = copy.deepcopy(bestIndividual)
					bestIndividual = copy.deepcopy(elite)
					# print('Elite si Elite est meilleur')
					# elite.print()
					# print('---------------')

					# print('bestIndividual si Elite est meilleur')
					# bestIndividual.print()
					# print('---------------')



			x = time.time()
			# newPop.pop.clear() # Erase all individual of pop
			# newPop.pop = copy.deepcopy(tprPop)
			# newPop.pop.append(elite) # Elitism

			tprPop.append(elite)

			y = time.time()
			z = round((y - x), 2)
			# print(f'------- Copy de la pop prend {z} sec.')





			# print('Elite FINAL')
			# elite.print()
			# print('---------------')




			# for i in newPop.pop:
			# 	i.print()



			# print(len(newPop.pop), '-------------------------------------------')
			# for i in newPop.pop:
			# 	i.print()
			# exit()


			# # the indi #1 (elite) could be affected by a mutation and his fitness could decrease
			# # So, the fitness is calculated again and compared to the previous elite value
			# if self.config["fitness_alg"] == "RMSE2":
			# 	if bestIndividual.fitness < elite.fitness:
			# 		elite = copy.deepcopy(bestIndividual)

			# if self.config["fitness_alg"] == "motif":
			# 	if indv.fitness > elite.fitness:
			# 		bestIndividual = copy.deepcopy(indv)
			# 		elite = copy.deepcopy(indv)


			
			

			# # the indi #1 (elite) could be affected by a mutation and his fitness could decrease
			# # So, the fitness is calculated again and compared to the previous elite value
			# if self.config["fitness_alg"] == "RMSE2":
			# 	if bestIndividual.fitness < elite.fitness:
			# 		elite = copy.deepcopy(bestIndividual)

			# if self.config["fitness_alg"] == "motif":
			# 	if best.fitness > elite.fitness:
			# 		bestIndividual = copy.deepcopy(best)
			# 		elite = copy.deepcopy(best)


		




			x = time.time()


			# replace the population by the newPop containing only the offspring
			if self.insertion_mtd == 'generational':
				self.P.pop.clear() # Erase all individual of pop
				self.P.pop = copy.deepcopy(tprPop)
				self.P.pop.append(elite) # Elitism

			# Add offspring in the previous generation to mix parent + offspring
			elif self.insertion_mtd == 'mix':
				self.P.pop += tprPop
			else:
				print('[ERROR] The population insertion method is wrong, '\
					  'choose between generational or mix, in the config.ini file')
				exit()

			y = time.time()
			z = round((y - x), 2)
			# print(f'------- ajoute de pop parent + enfant prend {z} sec.')



			###########################################
			###                                     ###
			###             REDUCTION               ###
			###                                     ###
			###########################################

			x = time.time()

			# Create a copy of the current population => FinalPop
			FinalPop = Population.Population(self.config, empty=True)



			# FinalPop = copy.deepcopy(self.P)
			# FinalPop.pop.clear() # Erase all individual of pop
			FinalPop.pop.append(elite) # Elitism
			FinalPop.pop.append(bestIndividual) # Elitism
			

			# [RANK method]
			if self.reduction_method == 'rank':
				# Reduce the population according to their fitness
				sorted_pop = sorted(self.P.pop, key=lambda x: x.fitness, reverse=True)
				FinalPop.pop = sorted_pop[0:self.pop_size]

			# [TOURNAMENT method]
			elif self.reduction_method == 'tournament':
				# Reduce the population according to a tournament selection
				while len(FinalPop.pop) < self.pop_size:
					selected_individual = self.pop_tournament_reduction() 
					FinalPop.pop.append(selected_individual)
			else:
				print('[ERROR] The reduction method is wrong, choose between '\
				      'tournament or rank, in the config.ini file')

			y = time.time()
			z = round((y - x), 2)
			# print(f'------- REDUCTION prend {z} sec.')



			# Save the best model of the current run
			self.save_best_model(archive_instance, elite, self.verbose)

			# Erase all individual from P.pop and replace them by individuals from FinalPop.pop
			self.P.pop.clear() 
			self.P.pop = copy.deepcopy(FinalPop.pop)





			for indv in self.P.pop:
				# Metrics
				avgFitness.append(indv.fitness) # fitness of each individual
				avgRuleCount.append(len(indv.rules)) # number of rules of the individual


			# # Metrics - std. dev of fitness/rule count for each run
			dev_fitness_count.append(np.std(avgFitness))
			dev_rule_count.append(np.std(avgRuleCount))
			
			# Metrics - Best fitness/rule count
			bestFitness = round(bestIndividual.fitness, 3)
			bestRuleCount = len(bestIndividual.rules)

			# Metrics - Average of fitness/rule count
			avgFitness = round(mean(avgFitness), 3)
			avgRuleCount = int(round(mean(avgRuleCount), 0))

			
			# # Metrics
			# self.best_fitness_list.append(bestFitness)
			# self.best_rule_count.append(bestRuleCount)
			# self.avg_fitness_list.append(avgFitness)
			# self.avg_rule_count.append(avgRuleCount)





			# count the number of mutations
			self.f_add_rule.append(dico_cnt_mut['cnt_add_rule'])
			self.f_remove_rule.append(dico_cnt_mut['cnt_remove_rule'])
			self.f_replace_rule.append(dico_cnt_mut['cnt_replace_rule'])
			self.f_replace_subtree.append(dico_cnt_mut['cnt_replace_subtree'])
			self.f_add_aa.append(dico_cnt_mut['cnt_add_aa'])
			self.f_replace_node.append(dico_cnt_mut['cnt_replace_node'])
			self.f_remove_subtree.append(dico_cnt_mut['cnt_remove_subtree'])

			




			# Time/Run
			stop_time = time.time()
			diff_time = stop_time - start_time
			# run_time.append(diff_time)

			print_string = f"\n\nRun {RUN}: -best {bestFitness} -brc {bestRuleCount} || -avg {avgFitness} -arc {avgRuleCount} "
			print(print_string)

			print('\n     [- Best individual -]\n')
			bestIndividual.print()

			log_string = f"{RUN};{bestFitness};{bestRuleCount};{avgFitness};{avgRuleCount};{len(self.P.pop)};{diff_time};{eval_time};{crossover_time};{mut}"
			
			# Log the evolution information
			archive_instance.saveEvo(log_string)

			# ajoute la meilleure fitness pour voir si on doit arreter l'algo ou non
			self.stop_fitness.append(bestFitness)


			

			


		###########################################
		###                                     ###
		###              GRAPHICS               ###
		###                                     ###
		###########################################

		# lance le model sur le jeu de donnees pour voir s'il y a une correlation
		evaluate_instance.regex_matching()


		# cree le dataframe avec les Scores en fonction des cest value
		df_eval = evaluate_instance.df
		df_mock = evaluate_instance.df_mock


		# df des metrics
		df_ref = pd.read_csv('output/evo.csv', sep=';')


		dico_mutation = {	'add_rule':self.f_add_rule,
							'remove_rule':self.f_remove_rule,
							'replace_rule':self.f_replace_rule,
							'replace_subtree':self.f_replace_subtree,
							'add_aa':self.f_add_aa,
							'replace_node':self.f_replace_node,
							'remove_subtree':self.f_remove_subtree }

		df_mut = pd.DataFrame.from_dict(dico_mutation)

		df_mut_T = df_mut.T # invert the df

		# nb_line, nb_col = 3, 2

		plt.rcParams['axes.spines.right'] = False
		plt.rcParams['axes.spines.top'] = False
		plt.rcParams['font.size'] = 10


		# self.plot_fitness(df_ref, dev_fitness_count)
		# self.plot_rule_count(df_ref, dev_rule_count)
		# self.plot_mutations(df_mut)


		fig, ax = plt.subplots(3, 2)

		# Best and Average Fitness
		ax[0,0].errorbar(df_ref['Run'], df_ref['average_fitness'], yerr=dev_fitness_count, label='Avg',linewidth=1, ecolor='black', elinewidth=0.3, color='#f39c12')
		sns.lineplot(ax=ax[0,0], data=df_ref, x='Run', y='best_fitness', color='#27ae60', linewidth=1, label='Best')
		ax[0,0].set_ylabel('Fitness')
		ax[0,0].set_xlabel('Generation')

		# Rule count
		ax[0,1].errorbar(df_ref['Run'], df_ref['average_rule_count'], yerr=dev_rule_count, label='Avg',linewidth=1, ecolor='black', elinewidth=0.3, color='#f39c12')
		sns.lineplot(ax=ax[0,1], data=df_ref, x='Run', y='best_rule_count', color='#27ae60', linewidth=1, label='Best')
		ax[0,1].set_ylabel('Number of rules')
		ax[0,1].set_xlabel('Generation')


		# mutations violin plot
		sns.violinplot(ax=ax[1,0], data=df_mut, palette="Set3", bw=.2, cut=1, linewidth=1, yticklabels=['+R', '-R', '*R', '*ST', '+AA', '*N', '-ST'])
		ax[1,0].set_ylabel('Number of mutations')

		# Mutation heatmap
		a = sns.heatmap(ax=ax[1,1], data=df_mut_T, yticklabels=['+R', '-R', '*R', '*ST', '+AA', '*N', '-ST'], cmap='RdYlGn')
		ax[1,0].set_xlabel('Generation')

		# Time
		sns.lineplot(ax=ax[2,0], data=df_ref, x='Run', y='Time', color='#2980b9', linewidth=1, label='Global')
		sns.lineplot(ax=ax[2,0], data=df_ref, x='Run', y='Eval_time', color='orange', linewidth=1, label='Eval.')
		sns.lineplot(ax=ax[2,0], data=df_ref, x='Run', y='Xover_time', color='green', linewidth=1, label='Xover.')
		sns.lineplot(ax=ax[2,0], data=df_ref, x='Run', y='Mut_time', color='#B233FF', linewidth=1, label='Mutation')
		

		ax[2,0].set_ylabel('Seconds')
		ax[2,0].set_xlabel('Generation')

		# Evaluation
		sns.regplot(ax=ax[2,1], x='CEST_value',y='Predicted_score', data=df_eval,  line_kws={'color':'lightblue', 'lw':0.5 }, scatter_kws={'color':'lightblue', "s": 1, 'alpha':0.5})
		sns.regplot(ax=ax[2,1], x='CEST_value',y='Predicted_score', data=df_mock,  line_kws={'color':'orange',  'lw':0.5 }, scatter_kws={'color':'orange', "s": 1, 'alpha':0.5})
		# sns.jointplot(x='CEST_value',y='Predicted_score', data=df, color='red', line_kws={'color':'blue', 'lw':0.5 }, scatter_kws={'color':'green', "s": 1, 'alpha':0.5}, kind="reg")
		# ax[2,1].scatter(evaluate.final_dico.keys(),evaluate.final_dico.values(), marker='x')
		ax[2,1].set_ylabel('Score')
		ax[2,1].set_xlabel('CEST value')
		# compute r2
		r, p = sp.stats.pearsonr(list(evaluate_instance.final_dico.keys()), list(evaluate_instance.final_dico.values()))
		rm, pm = sp.stats.pearsonr(list(evaluate_instance.final_dico_mock.keys()), list(evaluate_instance.final_dico_mock.values()))
		ax[2,1].text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p))
		ax[2,1].text(4, 5, '(Mock) r={:.2f}, p={:.2g}'.format(rm, pm))



		plt.tight_layout()
		plt.savefig("Results", dpi=300)
		plt.subplots_adjust(hspace = 0.4)

		plt.show()


		# for i in range(2):
		# 	for j in range(2):
		# 		ax[i,j].legend()
		# # for i in ax:
		# # 	i.legend()
		# plt.show()

		print('The Best individual is:')
		bestIndividual.print()	


	def first_evaluation(self, fitness_object, bestIndividual, verbose):
		'''
		Evaluate each individual according to the fitness function. The fitness
		of an individiual is modified internally after the function evaluate_individual().

		fitness_object: Instance of Fitness Class
		bestIndividual: The current best Individual

		return a new best individual
		'''

		if verbose == 'True':
			print('[INFO] First evaluation step')

		best_fitness_value = bestIndividual.fitness


		for individual in self.P.pop:

			# Calculate the fitness value of each individual
			# individual.fitness = fitness_object.evaluate_individual(individual, lock)
			_ = fitness_object.evaluate_individual(individual)

			individual_fitness_value = individual.fitness

			# Define the best individual
			if self.config["fitness_alg"] == "RMSE2" and individual_fitness_value < best_fitness_value:
				bestIndividual = copy.deepcopy(individual)
			elif self.config["fitness_alg"] == "motif":
				if individual_fitness_value > best_fitness_value:
					bestIndividual = copy.deepcopy(individual)
					best_fitness_value = bestIndividual.fitness

		return bestIndividual


	def plot_fitness(self, df, dev_fitness, name='Fitness_plot'):
		plt.errorbar(df['Run'], df['average_fitness'], yerr=dev_fitness, label='Avg',linewidth=2, ecolor='black', elinewidth=0.5, color='#f39c12')
		sns.lineplot(data=df, x='Run', y='best_fitness', color='#27ae60', linewidth=2, label='Best')

		plt.xlabel('Generation')
		plt.ylabel('Fitness')

		plt.legend(loc='lower right')
		plt.tight_layout()

		plt.savefig(name, dpi=300)
		plt.close()

	def plot_rule_count(self, df, dev_rule_count, name='Rule_count_plot'):
		plt.errorbar(df['Run'], df['average_rule_count'], yerr=dev_rule_count, label='Avg',linewidth=1, ecolor='black', elinewidth=0.3, color='#f39c12')
		sns.lineplot(data=df, x='Run', y='best_rule_count', color='#27ae60', linewidth=1, label='Best')
		plt.ylabel('Number of rules')
		plt.xlabel('Generation')
		plt.legend(loc='lower right')
		plt.tight_layout()

		plt.savefig(name, dpi=300)
		plt.close()

	def plot_mutations(self, df, name='Mutations_plot'):
		# sns.heatmap(data=df, xticklabels=['+R', '-R', '*R', '*ST', '+AA', '*N', '-ST'], cmap='RdYlGn')
		# sns.lineplot(data=df, linewidth=0.5, dashes=False, palette=['#8e44ad', '#2980b9', '#16a085', '#2ecc71', '#f1c40f', '#e67e22', '#c0392b'])
		sns.violinplot(data=df, palette="Set3", bw=.2, cut=1, linewidth=1, rotation=40)

		plt.xlabel('Generation')
		plt.ylabel('Mutation')


		plt.tight_layout()

		plt.savefig(name, dpi=300)
		plt.close()





	def pop_tournament_reduction(self):
		'''
		Tournament for reduction step. Draw randomly X individuals (defaut=3)
		then extract the best (highest fitness)

		return the best selected individual
		'''

		tournament = []

		while len(tournament) < self.reduction_tournament_size: # par defaut 3
			# Random draw of an individual and Avoid to draw the same individual twice in a tournament
			tournament = R.sample(self.P.pop, self.reduction_tournament_size)
			# Avoid individual without rules
			tournament = [indv for indv in tournament]

		# Sort individuals according to their fitness value
		selected_indi = max(tournament, key=lambda x: x.fitness)

		return selected_indi


	# a crossover between 2 rules, not between 2 individuals # NS
	def internal_crossover(self):

		pass


	


	def crossover_single(self, parentA, parentB):
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
		nb_max_rules = lenA
		rules = []

		if lenA < lenB:
			nb_max_rules = lenB

		list_pattern = []

		# we keep track of the rules we want to add to the offspring
		for j in range(nb_max_rules):

			if j < lenA:
				ruleA = parentA.rules[j]

				# if ruleA.status == 1:
					# rules.append(ruleA)
				# elif R.random() < self.crossRate:  # we give unused rules some chance to get selected
				if R.random() < 0.5:
					if ruleA.pattern not in list_pattern:	
						rules.append(ruleA)
						list_pattern.append(ruleA.pattern)

			if j < lenB:
				ruleB = parentB.rules[j]

				# if ruleB.status == 1:
					# rules.append(ruleB)
				# elif R.random() < self.crossRate:  # we give unused rules some chance to get selected
				if R.random() < 0.5:	
					if ruleB.pattern not in list_pattern:	
						rules.append(ruleB)
						list_pattern.append(ruleB.pattern)

		return rules



	# [ATTENTION]
	# Supprime les elements au milieu, qu'importe leur taille 
	# prendre en compte les status 0 . PE parce qu'on a pas de status 0
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
		on a subset of size 'self.tournamentSize' (default = 5)

		return: two best parents A and B
		'''

		tournament = R.sample(self.P.pop, self.tournamentSize)


		# while len(tournament) < self.tournamentSize:
		# 	# Random draw of an individual and Avoid to draw the same individual twice in a tournament
		# 	tournament = R.sample(self.P.pop, self.tournamentSize)
		# 	# Avoid individual without rules
		# 	tournament = [indv for indv in tournament if len(indv.rules) != 0]

		sorted_indices = list(reversed(sorted(tournament, key=lambda x: x.fitness)))
		parentA, parentB = sorted_indices[0], sorted_indices[1]

		return parentA, parentB

	# def sort_tournament(self, tournament):
	# 	'''
	# 	Function to sort each individual in the tournament according to their fitness value
	# 	tournament: Array with *self.tournamentSize* individuals

	# 	return: two best parents A and B

	# 	'''

	# 	fitness_values = [individual.fitness for individual in tournament]
	# 	sorted_indices = sorted(range(len(fitness_values)), key=lambda i: fitness_values[i])
		
	# 	parentA = copy.deepcopy(tournament[sorted_indices[-1]])
	# 	parentB = copy.deepcopy(tournament[sorted_indices[-2]])

	# 	return parentA, parentB

	def crossover_dual(self, parentA, parentB):
		'''	
		The crossover operator cuts each individual into two sub-parts, 
		each containing a list of regex. Then, a part of individual 1 is exchanged 
		with a part of individual 2, creating two offspring that are recombinations
		of the two parents.

		Parent A/B: an individual

		return two list of rules for the 2 futur offspring	
		'''
		
		if len(parentA.rules) == 0 or len(parentB.rules) == 0:
			print('[ERROR] Optimizer.py - Crossover, a parent is empty')
			print("PARENT A")
			parentA.print()
			print("PARENT B")
			parentB.print()
			exit()
		
		else:
			# define the limits where to cut the table of rules
			cutA = R.randint(0, len(parentA.rules)-1)
			cutB = R.randint(0, len(parentB.rules)-1)

			rules_off1 = copy.deepcopy(parentA.rules[:cutA]) + copy.deepcopy(parentB.rules[cutB:])
			rules_off2 = copy.deepcopy(parentB.rules[:cutB]) + copy.deepcopy(parentA.rules[cutA:])

		return rules_off1, rules_off2

	def run_crossover(self, individual_object, newPop):
		'''
		Run the crossover step

		individual_object: instance of the class Individual
		newPop: New population where to add the offspring
		'''

		# Create 2 new empty offspring
		offspring1 = individual_object.Individual(self.config) 
		offspring2 = individual_object.Individual(self.config) 

		d = time.time()
		# Selection of 2 parents with Tournament method (default pressur selection=5)
		parentA, x = self.tournament()
		parentB, x = self.tournament()
		e = time.time()
		ee=e-d

		f = time.time()
		# Crossover operator: obtain 2 new list of rules
		r1, r2 = self.crossover_dual(parentA, parentB)
		g = time.time()
		gg=g-f


		# Define the rules of the new offspring
		offspring1.rules = r1
		offspring2.rules = r2

		a = time.time()
		# remove duplicate/close rules
		# offspring1.remove_double_rules()
		# offspring2.remove_double_rules()
		b = time.time()
		cc = b-a


		# Sort rules according to their complexity
		# offspring1.complexitySort() 
		# offspring2.complexitySort() 
		
		# Resize the offspring so it doesn't exceed the maximum allowed count
		if len(offspring1.rules) > self.ruleCount:
			offspring1.rules = offspring1.rules[:self.ruleCount]
		if len(offspring2.rules) > self.ruleCount:
			offspring2.rules = offspring2.rules[:self.ruleCount]

		h= time.time()
		# Calculate the fitness of the offspring
		if self.config["fitness_alg"] == "motif":
			offspring1.compute_fitness()
			offspring2.compute_fitness()
		j= time.time()
		jj=j-h

		# Insertion of the new offspring in the population
		newPop.pop.append(offspring1)
		newPop.pop.append(offspring2)

		return ee,gg,cc,jj


	def run_mutations(self, indv, dico_cnt_mut):
		# add a new rule/regex
		if len(indv.rules) < self.ruleCount:
			if R.random() <= self.mAR:
				self.MUT_add_rule(indv, dico_cnt_mut)

		# remove a rule/regex
		if len(indv.rules) > 1:
			if R.random() <= self.mRR:
				self.MUT_remove_rule(indv, dico_cnt_mut)

		# replace a rule/regex
		if len(indv.rules) >= 1:
			if R.random() <= self.mReR:
				self.MUT_replace_rule(indv, dico_cnt_mut)

		# le pb si on active cette mutation (qui marche), on obtient des individus
		# qui ont X fois a peu de chose pres les meme regex. ils convergent tous 
		# vers la quasi meme solution.
		# pour eviter ca, on calcule le & d identite entre les regex, et sir le %
		# on supprime celle avec le plus petit score
		# Replace a subtree (modify a part/branch of the tree)
		if R.random() <= self.mReS:
			self.MUT_replace_subtree(R.choice(indv.rules), dico_cnt_mut)
		
		# mutation only with RMSE methods
		if self.config["fitness_alg"] == "RMSE2":
			# change weight of the rule
			if R.random() <= self.mCW:
				self.MUT_change_weight(R.choice(indv.rules), dico_cnt_mut)

		# Add X new AA in a leaf node
		if R.random() <= self.mAA:
			self.MUT_add_alphabet(R.choice(indv.rules), dico_cnt_mut)

		# Replace the value of one node (invert or replace)
		if R.random() <= self.mRN:
			self.MUT_replace_node(R.choice(indv.rules), dico_cnt_mut)


		# A Fverifer

		# remove from pattern
		if R.random() <= self.mRFP:
			# pick a random node (avoid None node and root)
			self.MUT_remove_from_pattern(R.choice(indv.rules), dico_cnt_mut)




	# Add a random rule mutation
	def MUT_add_rule(self, individual, dico):
		'''
		Mutation operator. Add a new regex to the list of rules of an individual

		individual: the target individual selected with a probability of *self.mAR*
		'''

		if len(individual.rules) >= self.ruleCount:
			pass
		else:
			# initialize a random weight
			weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

			# Create a new RE pattern
			if self.init_method == 'full':
				pattern_re, tree = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
			elif self.init_method == 'grow':
				pattern_re, tree = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)
			elif self.init_method == 'half':
				pattern_re, tree = regex.indi_half(self.depth_tree, self.min_braces, self.max_braces)


			rule = Rule.Rule(pattern_re, weight, 0, tree)
			individual.rules.append(rule)

			dico['cnt_add_rule'] += 1
	
	# Remove rule mutation
	def MUT_remove_rule(self, individual, dico):
		'''
		Mutation operator. Remove a rule to the list of rules of an individual

		individual: the target individual selected with a probability of *self.mRR*
		'''

		if len(individual.rules) == 0 or len(individual.rules) == 1:
			print("[ERROR] Optimizer.py (remove rule mutation) Not enough rules\n>>>")
			exit()
		else:
			tempRand = R.randint(0, len(individual.rules) - 1)
			del (individual.rules[tempRand])

			dico['cnt_remove_rule'] += 1

	# Replace a rule by another one
	def MUT_replace_rule(self, individual, dico):
		# initialize a random weight
		weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

		tempRand = R.randint(0, len(individual.rules) - 1)
		del (individual.rules[tempRand])

		# Create a new RE pattern
		if self.init_method == 'full':
			pattern_re, tree = regex.indi_full(self.depth_tree, self.min_braces, self.max_braces)
		elif self.init_method == 'grow':
			pattern_re, tree = regex.indi_grow(self.depth_tree, self.min_braces, self.max_braces)
		elif self.init_method == 'half':
			pattern_re, tree = regex.indi_half(self.depth_tree, self.min_braces, self.max_braces)

		rule = Rule.Rule(pattern_re, weight, 0, tree)
		individual.rules.append(rule)

		dico['cnt_replace_rule'] += 1

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
	def MUT_replace_subtree(self, rule, dico):
		# Draw a random node
		random_node = self.pick_a_node(rule.tree_shape)

		# Define the layer of the random node
		for num_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = num_layer

		# Leaf case, node in the last layer
		if layer == self.last_layer:
			# Change the leaf and modify tree shape and pattern of the rule
			# if the leaf is a list
			if type(rule.tree_shape[random_node]) == list:

				if ('^' in rule.tree_shape[random_node]): # [^] case
					rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint( int(len(regex.ALPHABET)/2), len(regex.ALPHABET)-1))
					rule.tree_shape[random_node].insert(0,'^')

				else: # [] case
					rule.tree_shape[random_node] = R.sample(regex.ALPHABET, R.randint(1, int(len(regex.ALPHABET)/2)+1 ))

			# if the leaf is a single element
			else:
				rule.tree_shape[random_node] = R.choice(regex.LAST)

		# Operator node case
		else:
			# Generate a new sub-tree according to a initialization method
			tree_method = np.random.choice(['full', 'grow'], p=[0.5,0.5])

			if tree_method == 'full':
				new_regex, new_tree = regex.indi_full( (self.depth_tree - layer)+1, self.min_braces, self.max_braces) 
			elif tree_method == 'grow':
				new_regex, new_tree = regex.indi_grow( (self.depth_tree - layer)+1, self.min_braces, self.max_braces) 

			rule.tree_shape[random_node] = new_tree[0]
	
			# Create list of children of the random node
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

		dico['cnt_replace_subtree'] += 1


	# Add some new AA in leaf node /!\ dont respect tree shape
	def MUT_add_alphabet(self, rule, dico):
		'''
		Add 1 to 4 new AA in a leaf node. Ex:
		Initial: (H{4}|(P+))|(O+)
		Tree shape: ['|', '|', '+', '{4}', '+', 'O', None, 'HOHAA', None, 'P', None, None, None, None, None]
		New:(HOHAA{4}|(P+))|(O+)
			  ****

		Rule: The rule/regex impacted by the mutation
		random_node: the node impacted by the mutation, randomly selected
		'''

		# Draw a random node
		random_node = self.pick_a_node(rule.tree_shape)

		# Define the layer of the random node
		for key_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = key_layer

		if layer == self.last_layer: # Modify only leaf node (no list)
			if type(rule.tree_shape[random_node]) != list:
				nbr_new_aa = R.randint(1, 4) # Number of AA to add

				for i in range(nbr_new_aa):
					aa = R.choice(regex.LAST)
					rule.tree_shape[random_node] += aa # Add new AA

				# Aply modification on the pattern of the rule
				rule.pattern = regex.tree2regex(rule.tree_shape)
		dico['cnt_add_aa'] += 1

	# Change/invert one node
	def MUT_replace_node(self, rule, dico):


		# pick a random node (avoid None node)
		random_node = self.pick_a_node(rule.tree_shape)
		
		if rule.tree_shape[random_node] == 'cat':
			rule.tree_shape[random_node] = '|'

		elif rule.tree_shape[random_node] == '|':
			rule.tree_shape[random_node] = 'cat'
		
		elif '{' in rule.tree_shape[random_node]:
			rule.tree_shape[random_node] = '{' + str(R.randint(self.min_braces, self.max_braces)) + '}'

		elif rule.tree_shape[random_node] == '[]':
			try:
				rule.tree_shape[random_node] = '[^]'
				rule.tree_shape[(random_node*2) + 1].insert(0,'^')
				
			except AttributeError:
				print('[WARNING] optimizer.py - mutation replace_node doesn\'t find the child of [^]')
				print(rule.tree_shape,rule.tree_shape[random_node], random_node)
				print(rule.pattern)
				pass

		elif rule.tree_shape[random_node] == '[^]':
			rule.tree_shape[random_node] = '[]'
			if ('^' in rule.tree_shape[(random_node*2) + 1]):
				rule.tree_shape[(random_node*2) + 1].remove('^')
		else:
			pass

		rule.pattern = regex.tree2regex(rule.tree_shape)
		dico['cnt_replace_node'] += 1

	# Alter patterns mutation (remove letter)    ------> DONT WORK
	def MUT_remove_from_pattern(self, rule, dico):
		copysave = copy.deepcopy(rule.tree_shape)



		# print("\n\n>>> AVANT",rule.tree_shape)
		# print(">>> AVANT",rule.pattern)
		random_node = self.pick_a_node(rule.tree_shape)
		# print('------------------------------------------------------')
		# print('Noeud supprime', random_node, rule.tree_shape[random_node])
		# self.deleteTree(random_node, rule)



		# Define the layer of the random node
		for key_layer, list_nodes in self.dict_layer.items():
			if random_node in list_nodes:
				layer = key_layer

		brother = self.I_am_ur_brother(random_node)
		child = self.create_subtree_list(random_node)
		parent = self.I_am_ur_parent(random_node)


		# print('frere', brother)
		# print('enfant', child)
		# print('parent', parent, rule.tree_shape[parent])
		# rule.tree_shape[random_node] = None

		try:
			if '[' in rule.tree_shape[random_node]:
	
				# rule.tree_shape[random_node] = None
				rule.tree_shape[parent] = 'cat'
	
			for i in child:
				rule.tree_shape[i] = None
	
	
	
			if rule.tree_shape[parent] == '|':
				rule.tree_shape[parent] = 'cat'
			if rule.tree_shape[parent] == '+':
				rule.tree_shape[parent] = None
	
			if '{' in rule.tree_shape[parent]:
				rule.tree_shape[parent] = None
				x = self.I_am_ur_parent(parent)
				if rule.tree_shape[x] == '|':
					rule.tree_shape[x] = 'cat'
	
			if rule.tree_shape[parent] == '[^]':
				rule.tree_shape[parent] = None
			if rule.tree_shape[parent] == '[]':
				rule.tree_shape[parent] = None
		except:
			rule.tree_shape = copysave
			rule.pattern = regex.tree2regex(rule.tree_shape)



		rule.pattern = regex.tree2regex(rule.tree_shape)

		if '(|' in rule.pattern or '|)' in rule.pattern:
			rule.tree_shape = copysave
			rule.pattern = regex.tree2regex(rule.tree_shape)
		else:
			dico['cnt_remove_subtree'] += 1












		# # Define last layer, and leaf nodes
		# # right_subtree =  self.create_subtree_list(1)
		# # left_subtree = self.create_subtree_list(2)
		# # right_leaf = 0
		# # left_leaf = 0

		# # for i in right_subtree:
		# # 	if rule.tree_shape[i] == None:
		# # 		pass
		# # 	else:
		# # 		right_leaf = i

		# # for i in left_subtree:
		# # 	if rule.tree_shape[i] == None:
		# # 		pass
		# # 	else:
		# # 		left_leaf = i

		# if layer == self.last_layer: # Never delete a leaf node
		# 	pass
		# # elif random_node == right_leaf: # Never delete a leaf node
		# # 	pass
		# # elif random_node == left_leaf: # Never delete a leaf node
		# # 	pass
		# else:

		# 	# Define list of node in child subtree
		# 	child = self.create_subtree_list(random_node)
		# 	# print('j influence', child)

		# 	# Replace all child nodes by None value
		# 	for i, node_index in enumerate(child):
		# 		rule.tree_shape[node_index] = None

		# 	# Define parent and brother of the random node
		# 	parent = self.I_am_ur_parent(random_node)
		# 	brother = self.I_am_ur_brother(random_node)

		# 	# print('parent',parent, rule.tree_shape[parent])
		# 	# print('brother',brother, rule.tree_shape[brother])

		# 	# Define list of node in brother and parent subtree
		# 	brother_tree = self.create_subtree_list(brother)
		# 	parent_tree = self.create_subtree_list(parent)
		# 	# print('parent_tree',parent_tree )
		# 	# print('brother_tree',brother_tree)

		# 	# Replace parent nodes by brother nodes
		# 	for i, j in enumerate(parent_tree):
		# 		try:
		# 			rule.tree_shape[parent_tree[i]] = rule.tree_shape[brother_tree[i]]
		# 		except IndexError:
		# 			rule.tree_shape[parent_tree[i]] = None

		
		# # Checkpoint
		# if rule.tree_shape[0] != 'cat' and rule.tree_shape[0] != '|':
		# 	rule.tree_shape[1] = rule.tree_shape[0]
		# 	rule.tree_shape[0] = 'cat'

		# try:
		# 	rule.pattern = regex.tree2regex(rule.tree_shape)
		# except:
		# 	# Error of compilation, tree shape is erroneous. remove changing
		# 	rule.tree_shape = copysave
		# 	rule.pattern = regex.tree2regex(rule.tree_shape)


		

	

		


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




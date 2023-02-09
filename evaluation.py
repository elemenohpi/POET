#!/usr/bin/env python3
# @autor: N. Scalzitti
# @date: 11/01/2022

import re
import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


class Eval:
	def __init__(self, config):
		self.config = config
		self.model = config['output_model']
		self.dico_regex = {}
		self.final_dico = {}
		self.final_dico_mock = {}
		self.training_set = 'data/Training/' + config['learn_data'] + '.csv'
		self.mock_data = 'data/mock.csv'

		self.df = ''
		self.df_mock = ''



	def sort_dict_by_value(self, dico):
	    return dict(sorted(dico.items(), key=lambda item: item[0]))

	def model_to_dict(self):
		with open(self.model, 'r') as csv_file:
			file = csv.reader(csv_file, delimiter=',')
			next(file) # saute le header
			for line in file:
				regex = line[2]
				score = float(line[3])
				self.dico_regex[regex] = score

		return self.dico_regex

	def regex_matching(self):
		self.model_to_dict()

		data = {'CEST_value': [], 'Predicted_score': [],}
		data_mock = {'CEST_value': [], 'Predicted_score': [],}

		with open(self.training_set, 'r') as csv_file:
			file = csv.reader(csv_file, delimiter=',')
			next(file) # saute le header
			for line in file:
				peptide = line[0]
				cest = line[1]
				final_score = 0

				for regex, score in self.dico_regex.items():
					object_find = re.finditer(regex, peptide)
					
					for match in object_find:
						motif = match.group()
						if len(motif) >= 3:
							final_score += score # ne prend pas en compte la taille des motifs trouve

				data['Predicted_score'].append(float(final_score))
				data['CEST_value'].append(float(cest))

				self.final_dico[float(cest)] = final_score

		with open(self.mock_data, 'r') as csv_file:
			file = csv.reader(csv_file, delimiter=',')
			next(file) # saute le header
			for line in file:
				peptide = line[0]
				cest = line[1]
				final_score = 0

				for regex, score in self.dico_regex.items():
					object_find = re.finditer(regex, peptide)
					
					for match in object_find:
						motif = match.group()
						if len(motif) >= 3:
							final_score += score

				data_mock['Predicted_score'].append(float(final_score))
				data_mock['CEST_value'].append(float(cest))
				self.final_dico_mock[float(cest)] = final_score

		self.df = pd.DataFrame.from_dict(data)
		self.df_mock = pd.DataFrame.from_dict(data_mock)

	def plot_evaluation(self):
		sns.regplot(x=list(self.final_dico.keys()),y=list(self.final_dico.values()), color='red')
		plt.scatter(self.final_dico.keys(),self.final_dico.values(), alpha=0.8)

		sns.regplot(x=list(self.final_dico_mock.keys()),y=list(self.final_dico_mock.values()), color='orange')
		plt.scatter(self.final_dico_mock.keys(),self.final_dico_mock.values(), color='lightblue', alpha=0.8)
		
		plt.ylabel('Score')
		plt.xlabel('CEST value')

		plt.show()


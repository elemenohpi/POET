#!/usr/bin/env python3
# @autor: N. Scalzitti
# @date: 01/20/2022

import matplotlib.pyplot as plt
from statistics import mean


MOTIF_BANK = {}

def open_file(filename):
	all_seq = []
	all_cest = []
	all_class = []
	with open(f"data/Training/{filename}.csv", 'r') as file_R1:
		for line in file_R1:
			line = line.strip().split(',')
			sequence = line[0]
			cest = line[1]
			classe = line[2]
			all_seq.append(sequence)
			all_cest.append(cest)
			all_class.append(classe)

	return all_seq, all_cest, all_class

# Create a motif bank for the training step

# - standard_motif_bank = extract motif with best CEST value
# - average_motif_bank =  compute the average CEST value
# - occ_motif_bank = extract the number of occurence in each classe

def standard_motif_bank(motif_size, alphabet):
	'''
	Create a Dict() containing motifs and their class for the training step
	motif_size: size of extracted motifs (2-6)
	alphabet: Name of the alphabet (ex: Classical, Zappo...)

	return a Dict() with key=motif, value = [occurence, best CEST, classe]

	'''

	list_seq, list_cest, list_class = open_file(alphabet)

	for i, seq in enumerate(list_seq):
		if i != 0: # header
			for j in range(len(seq)):
				try:
					motif = seq[j:j+motif_size]
					cest = float(list_cest[i])
					classe = int(list_class[i])

					if len(motif) == motif_size:
						if '-' in motif:
							pass
						else:
							if motif not in MOTIF_BANK.keys(): # unknow motif
								MOTIF_BANK[motif] = [1, cest, classe]
							else: # known motif
								MOTIF_BANK[motif][0] += 1

								# change the CEST value if it is better
								if cest > MOTIF_BANK[motif][1]:
									MOTIF_BANK[motif][1] = cest

									if cest > 12.5:
										MOTIF_BANK[motif][1] = 1
				except:
					pass

	return MOTIF_BANK

def average_motif_bank(motif_size, alphabet):
	'''
	Calculate the average CEST, if > 12.5 then group 1 otherwise group 0
	motif_size: size of extracted motifs (2-6)
	alphabet: Name of the alphabet (ex: Classical, Zappo...)

	return a Dict() with key=motif, value = [occurence, CEST mean, classe]

	'''

	# Extract information of the training file
	list_seq, list_cest, list_class = open_file(alphabet)
	TPR_DICO = {}

	for i, seq in enumerate(list_seq):

		if i != 0: # header
			for j in range(len(seq)): # create sliding windows

				try:
					motif = seq[j:j+motif_size]
					cest = float(list_cest[i])
					classe = int(list_class[i])

					if len(motif) == motif_size:
						if '-' in motif:
							pass
						else:
							if motif not in TPR_DICO.keys():  # new motif
								TPR_DICO[motif] = [1, [cest] ]
							else: # Known motif
								# Add a new value of CEST associated with the motif in the training set
								TPR_DICO[motif][1].append(cest)
								TPR_DICO[motif][0] += 1
				except:
					pass

	# Compute the mean of CEST value 

	for k, v in TPR_DICO.items():
		MOTIF_BANK[k] = []
		MOTIF_BANK[k].append(v[0]) # Number of motifs in all dataset
		MOTIF_BANK[k].append(round(mean(v[1]),2)) # compute the CEST mean
		
		# Change the classe
		if mean(v[1]) >= 12.5:
			MOTIF_BANK[k].append(1)
		else:
			MOTIF_BANK[k].append(0)

	return MOTIF_BANK

def occ_motif_bank(motif_size, alphabet):
	'''
	Define the classe of the motif. Counting the number of occurence of classe 1 and classe 0
	if nbr classe 1 > nbr classe 0, then motif is classe 1, otherwise it is classe 0.
	if it is the same number, it compute the difference between CEST value and 12.5.
	Ex, the CEST value of motif class 1 is 28 and the same motif have another CEST value of 3 (classe 0)
	we compute :|12.5-3| = 9.5 and |12.5 - 28| = 15.5, then 15.5 > 9.5 so it classe 1

	motif_size: size of extracted motifs (2-6)
	alphabet: Name of the alphabet (ex: Classical, Zappo...)

	return a Dict() with key=motif, value = [occurence, CEST , classe]

	'''
	list_seq, list_cest, list_class = open_file(alphabet)
	TPR_DICO = {}

	for i, seq in enumerate(list_seq):

		if i != 0: # header
			for j in range(len(seq)): # create sliding windows

				try:
					motif = seq[j:j+motif_size]
					cest = float(list_cest[i])
					classe = int(list_class[i])

					if len(motif) == motif_size:
						if '-' in motif:
							pass
						else:
							if motif not in TPR_DICO.keys():  # new motif
								if classe == 0:
									TPR_DICO[motif] = [1, cest,0, 1, 0]
								else:
									TPR_DICO[motif] = [1, 12.5,cest, 0, 1]

							else:
								# add 1 to the motif counter
								if classe == 0:
									TPR_DICO[motif][3] += 1
									TPR_DICO[motif][0] += 1

									if cest < TPR_DICO[motif][1]:
										TPR_DICO[motif][1] = cest
								else:
									TPR_DICO[motif][4] += 1
									TPR_DICO[motif][0] += 1

									if cest > TPR_DICO[motif][2]:
										TPR_DICO[motif][2] = cest
				except:
					pass

	for k, v in TPR_DICO.items():

		MOTIF_BANK[k] = []
		MOTIF_BANK[k].append(v[0]) # Number of motifs in all dataset

		if TPR_DICO[k][3] > TPR_DICO[k][4]:
			MOTIF_BANK[k].append(v[1])
			MOTIF_BANK[k].append(0)
		elif TPR_DICO[k][4] > TPR_DICO[k][3]:
			MOTIF_BANK[k].append(v[2])
			MOTIF_BANK[k].append(1)

		elif TPR_DICO[k][3] == TPR_DICO[k][4]:

			classe0 = abs(12.5 - v[1])
			classe1 = abs(12.5 - v[2])

			if classe0 > classe1:
				MOTIF_BANK[k].append(v[1])
				MOTIF_BANK[k].append(0)
			else:
				MOTIF_BANK[k].append(v[2])
				MOTIF_BANK[k].append(1)

	return MOTIF_BANK



def plot_motif_occurence(motif_size, alphabet):
	list_seq, list_cest, list_class = open_file(alphabet)

	for i, seq in enumerate(list_seq):
		if i != 0: # header
			for j in range(len(seq)):
				try:
					motif = seq[j:j+motif_size]
					cest = float(list_cest[i])
					classe = int(list_class[i])

					if len(motif) == motif_size:
						if '-' in motif:
							pass
						else:
							if motif not in MOTIF_BANK.keys():  # new motif
								if classe == 0:
									MOTIF_BANK[motif] = [1, cest, 1, 0]
								else:
									MOTIF_BANK[motif] = [1, cest, 0, 1]

							else:
								# add 1 to the motif counter
								if classe == 0:
									MOTIF_BANK[motif][2] += 1
								else:
									MOTIF_BANK[motif][3] += 1

								if cest > MOTIF_BANK[motif][1]: # better CEST
									MOTIF_BANK[motif][1] = cest

				except:
					pass


	# Calculer la moyenne des CEST, si > 12.5 alors groupe 1 sinon groupe 0
	# si nbr groupe 1 sup nbr groupe 0, alors groupe 1 sinon groupe 0, si == alors calculer le plus
	# grand ecarts avec 12.5, par ex, le motif du groupe 1 a un CEST de 28 et le meme motif mais 
	# dans le groupe 0 a un CEST de 3 --> |12.5-3| = 9.5 et |12.5 - 28| = 15.5, 15.5 est > donc
	# c'est un groupe 1

	l_motif = []
	l_pos = []
	l_neg = []


	for motif, value in MOTIF_BANK.items():
		l_motif.append(motif)
		l_pos.append(value[3])
		l_neg.append(value[2])


	p1 = plt.bar(l_motif, l_pos, 0.35, label='Classe 1', color='#3ea100')
	p2 = plt.bar(l_motif, l_neg, 0.35, bottom = l_pos, label='Classe 0', color='#ab0300')

	plt.ylabel('Number')
	plt.xticks(rotation=45, size=7)

	# plt.xticks(ind, ('T1', 'T2', 'T3', 'T4', 'T5'))
	# plt.yticks(np.arange(0, 81, 10))
	plt.legend()
 
	plt.show()




	return MOTIF_BANK


d=occ_motif_bank(4, 'Zappo')
print(d)


#!/usr/bin/env python3
# @autor: N. Scalzitti
# @date: 12/16/2022

# Create files for the differents Alphabet

INPUT = 'data/Training/Classical.csv'





# Lenckowski, 2007, good for alignment
lenckowski_dico = {'A':['A', 'C', 'K', 'P'],
			  	   'D':['D', 'H', 'M', 'F', 'S'],
			  	   'G':['G', 'I', 'V'], 
			  	   'N':['N', 'T', 'W', 'Y'], 
			  	   'R':['R', 'Q', 'E', 'L'],
			  	   '-':['-']}

# the alphabet based on the scoring schema from Cannata, 2002
cannata_dico = {   'I':['I', 'M', 'V', 'L', 'F', 'Y'],
			  	   'P':['P', 'G', 'A', 'T', 'S'],
			  	   'H':['H', 'N', 'Q', 'E', 'D', 'R', 'K'], 
			  	   'W':['W'], 
			  	   'C':['C'],
			  	   '-':['-']}

# Liu 2002, (the alphabet based on the residue pair counts for the MJ matrix
liu1_dico = {	   'I':['I', 'M', 'V', 'L', 'F'],
			  	   'A':['A', 'C', 'W'],
			  	   'Y':['Y', 'Q', 'H', 'P', 'G', 'T', 'S', 'N'], 
			  	   'R':['R', 'K'], 
			  	   'D':['D', 'E'],
			  	   '-':['-']}

# Liu 2002, (the alphabet based on the residue pair counts for the BLOSUM matrix
liu2_dico = {      'I':['I', 'M', 'V', 'L'],
			  	   'F':['F', 'W', 'Y'],
			  	   'P':['P', 'C', 'A', 'S', 'T'], 
			  	   'N':['N', 'H', 'Q', 'E', 'D', 'R', 'K'], 
			  	   'G':['G'],
			  	   '-':['-']}


# Murphy, 2000, alphabet based on the correlation coefficient
liu2_dico = {      'I':['I', 'M', 'V', 'L'],
			  	   'F':['F', 'W', 'Y'],
			  	   'P':['P', 'C', 'A', 'S', 'T'], 
			  	   'N':['N', 'H', 'Q', 'E', 'D', 'R', 'K'], 
			  	   'G':['G'],
			  	   '-':['-']}
# ILMVC, ASGPT, FYW, EDQN, KRH

def create_dico(alphabet_name):
	dico = {}

	with open(f'../data/Conversion/{alphabet_name}.csv') as f:
		for line in f:
			line = line.strip().split(";")
			value = line[0]
			key = line[1]
		
			if key not in dico.keys():
				dico[key] = []
				dico[key].append(value)
			else:
				dico[key].append(value)

	return dico

def create_converted_file(output_file, dico):
	'''
	Create a new dataset with converted sequence, depending of the translation dictionnary
	'''
	with open(INPUT, 'r') as file_R1:
		with open(output_file, 'w') as file_W1:
			for i, line in enumerate(file_R1):
				if i !=0:

					line = line.strip().split(',')
					seq = line[0]
					score = line[1]
					classe = line[2]

					print('>>>', seq, score)

					converted_seq = ''

					for aa in seq:
						for k,v in dico.items():
							if aa in v:
								converted_seq += k

					print(converted_seq)

					file_W1.write(f'{converted_seq};{score};{classe}\n')
				else:
					file_W1.write("sequence;fitness;Classe\n")



if __name__ == '__main__':
	# Physiochemical properties from Zappo
	zappo_dico = create_dico("Zappo")


	create_converted_file("../data/Training/Zappo.csv", zappo_dico)


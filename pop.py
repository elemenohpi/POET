# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import individual as Individual

class Population:
	# constructor for random initialization
	def __init__(self):
		self.size = settings.population_size
		self.TT = settings.TT
		self.pop = []

		print ("Initializing a population with size of " +str(self.size) + "...\n")

		if self.is_numeric():
			print ("Translation Table supports the NUMERIC mode - Generating formula...\n")
			self.populate_formulas()
		else:
			print ("Translation Table supports the PATTERN mode - Generating rules...\n")
			self.populate_rules()

	# Randomly initializes the population with rules
	def populate_rules(self):
		for i in range(self.size):
			indv = Individual.Individual()
			indv.init_pattern()
			self.pop.append(indv)

	# Randomly initializes the population with formulas
	def populate_formulas(self):
		raise("This feature is not coded yet.")
		for i in range(self.size):
			indv = Individual.Individual()
			indv.init_formula()
			self.pop.append(indv)

	# Uses the input files to create the population
	def populate_preset(self, population):
		raise ValueError('Uncharted terretories... Exiting')
		pass

	# Returns True if all the codes are numeric, otherwise returns False
	def is_numeric(self):
		codes = self.TT['code']
		for i in range(codes.size):
			try:
				float(codes[i])
			except ValueError:
				return False
		return True
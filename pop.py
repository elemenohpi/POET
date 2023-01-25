# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import individual as Individual


class Population:
	# constructor for random initialization
	def __init__(self, config):
		self.config = config
		self.size = int(config["population_size"])
		self.pop = []
		self.init_method = str(config["init_pop_met"])

		print(f"Initializing a population with size of {self.size} and with method {self.init_method.upper()}\n")
		print("Translation Table supports the PATTERN mode - Generating rules...\n")
		self.populate_rules()
			

	# Randomly initializes the population with rules
	def populate_rules(self):
		
		for i in range(self.size):
			indv = Individual.Individual(self.config)
			indv.init_pattern(self.init_method)
			# indv.print()
			self.pop.append(indv)

	# Randomly initializes the population with formulas
	def populate_formulas(self):
		raise ("This feature is not coded yet.")
		for i in range(self.size):
			indv = Individual.Individual()
			indv.init_formula()
			self.pop.append(indv)

	# Uses the input files to create the population
	def populate_preset(self, population):
		raise ValueError('Uncharted terretories... Exiting')
		pass

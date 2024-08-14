# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import individual as Individual


class Population:
    # constructor for random initialization
    def __init__(self, config):
        self.config = config
        self.size = int(config["population_size"])
        self.pop = []

        print("Initializing a population with size of " + str(self.size) + "...\n")

        print("Translation Table supports the PATTERN mode - Generating rules...\n")
        self.populate_rules()

    # Randomly initializes the population with rules
    def populate_rules(self):
        def create_and_append_individual(init_method):
            indv = Individual.Individual(self.config)
            init_method(indv)
            self.pop.append(indv)

        if self.config['rules'] is None:
            for _ in range(self.size):
                create_and_append_individual(lambda indv: indv.init_pattern())
        else:
            half_size = int(self.size / 2)
            for _ in range(half_size):
                create_and_append_individual(lambda indv: indv.makeFromFile(self.config['rules']))
            for _ in range(self.size - half_size):
                create_and_append_individual(lambda indv: indv.init_pattern())

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

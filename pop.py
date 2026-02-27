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
        for i in range(self.size):
            indv = Individual.Individual(self.config)
            indv.init_pattern()
            self.pop.append(indv)

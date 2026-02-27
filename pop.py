import individual as Individual


class Population:
    # constructor for random initialization
    def __init__(self, config):
        self.config = config
        self.size = int(config["population_size"])
        self.pop = []
        self.matching_mode = config.get("matching_mode", "substring")

        print("Initializing a population with size of " + str(self.size) + "...\n")

        if self.matching_mode == "regex":
            import regex_tree

            print("Using REGEX matching mode - Generating regex rules...\n")
            # Ensure regex_tree module is initialized before creating individuals
            if not regex_tree._initialized:
                alphabet_file = config.get(
                    "regex_alphabet", "data/translation/regex_alphabet.csv"
                )
                regex_tree.init(alphabet_file)
        else:
            print("Using SUBSTRING matching mode - Generating pattern rules...\n")

        self.populate_rules()

    # Randomly initializes the population with rules
    def populate_rules(self):
        for i in range(self.size):
            indv = Individual.Individual(self.config)
            if self.matching_mode == "regex":
                indv.init_regex_pattern(self.config)
            else:
                indv.init_pattern()
            self.pop.append(indv)

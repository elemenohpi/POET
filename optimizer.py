import fitness as F
import pandas as pd
import archivist as Archivist
import copy
import random as R
import individual as I
import rule as Rule
import regex_tree


class Optimizer:
    def __init__(self, config, population):
        self.config = config
        self.P = population
        self.runs = int(config["runs"])
        self.tournamentSize = int(config["tournament_size"])
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
        self.mCWmin = 0
        self.mCWmax = 1

        # Matching mode
        self.matching_mode = config.get("matching_mode", "substring")

        if self.matching_mode == "substring":
            # Substring-mode mutation rates
            self.mATP = float(config["mut_add_to_pattern"])
            self.mRFP = float(config["mut_remove_from_pattern"])
            codes = pd.read_csv("data/translation/amino_to_amino.csv")
            self.codes = codes["code"].tolist()
        elif self.matching_mode == "regex":
            # Regex-mode config and mutation rates
            self.depth_tree = int(config.get("max_depth_tree", "4"))
            self.min_braces = int(config.get("min_braces", "1"))
            self.max_braces = int(config.get("max_braces", "3"))
            self.init_method = config.get("init_method", "half")
            self.mReR = float(config.get("mut_replace_rule", "0.1"))
            self.mReS = float(config.get("mut_replace_subtree", "0.1"))
            self.mAA = float(config.get("mut_add_aa", "0.1"))
            self.mRN = float(config.get("mut_replace_node", "0.1"))
            self.mRFP = float(config.get("mut_remove_from_pattern", "0.1"))
            # Tree helper state
            self.maxnodes = (2**self.depth_tree) - 1
            self.dict_layer: dict[int, list[int]] = {}
            self.last_layer = self._build_layer_dict()
            # Initialize the regex_tree module
            alphabet_file = config.get(
                "regex_alphabet", "data/translation/regex_alphabet.csv"
            )
            regex_tree.init(alphabet_file)
        else:
            raise ValueError(
                "Invalid matching_mode '{}': expected 'substring' or 'regex'".format(
                    self.matching_mode
                )
            )

    def optimize(self):
        fitness = F.Fitness(self.config)
        arch = Archivist.Archivist(self.config)

        num_workers = int(self.config.get("workers", 0))
        fitness.start_workers(num_workers)

        log_string = (
            "best fitness, best test, best rule count, best unused rule count, average fitness, "
            "average test, average rule count, average unused rule count "
        )
        arch.saveEvo(log_string)

        try:
            self._run_generations(fitness, arch)
        finally:
            fitness.stop_workers()

    def _run_generations(self, fitness, arch):
        for i in range(self.runs):
            # Evaluate all individuals (parallel or sequential)
            fitness.measureBatch(self.P.pop)

            avgFitness = 0.0
            avgTest = 0.0
            avgRuleCount = 0.0
            avgUsedRulesCount = 0.0
            bestIndividual = self.P.pop[0]

            for j in self.P.pop:
                avgFitness += j.fitness
                avgTest += j.test
                avgRuleCount += len(j.rules)
                avgUsedRulesCount += j.usedRulesCount
                if j.fitness < bestIndividual.fitness:
                    bestIndividual = j

            avgFitness = round(avgFitness / float(len(self.P.pop)), 3)
            avgTest = round(avgTest / float(len(self.P.pop)), 3)
            avgRuleCount = round(avgRuleCount / float(len(self.P.pop)), 0)
            avgUsedRulesCount = round(avgUsedRulesCount / float(len(self.P.pop)), 0)
            bestFitness = round(bestIndividual.fitness, 3)
            bestTest = round(bestIndividual.test, 3)
            bestRuleCount = len(bestIndividual.rules)
            bestUsedRulesCount = bestIndividual.usedRulesCount

            # Log the outcome before doing the changes to the population / generating a new population
            print_string = "{}: -b {} -bt {} -rc {} -urc {} ||| -a {} -at {} -arc {} -aurc {}".format(
                i,
                bestFitness,
                bestTest,
                bestRuleCount,
                bestUsedRulesCount,
                avgFitness,
                avgTest,
                avgRuleCount,
                avgUsedRulesCount,
            )

            log_string = "{},{},{},{},{},{},{},{},{}".format(
                i,
                bestFitness,
                bestTest,
                bestRuleCount,
                bestUsedRulesCount,
                avgFitness,
                avgTest,
                avgRuleCount,
                avgUsedRulesCount,
            )
            # Print the evolutionary log
            print(print_string, flush=True)

            # Log the evolution
            arch.saveEvo(log_string)

            # Elitism
            new_pop = [bestIndividual]

            # Save the best model
            data = []
            for rule in bestIndividual.rules:
                data.append(
                    [rule.pattern, rule.weight, rule.status, rule.match_direction]
                )
            df = pd.DataFrame(
                data, columns=["pattern", "weight", "status", "match_direction"]
            )
            arch.saveModel(df)

            # Select Parents (Tournament Selection) and Crossover
            for k in range(len(self.P.pop) - 1):
                tournament = []
                offspring = I.Individual(self.config)

                for j in range(self.tournamentSize):
                    tournament.append(self.P.pop[R.randint(0, len(self.P.pop) - 1)])

                tournament = self.sort_tournament(tournament)

                # We got two best parents
                parentA = copy.deepcopy(tournament[0])
                parentB = copy.deepcopy(tournament[1])

                # Do the crossover magic - Cluster crossover
                # Efficiency thing. Find the greater rule length
                lenA = len(parentA.rules)
                lenB = len(parentB.rules)
                maxLen = max(lenA, lenB)

                # we keep track of the rules we want to add to the offspring
                rules = []
                for j in range(maxLen):
                    if j < lenA:
                        ruleA = parentA.rules[j]
                        if ruleA.status == 1:
                            rules.append(ruleA)
                        elif (
                            R.random() < self.crossRate
                        ):  # we give unused rules some chance to get selected
                            rules.append(ruleA)
                    if j < lenB:
                        ruleB = parentB.rules[j]
                        if ruleB.status == 1:
                            rules.append(ruleB)
                        elif (
                            R.random() < self.crossRate
                        ):  # we give unused rules some chance to get selected
                            rules.append(ruleB)

                offspring.rules = rules
                offspring.bubbleSort()

                # Resize the offspring so it doesn't exceed the maximum allowed count
                while len(offspring.rules) > self.ruleCount:
                    countGreens = 0
                    for index in range(len(offspring.rules) - 1, -1, -1):
                        if countGreens >= index:
                            del offspring.rules[index]
                            break
                        else:
                            if offspring.rules[index].status == 0:
                                del offspring.rules[index]
                                break
                            else:
                                countGreens += 1

                new_pop.append(offspring)

            new_pop[0] = copy.deepcopy(new_pop[0])
            self.P.pop = new_pop

            # Mutations

            # We keep a copy of the elite
            elite = copy.deepcopy(self.P.pop[0])

            if self.matching_mode == "substring":
                self._mutate_substring(self.P.pop)
            else:
                self._mutate_regex(self.P.pop)

            zeroFitness, testData = fitness.measureTotal(self.P.pop[0])

            # Check if elite got worse
            if elite.fitness < zeroFitness:
                self.P.pop[0] = elite

    def sort_tournament(self, t):
        t.sort(key=lambda x: x.fitness)

        if self.config["diversity_selection"] == "False":
            return t[:2]
        if self.config["diversity_selection"] != "True":
            raise ValueError(
                "Unknown value for diversity_selection: '{}'".format(
                    self.config["diversity_selection"]
                )
            )
        patterns = {rule.pattern for rule in t[0].rules}

        highest_closeness = 0
        highest_fitness = 0
        for i in range(1, len(t)):
            t[i].closeness = 0
            t[i].relative_fitness = 0
            for rule in t[i].rules:
                if rule.pattern in patterns:
                    t[i].closeness += 1
            if t[i].closeness > highest_closeness:
                highest_closeness = t[i].closeness
            if t[i].fitness > highest_fitness:
                highest_fitness = t[i].fitness

        if highest_closeness == 0:
            return t[:2]

        w = float(self.config["diversity_weight"])
        for i in range(1, len(t)):
            t[i].relative_fitness = (
                t[i].fitness / highest_fitness + w * t[i].closeness / highest_closeness
            )
        t2 = t[1:]
        t2.sort(key=lambda x: x.relative_fitness)
        return [t[0], t2[0]]

    # ── Mutation dispatch ──────────────────────────────────────────────────

    def _mutate_substring(self, pop):
        """Apply substring-mode mutations to all individuals."""
        for indv in pop:
            needs_sort = False
            if R.random() <= self.mAR:
                self.mut_add_rule(indv)
                needs_sort = True
            if R.random() <= self.mRR:
                self.mut_remove_rule(indv)
            for rule in indv.rules[:]:
                if R.random() <= self.mCW:
                    self.mut_change_weight(rule)
                if R.random() <= self.mATP:
                    self.mut_add_to_pattern(rule)
                    needs_sort = True
                if R.random() <= self.mRFP:
                    self.mut_remove_from_pattern(rule)
                    needs_sort = True
                    if rule.pattern == "":
                        indv.rules.remove(rule)
            if needs_sort:
                indv.bubbleSort()

    def _mutate_regex(self, pop):
        """Apply regex-mode tree mutations to all individuals."""
        for indv in pop:
            needs_sort = False
            # Add a new regex rule
            if len(indv.rules) < self.ruleCount:
                if R.random() <= self.mAR:
                    self.regex_mut_add_rule(indv)
                    needs_sort = True
            # Remove a rule
            if len(indv.rules) > 1:
                if R.random() <= self.mRR:
                    self.mut_remove_rule(indv)
            # Replace a rule with a new regex
            if len(indv.rules) >= 1:
                if R.random() <= self.mReR:
                    self.regex_mut_replace_rule(indv)
                    needs_sort = True
            # Replace a subtree in a rule
            if R.random() <= self.mReS:
                rule = R.choice(indv.rules) if indv.rules else None
                if rule is not None:
                    self.regex_mut_replace_subtree(rule)
                    needs_sort = True
            # Change weight
            for rule in indv.rules[:]:
                if R.random() <= self.mCW:
                    self.mut_change_weight(rule)
            # Add amino acids to a leaf node
            if R.random() <= self.mAA:
                rule = R.choice(indv.rules) if indv.rules else None
                if rule is not None:
                    self.regex_mut_add_alphabet(rule)
                    needs_sort = True
            # Replace/invert a node
            if R.random() <= self.mRN:
                rule = R.choice(indv.rules) if indv.rules else None
                if rule is not None:
                    self.regex_mut_replace_node(rule)
                    needs_sort = True
            # Remove from pattern (tree pruning)
            if R.random() <= self.mRFP:
                rule = R.choice(indv.rules) if indv.rules else None
                if rule is not None:
                    self.regex_mut_remove_from_pattern(rule)
                    needs_sort = True
            # Purge rules whose pattern became None after tree mutations
            indv.rules = [r for r in indv.rules if r.pattern]
            if needs_sort:
                indv.bubbleSort()

    # ── Substring-mode mutations ────────────────────────────────────────────

    # Add a random rule mutation
    def mut_add_rule(self, individual):
        if len(individual.rules) >= self.ruleCount:
            return
        pattern = ""
        weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
        # Add these many rules
        for i in range(R.randint(1, self.ruleSize)):
            # Rule size is calculated randomly, and now we need to select a random combination of codes with a
            # specified size
            randomchar = self.codes[R.randint(0, (len(self.codes) - 1))]
            pattern += randomchar
        rule = Rule.Rule(pattern, weight, 0)
        individual.rules.append(rule)

    # Remove rule mutation
    def mut_remove_rule(self, individual):
        if len(individual.rules) == 0:
            return
        tempRand = R.randint(0, len(individual.rules) - 1)
        del individual.rules[tempRand]

    # Add to weight mutation
    def mut_change_weight(self, rule):
        optRand = R.randint(0, 1)
        weightRand = round(R.uniform(self.mCWmin, self.mCWmax), 2)
        if optRand == 0:
            # Addition
            rule.weight += weightRand
        elif optRand == 1:
            # Substraction
            rule.weight -= weightRand

    # Alter patterns mutation (add letter)
    def mut_add_to_pattern(self, rule):
        if len(rule.pattern) >= self.ruleSize:
            return
        pattern = rule.pattern
        randomchar = self.codes[R.randint(0, len(self.codes) - 1)]
        if len(pattern) == 0:
            pattern = randomchar
        else:
            insPos = R.randint(0, len(pattern))
            pattern = pattern[0:insPos] + randomchar + pattern[insPos : (len(pattern))]
        rule.pattern = pattern

    # Alter patterns mutation (remove letter)
    def mut_remove_from_pattern(self, rule):
        if len(rule.pattern) == 0:
            return
        if len(rule.pattern) == 1:
            rule.pattern = ""
            return
        pattern = rule.pattern
        insPos = R.randint(0, len(pattern) - 1)
        pattern = pattern[0:insPos] + pattern[insPos + 1 : (len(pattern))]
        rule.pattern = pattern

    # ── Regex-mode mutations ────────────────────────────────────────────────

    def regex_mut_add_rule(self, individual):
        """Add a new tree-based regex rule."""
        if len(individual.rules) >= self.ruleCount:
            return
        weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
        if self.init_method == "full":
            pattern_re, tree = regex_tree.indi_full(
                self.depth_tree, self.min_braces, self.max_braces
            )
        elif self.init_method == "grow":
            pattern_re, tree = regex_tree.indi_grow(
                self.depth_tree, self.min_braces, self.max_braces
            )
        else:
            pattern_re, tree = regex_tree.indi_half(
                self.depth_tree, self.min_braces, self.max_braces
            )
        if pattern_re is not None:
            rule = Rule.Rule(pattern_re, weight, 0, tree_shape=tree)
            individual.rules.append(rule)

    def regex_mut_replace_rule(self, individual):
        """Replace a random rule with a new tree-based regex."""
        if len(individual.rules) == 0:
            return
        idx = R.randint(0, len(individual.rules) - 1)
        del individual.rules[idx]
        weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
        if self.init_method == "full":
            pattern_re, tree = regex_tree.indi_full(
                self.depth_tree, self.min_braces, self.max_braces
            )
        elif self.init_method == "grow":
            pattern_re, tree = regex_tree.indi_grow(
                self.depth_tree, self.min_braces, self.max_braces
            )
        else:
            pattern_re, tree = regex_tree.indi_half(
                self.depth_tree, self.min_braces, self.max_braces
            )
        if pattern_re is not None:
            rule = Rule.Rule(pattern_re, weight, 0, tree_shape=tree)
            individual.rules.append(rule)

    def regex_mut_replace_subtree(self, rule):
        """Replace a random subtree in the rule's tree with a new one."""
        if rule.tree_shape is None or len(rule.tree_shape) == 0:
            return
        random_node = self._pick_a_node(rule.tree_shape)
        if random_node is None:
            return

        # Find the layer of the selected node
        layer = 0
        for num_layer, nodes in self.dict_layer.items():
            if random_node in nodes:
                layer = num_layer

        if layer == self.last_layer:
            # Leaf node: replace value
            if isinstance(rule.tree_shape[random_node], list):
                if "^" in rule.tree_shape[random_node]:
                    rule.tree_shape[random_node] = R.sample(
                        regex_tree.ALPHABET,
                        R.randint(
                            len(regex_tree.ALPHABET) // 2, len(regex_tree.ALPHABET) - 1
                        ),
                    )
                    rule.tree_shape[random_node].insert(0, "^")
                else:
                    rule.tree_shape[random_node] = R.sample(
                        regex_tree.ALPHABET,
                        R.randint(1, regex_tree.MAX_IN_SQUARE + 1),
                    )
            else:
                rule.tree_shape[random_node] = R.choice(regex_tree.LAST)
        else:
            # Interior node: generate new subtree
            import numpy as np

            tree_method = np.random.choice(["full", "grow"], p=[0.5, 0.5])
            new_depth = (self.depth_tree - layer) + 1
            if tree_method == "full":
                _, new_tree = regex_tree.indi_full(
                    new_depth, self.min_braces, self.max_braces
                )
            else:
                _, new_tree = regex_tree.indi_grow(
                    new_depth, self.min_braces, self.max_braces
                )

            rule.tree_shape[random_node] = new_tree[0]
            # Build child index list for the subtree rooted at random_node
            child = [random_node]
            i = 0
            while (i * 2) + 1 != self.maxnodes:
                if i in child:
                    child.append((i * 2) + 1)
                    child.append((i * 2) + 2)
                i += 1
            # Replace nodes
            for idx, node_index in enumerate(child):
                try:
                    rule.tree_shape[node_index] = (
                        new_tree[idx] if idx < len(new_tree) else None
                    )
                except IndexError:
                    rule.tree_shape[node_index] = None

        rule.pattern = regex_tree.tree2regex(rule.tree_shape)

    def regex_mut_add_alphabet(self, rule):
        """Add 1-4 amino acids to a random leaf node."""
        if rule.tree_shape is None or len(rule.tree_shape) == 0:
            return
        random_node = self._pick_a_node(rule.tree_shape)
        if random_node is None:
            return

        # Find layer
        layer = 0
        for key_layer, nodes in self.dict_layer.items():
            if random_node in nodes:
                layer = key_layer

        if layer == self.last_layer:
            if not isinstance(rule.tree_shape[random_node], list):
                nbr_new_aa = R.randint(1, 4)
                for _ in range(nbr_new_aa):
                    rule.tree_shape[random_node] += R.choice(regex_tree.LAST)
                rule.pattern = regex_tree.tree2regex(rule.tree_shape)

    def regex_mut_replace_node(self, rule):
        """Replace/invert a single operator or character class node."""
        if rule.tree_shape is None or len(rule.tree_shape) == 0:
            return
        random_node = self._pick_a_node(rule.tree_shape)
        if random_node is None:
            return

        node_val = rule.tree_shape[random_node]
        if node_val == "cat":
            rule.tree_shape[random_node] = "|"
        elif node_val == "|":
            rule.tree_shape[random_node] = "cat"
        elif isinstance(node_val, str) and "{" in node_val:
            rule.tree_shape[random_node] = (
                "{" + str(R.randint(self.min_braces, self.max_braces)) + "}"
            )
        elif node_val == "[]":
            rule.tree_shape[random_node] = "[^]"
            child = (random_node * 2) + 1
            if child < len(rule.tree_shape) and isinstance(
                rule.tree_shape[child], list
            ):
                rule.tree_shape[child].insert(0, "^")
        elif node_val == "[^]":
            rule.tree_shape[random_node] = "[]"
            child = (random_node * 2) + 1
            if child < len(rule.tree_shape) and isinstance(
                rule.tree_shape[child], list
            ):
                if "^" in rule.tree_shape[child]:
                    rule.tree_shape[child].remove("^")

        rule.pattern = regex_tree.tree2regex(rule.tree_shape)

    def regex_mut_remove_from_pattern(self, rule):
        """Remove/prune a subtree from the rule's tree structure."""
        if rule.tree_shape is None or len(rule.tree_shape) == 0:
            return
        copysave = copy.deepcopy(rule.tree_shape)

        random_node = self._pick_a_node(rule.tree_shape)
        if random_node is None:
            return

        # Find layer
        layer = 0
        for key_layer, nodes in self.dict_layer.items():
            if random_node in nodes:
                layer = key_layer

        parent = self._parent_of(random_node)
        child = self._subtree_list(random_node)

        try:
            if (
                isinstance(rule.tree_shape[random_node], str)
                and "[" in rule.tree_shape[random_node]
            ):
                rule.tree_shape[parent] = "cat"
            for c in child:
                rule.tree_shape[c] = None
            if rule.tree_shape[parent] == "|":
                rule.tree_shape[parent] = "cat"
            if rule.tree_shape[parent] == "+":
                rule.tree_shape[parent] = None
            if (
                isinstance(rule.tree_shape[parent], str)
                and "{" in rule.tree_shape[parent]
            ):
                rule.tree_shape[parent] = None
                x = self._parent_of(parent)
                if rule.tree_shape[x] == "|":
                    rule.tree_shape[x] = "cat"
            if rule.tree_shape[parent] == "[^]":
                rule.tree_shape[parent] = None
            if rule.tree_shape[parent] == "[]":
                rule.tree_shape[parent] = None
        except (IndexError, TypeError):
            rule.tree_shape = copysave
            rule.pattern = regex_tree.tree2regex(rule.tree_shape)
            return

        rule.pattern = regex_tree.tree2regex(rule.tree_shape)
        if rule.pattern is not None and ("(|" in rule.pattern or "|)" in rule.pattern):
            rule.tree_shape = copysave
            rule.pattern = regex_tree.tree2regex(rule.tree_shape)

    # ── Tree helper methods (regex mode) ────────────────────────────────────

    def _build_layer_dict(self) -> int:
        """Build the dict mapping layer number → list of node indices.

        Returns the last_layer number.
        """
        node = 0
        last_layer = 0
        for j, layer in enumerate(range(self.depth_tree - 1, -1, -1)):
            if j == 0:
                last_layer = layer + 1
            for _ in range(2 ** (self.depth_tree - (layer + 1))):
                node += 1
                key = self.depth_tree - layer
                if key not in self.dict_layer:
                    self.dict_layer[key] = []
                self.dict_layer[key].append(node - 1)
        return last_layer

    def _pick_a_node(self, tree_shape) -> int | None:
        """Pick a random non-None, non-root node from the tree."""
        if len(tree_shape) <= 1:
            return None
        for _ in range(100):  # avoid infinite loops
            idx = R.randint(0, len(tree_shape) - 1)
            if tree_shape[idx] is not None and idx != 0:
                return idx
        return None

    def _parent_of(self, index: int) -> int:
        if index == 0:
            return 0
        if index % 2 == 0:
            return (index - 2) // 2
        return (index - 1) // 2

    def _brother_of(self, index: int) -> int:
        if index == 0:
            return 0
        if index % 2 == 0:
            return index - 1
        return index + 1

    def _subtree_list(self, root_node: int) -> list[int]:
        """Get list of all node indices in the subtree rooted at root_node."""
        nodes = [root_node]
        i = 0
        while (i * 2) + 1 != self.maxnodes:
            if i in nodes:
                nodes.append((i * 2) + 1)
                nodes.append((i * 2) + 2)
            i += 1
        return nodes

    # ── Shared utilities ────────────────────────────────────────────────────

    def removeExtra(self, indv):
        seen = {}
        to_remove = set()
        for i, rule in enumerate(indv.rules):
            if rule.pattern in seen:
                prev_idx = seen[rule.pattern]
                if rule.status == 0:
                    to_remove.add(i)
                elif indv.rules[prev_idx].status == 0:
                    to_remove.add(prev_idx)
                    seen[rule.pattern] = i
                else:
                    to_remove.add(i)
            else:
                seen[rule.pattern] = i
        indv.rules = [r for i, r in enumerate(indv.rules) if i not in to_remove]

import rule as Rule
import random as R
import re
import pandas as pd

import regex_tree
import pattern_engine as PE

_codes_cache = None


def _get_codes():
    global _codes_cache
    if _codes_cache is None:
        _codes_cache = pd.read_csv("data/translation/amino_to_amino.csv")[
            "code"
        ].tolist()
    return _codes_cache


class Individual:
    # Constructor
    def __init__(self, config):
        self.rules = []
        self.usedRules = {}
        self.usedRulesCount = 0
        self.ruleSize = int(config["maximum_rule_size"])
        self.maxRuleCount = int(config["maximum_rule_count"])
        self.minWeight = float(config["rule_weight_min"])
        self.maxWeight = float(config["rule_weight_max"])
        self.fitness = 0
        self.test = 0
        self.extra = {}
        self.matching_mode = config.get("matching_mode", "substring")

    def remove_unexpressed(self):
        self.rules = [rule for rule in self.rules if rule.status != 0]

    def makeFromFile(self, file):
        self.rules = []
        self.usedRules = {}
        self.fitness = 0
        tempIndv = pd.read_csv(file)
        tmpPatterns = tempIndv["pattern"]
        tmpWeights = tempIndv["weight"]
        try:
            tmpStatus = tempIndv["status"]
        except:
            print(
                "Pro-Predictor: model {} does not have status column, putting 0 for all.".format(
                    file
                )
            )
            tmpStatus = [0] * len(tmpPatterns)
        for i in range(len(tmpPatterns)):
            pattern = tmpPatterns[i]
            if not isinstance(pattern, str) or pattern != pattern:
                continue
            rule = Rule.Rule(pattern, tmpWeights[i], tmpStatus[i])
            self.rules.append(rule)

    def init_pattern(self, gap_chance=0.15, config=None):
        """Initialize random substring rules.

        *gap_chance* controls the probability that any one character position
        in a newly created pattern is a gap ('_') instead of a concrete amino
        acid.  Set to 0.0 to disable initial gaps entirely.

        When *config* is provided and experimental features are enabled, the
        element-based engine is used to generate richer initial patterns.
        """
        codes = _get_codes()

        # Detect experimental features from config
        _bool = lambda k, d="False": (
            config.get(k, d).strip().lower() == "true" if config else False
        )
        use_exp = _bool("exp_char_classes") or _bool("exp_variable_gaps")
        use_pw = _bool("exp_weighted_positions")
        class_chance = 0.15 if _bool("exp_char_classes") else 0.0
        _mcs = int(config.get("max_class_size", "0")) if config else 0
        max_class_size = _mcs if _mcs > 0 else None

        for _i in range(R.randint(1, int(self.maxRuleCount / 3))):
            weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
            pat_len = R.randint(1, self.ruleSize)

            if use_exp:
                elements = []
                for _j in range(pat_len):
                    elements.append(
                        PE.random_element(
                            codes, gap_chance, class_chance, max_class_size
                        )
                    )
                # Ensure at least one concrete element
                if PE.is_all_wildcard(elements):
                    idx = R.randint(0, len(elements) - 1)
                    elements[idx] = ("c", R.choice(codes))
                # Gaps are only allowed in the middle - strip leading/trailing
                elements = PE.strip_edge_gaps(elements)
                if not elements:
                    elements = [("c", R.choice(codes))]
                pattern = PE.render_elements(elements)
                rule = Rule.Rule(pattern, weight, 0)
                if use_pw and not PE.has_variable_length(elements):
                    rule.position_weights = [1.0] * len(elements)
            else:
                pattern = ""
                for _j in range(pat_len):
                    if R.random() < gap_chance:
                        pattern += "_"
                    else:
                        pattern += codes[R.randint(0, len(codes) - 1)]
                # Ensure at least one concrete amino acid (no all-gap patterns)
                if all(ch == "_" for ch in pattern):
                    idx = R.randint(0, len(pattern) - 1)
                    pattern = (
                        pattern[:idx]
                        + codes[R.randint(0, len(codes) - 1)]
                        + pattern[idx + 1 :]
                    )
                # Gaps are only allowed in the middle - strip leading/trailing
                pattern = pattern.strip("_")
                if not pattern:
                    pattern = codes[R.randint(0, len(codes) - 1)]
                rule = Rule.Rule(pattern, weight, 0)

            self.rules.append(rule)
        self.bubbleSort()

    def init_regex_pattern(self, config):
        """Initialize rules with tree-based regex patterns (regex mode)."""
        depth = int(config.get("max_depth_tree", "4"))
        min_braces = int(config.get("min_braces", "1"))
        max_braces = int(config.get("max_braces", "3"))
        init_method = config.get("init_method", "half")

        for _ in range(R.randint(1, int(self.maxRuleCount / 4))):
            weight = round(R.uniform(self.minWeight, self.maxWeight), 2)

            if init_method == "half":
                pattern_re, tree = regex_tree.indi_half(depth, min_braces, max_braces)
            elif init_method == "grow":
                pattern_re, tree = regex_tree.indi_grow(depth, min_braces, max_braces)
            elif init_method == "full":
                pattern_re, tree = regex_tree.indi_full(depth, min_braces, max_braces)
            else:
                raise ValueError(
                    "Invalid init_method '{}': expected grow, full, or half".format(
                        init_method
                    )
                )

            if pattern_re is not None and self._check_regex(pattern_re):
                rule = Rule.Rule(pattern_re, weight, 0, tree_shape=tree)
                if rule.tree_shape[0] not in ("cat", "|"):
                    continue
                self.rules.append(rule)
        self.bubbleSort()

    @staticmethod
    def _check_regex(pattern):
        """Validate that a pattern compiles as a regex."""
        if pattern is None:
            return False
        try:
            re.compile(pattern)
            return True
        except re.error:
            return False

    # for i in self.rules:
    # 	print(str(i.pattern) + " => " + str(i.weight))

    def bubbleSort(self):
        self.rules.sort(key=lambda r: len(r.pattern) if r.pattern else 0, reverse=True)

    def print(self):
        for kh, rule in enumerate(self.rules):
            print("{}- {} - {} - {}".format(kh, rule.pattern, rule.weight, rule.status))

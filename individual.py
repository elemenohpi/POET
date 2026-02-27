import rule as Rule
import random as R
import pandas as pd

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

    def init_pattern(self):
        codes = _get_codes()
        for i in range(R.randint(1, int(self.maxRuleCount / 3))):
            pattern = ""
            weight = round(R.uniform(self.minWeight, self.maxWeight), 2)
            for j in range(R.randint(1, self.ruleSize)):
                pattern += codes[R.randint(0, len(codes) - 1)]
            rule = Rule.Rule(pattern, weight, 0)
            self.rules.append(rule)
        self.bubbleSort()

    # for i in self.rules:
    # 	print(str(i.pattern) + " => " + str(i.weight))

    def bubbleSort(self):
        self.rules.sort(key=lambda r: len(r.pattern), reverse=True)

    def print(self):
        for kh, rule in enumerate(self.rules):
            print("{}- {} - {} - {}".format(kh, rule.pattern, rule.weight, rule.status))

# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

import configparser

import random as R

global rule_size, TT, default_learn, default_unseen, population_size, alphabet_size, maximum_rules_count, rule_value_min, rule_value_max

config = configparser.ConfigParser()
config.read('config.ini')

# Parameter Initialization

# Threading
enable_threading = config["THREADING"]["enable_threading"]

# LEARN
default_learn = config['LEARN']['learn_data']	
default_unseen = config['LEARN']['unseen_data']

# PREDICTOR
prediction_model = config['PREDICTOR']['prediction_model']

# TRANSLATION
TT = "data/translation/" + config['TRANSLATION']['translation_table']
dic = {}
alphabet_size = int(config['TRANSLATION']['alphabet_size'])				

# GENETIC PROGRAMMING
population_size = int(config['GP']['population_size'])
rule_size = int(config['GP']['maximum_rule_size'])
maximum_rule_count = int(config['GP']['maximum_rule_count'])
rule_weight_min = float(config['GP']['rule_weight_min'])
rule_weight_max = float(config['GP']['rule_weight_max'])
runs = int(config['GP']['runs'])
pattern_mode = int(config['GP']['pattern_mode'])
tournament_size = int(config['GP']['tournament_size'])
cross_rate = float(config['GP']['crossover_unused_selection_chance'])

# output_name
pop_log_interval = int(config['OUTPUT']['pop_log_interval'])
output_file_name = config['OUTPUT']['output_name']

# GENERAL
seed = int(config['GENERAL']['seed'])
debug = int(config['GENERAL']['debug'])

# MUTATION
mut_add_rule = float(config['MUTATION']['mut_add_rule'])
mut_remove_rule = float(config['MUTATION']['mut_remove_rule'])
mut_change_weight = float(config['MUTATION']['mut_change_weight'])
mut_add_to_pattern = float(config['MUTATION']['mut_add_to_pattern'])
mut_remove_from_pattern = float(config['MUTATION']['mut_remove_from_pattern'])
# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

class Rule:
	def __init__(self, pattern, weight, status, tree):
		self.pattern = pattern
		self.weight = weight
		self.status = status

		self.tree_shape = tree
		self.complexity = 0
		self.score = 0

		self.db_motif = {}


	def clear_pattern(self):
		self.tree_shape = []
		self.pattern = ''

	





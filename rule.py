# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University

class Rule:
	def __init__(self, pattern, weight, status):
		self.pattern = pattern
		self.weight = weight
		self.status = status

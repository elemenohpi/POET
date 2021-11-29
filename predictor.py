# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University
import settings
import fitness as F
import random as R
import individual as I
import copy

class Sequence:
	def __init__(self, pattern, fitness):
		self.pattern = pattern
		self.fitness = fitness
		pass

# all I know is that iter is a lie... #ToDo:: 
class Predictor:

	def __init__(self, count, size, iter):
		self.count = count
		self.size = size
		self.iter = iter
		self.pop = []
		self.output = []
		self.popSize = 10
		self.codes = settings.TT['key']
		self.hydroDic = {
			"I":-0.31,
			"L":-0.56,
			"F":-1.13,
			"V":0.07,
			"M":-0.23,
			"P":0.45,
			"W":-1.85,
			"J":0.17,
			"T":0.14,
			"E":2.02,
			"Q":0.58,
			"C":-0.24,
			"Y":-0.94,
			"A":0.17,
			"S":0.13,
			"N":0.42,
			"D":1.23,
			"R":0.81,
			"G":0.01,
			"H":0.96,
			"K":0.99
		}
		pass

	def predict(self):
		# self.populate()
		# self.sort()
		# for i in range(self.count):
		# 	print(self.pop[i].pattern+" -> "+str(self.pop[i].fitness))
		# pass

		objF = F.Fitness()
		individual = I.Individual()
		individual.makeFromFile(settings.prediction_model)

		print("Model Rules:")
		rules = []
		for rule in individual.rules:
			if rule.status == 1:
				rules.append(rule)

		rules = self.sort_rules(rules)

		for rule in rules:
			print (round(rule.weight, 2), " ", rule.pattern, " ", rule.status)

		self.populate()

		self.sort()
		for i in range(self.count):
			print(repr(i+1)+":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness,2))

		for _ in range(self.iter):
			self.sort()
			# for i in range(self.count):
			# 	print(repr(i+1)+":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness,2))
			
			improvement = False

			# mutate
			for myi, seq in enumerate(self.pop):
				randomsite = R.randint(0, self.size-1)
				# print("randomsite", randomsite)
				rand = R.randint(0, len(self.codes) - 1)
				newseq = copy.copy(seq)
				# print("seq", seq.pattern)
				# print(newseq.pattern)
				# exit()
				newamino = self.codes[rand]
				# print("newamino", newamino)
				newseq.pattern = newseq.pattern[:randomsite] + newamino + newseq.pattern[randomsite+1:]
				# print(seq.pattern)
				# print(newseq.pattern)
				tempF = objF.predict(newseq.pattern, individual)
				# print("tempF", tempF)
				if tempF > seq.fitness:
					newseq.fitness = tempF
					# print(seq.pattern + " -> " + newseq.pattern + " | fitness: " + repr(seq.fitness) + " -> " + repr(newseq.fitness) )
					self.pop[myi] = newseq
					improvement = True

			if improvement:
				for i in range(self.count):
					print(repr(i+1)+":\t", self.pop[i].pattern, "\t", round(self.pop[i].fitness,2))

		
	
	def sort_rules(self, rules):
		sortedR = []
		for rule in rules:
			placed = False
			for index, r in enumerate(sortedR):
				if rule.weight < r.weight:
					sortedR.insert(index, rule)
					placed = True
					break
			if not placed:
				sortedR.append(rule) 
		sortedR = reversed(sortedR)
		return sortedR


		# self.size
		# self.count
		

	def populate(self):
		print ("Populating the sequence pool...\n")
		# Read the best model and turn it to an individual
		individual = I.Individual()
		individual.makeFromFile(settings.prediction_model)

		objF = F.Fitness()
		for i in range(self.popSize):
			pattern = ""
			if self.size <= 0:
				self.size = R.randint(8, 12)
			for j in range(self.size - 1):
				rand = R.randint(0, len(self.codes) - 1)
				pattern += self.codes[rand]
			if self.hydrophobic(pattern):
				tempF = 0
			else:
				tempF = objF.predict(pattern, individual)
			sequence = Sequence(pattern, tempF)
			self.pop.append(sequence)

	def hydrophobic(self, pattern):
		temp = 0
		for char in pattern:
			temp += self.hydroDic[char]
		if temp >= 0:
			return False
		else:
			return True

	def sort(self):
		print ("Sorting the sequences...\n")
		n = len(self.pop)
		# Traverse through all array elements
		for i in range(n):
			# Last i elements are already in place
			for j in range(0, n-i-1):
				# traverse the array from 0 to n-i-1
				# Swap if the element found is greater
				# than the next element
				if self.pop[j].fitness < self.pop[j+1].fitness:
					self.pop[j], self.pop[j+1] = self.pop[j+1], self.pop[j]

# Restrict the amino acid search
# impirical research
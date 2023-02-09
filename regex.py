#!/usr/bin/env python3
# @autor: N. Scalzitti
# @date: 11/01/2022

### LIBRAIRIES ###
import random
import numpy as np
import csv
import re
import eletility

configparser = eletility.ConfigParser()
config = configparser.read("config.ini")

ALPHABET_FILE = config["learn_data"]

### CONSTANTES ###

# TODO: {m,n} implementer ca ?

# OPCODE = { 'A':0, 'C':0, 'G':0, 'T':0, '.':0,
# 		   '[':1, '{':1, '+':1, '^':1, # '*':1,  '?':1
# 		   '|':2, 'cat':2 }

# ALPHABET = ["A", "C", "G", "T"]
# LAST = ['A', 'C', 'G', 'T', '.']
# OPERATOR = ['+', '{', '[', '|', 'cat', '^'] # remove '*' and '?'

def open_file(file):
	'''
	Create the alphabet, leaf node etc for individuals
	file: the file containing the specificity of the alphabet
	'''

	opcode = {} # Dictionnary containing all values for leaf nodes and all operators 
	alphabet = [] # all values for leaf nodes except '.'
	last = [] # all values for leaf nodes
	operator = [] # all operators

	with open(file,'r') as data:
		for line in data:
			line = line.strip().split(';')
			cara = line[0]
			arity = int(line[1])
			opcode[cara] = arity

	for k, v in opcode.items():
		
		if v == 1 or v == 2:
			operator.append(k)

		if v == 0 :
			last.append(k)

		if v == 0 and k != '.':
			alphabet.append(k)

	return opcode, operator, last, alphabet

# alphabet_file = "data/translation/alphabet.csv"
alphabet_file = f"data/Alphabet/{ALPHABET_FILE}.csv" # pptes physico-chimiques des AA

# Create the alphabet and the available operators
OPCODE, OPERATOR, LAST, ALPHABET = open_file(alphabet_file)
MAX_IN_SQUARE = int(len(LAST)/2)



### Individual generation  ###
def indi_grow(depth, min_braces, max_braces):
	'''
	Create a new individual with GROW initialization method (heterogeneous branch)

	depth: 		The depth of the tree shape (size of the regex)
	min_braces: minimal value in braces
	max_braces: maximal value in braces

	return a regular expression
	'''

	indi = Grow(depth, min_braces, max_braces)
	indi.build_grow_tree()

	return tree2regex(indi.tree), indi.tree

def indi_full(depth, min_braces, max_braces):
	'''
	Create a new individual with FULL initialization method (homogeneous branch)

	depth: 		The depth of the tree shape (size of the regex)
	min_braces: minimal value in braces
	max_braces: maximal value in braces

	return a regular expression
	'''

	indi = Full(depth, min_braces, max_braces)
	indi.build_full_tree()

	return tree2regex(indi.tree), indi.tree

def indi_half(depth, min_braces, max_braces):
	'''
	Create a new individual with RAMPED HALF and HALF initialization method (homogeneous branch)

	depth: 		The depth of the tree shape (size of the regex)
	min_braces: minimal value in braces
	max_braces: maximal value in braces

	return a regular expression
	'''

	# Probability of 50% to iniate a individual with the Grow or Full method
	x = np.random.choice([0,1], p=[0.5, 0.5])
	if x == 1:
		return indi_grow(depth, min_braces, max_braces)
	else:
		return indi_full(depth, min_braces, max_braces)


# functions to build regex/individuals
def build_regex(array):
	'''
	Transform the array in readable regex
	array: the array of the regex

	return the regex
	'''

	for i, node in enumerate(array):
		if isinstance(node,list):
			node = ''.join(node)
			array[i] = "[" + node + "]"

	regex = ''.join(array[1:-1]).replace("cat","")

	return regex

def explore_tree(tree, arr, index_node=0):
	'''
	Visit each node og a tree to obtain a readable array
	Inspired by: https://en.wikipedia.org/wiki/Tree_traversal

	tree: the tree shape to explore
	arr: the tree shape is converted in array to obtain a regex
	'''

	# The node contain a value
	if tree[index_node] != None:

		if tree[index_node] in OPERATOR:
			arr.append('(')
		try:
			explore_tree(tree, arr, index_node=(index_node*2)+1)
		except:
			pass

		# Node is [] or [^]
		if isinstance(tree[index_node], list):
			arr.append(tree[index_node])
		else:
			if tree[index_node] == "[]":
				pass
			elif tree[index_node] == "[^]":
				pass
			else:
				 arr.append(tree[index_node])

		try:
			explore_tree(tree, arr, index_node=(index_node*2)+2)
		except:
			pass

		if tree[index_node] in OPERATOR:
			arr.append(')')

def tree2regex(tree):
	'''
	Convert a regex tree shape by an array
	tree: the tree shape of the regex

	return the regex
	'''
	array=[]
	explore_tree(tree, array, index_node=0) # build the array to transform a tree shape in readable regex
	regex = build_regex(array)

	try:
		re.compile(regex)
	except:
		return None

	return regex


# Function to test the construction of individual tree shape
def test():
	'''
	To test the different methods
	'''
	g = Grow(4,1,3)
	g.build_grow_tree()
	print('Grow: ', g.tree)
	regex = tree2regex(g.tree)
	print('Grow: ',regex)

	# f = Full(4,1,3)
	# f.build_full_tree()
	# print('Full: ',f.tree)
	# regex = tree2regex(f.tree)
	# print('Full: ',regex)


			###########################################
			###                                     ###
			###             FULL INDI               ###
			###                                     ###
			###########################################

class Full:

	def __init__(self, depth, min_braces, max_braces):
		self.depth = depth 				 # Depth of the tree
		self.max_nodes = (2**depth)-1	 # number of nodes
		self.tree = [0] * self.max_nodes # General structure of the tree
		self.min_braces = min_braces
		self.max_braces = max_braces
		self.leaves = [] 				 # nodes of the last layer
		self.mid_layer = [] 			 # nodes of the penultimate layer



	def define_leaves_in_tree(self):
		'''
		Define which node belong to which layer
		'''

		# print(type(2**(self.depth-2)), 2**(self.depth-2))
		if self.depth == 1:
			self.leaves.append(0)
		elif self.depth == 0:
			self.leaves.append(0)
		else:
			# Nodes in last layer
			for i in range(2**(self.depth-1)):
				self.leaves.append((self.max_nodes - i)-1)
			
			# Before last layer nodes
			for i in range(2**(self.depth-2)):
				self.mid_layer.append( ((self.max_nodes - 2**(self.depth-1)) - i)-1)

		self.mid_layer = list(reversed(self.mid_layer))
		self.leaves = list(reversed(self.leaves))


	def add_extra_leaf(self, extra, proba=[0.3, 0.7]):
		'''
		Probability of (default) 30% to replace a 'None' leaf by a leaf in LAST list
		Increase the size and complexity of the Regex

		extra: Activate or not the function
		return a random leaf (30%) or a None leaf (70%) is the function is activate (extra=True), otherwise
		return None
		'''
		if extra:
			x = np.random.choice([0,1], p=proba)
			if x == 0:
				return random.choice(LAST)
			else:
				return None
		else:
			return None

	def add_leaf(self, i):

		return random.choice(LAST)

	def add_operator_arity1(self, i, extra=False):
		'''
		Add a new operator with an arity of 1 or 2. That means the node i is in the second to last layer

		i: The current node in the tree
		extra: activate or not the function to extend the regex
		'''

		# Generate a new random node
		self.tree[i] = random.choice(OPERATOR) #['[]',"+",'[^]', 'cat', '|', '{']

		if self.tree[i] == '[':
			self.tree[i] = '[]'
			self.tree[(i*2) + 1] = random.sample(ALPHABET, random.randint(1, MAX_IN_SQUARE+1)) # child1=List
			self.tree[(i*2) + 2] = None # self.add_extra_leaf(extra) # or None

		elif self.tree[i] == '^':
			self.tree[i] = '[^]'
			self.tree[(i*2) + 1] = ['^']
			self.tree[(i*2) + 1]+=random.sample(ALPHABET, random.randint(MAX_IN_SQUARE, len(ALPHABET)-1)) #random.randint(1, len(ALPHABET)-1)
			self.tree[(i*2) + 2] = None # self.add_extra_leaf(extra) # or None

		elif self.tree[i] == '+':
			self.tree[(i*2) + 1] = random.choice(LAST) # ---------> pourrait egalement etre une liste
			self.tree[(i*2) + 2] = None # self.add_extra_leaf(extra) # or None

		elif self.tree[i] == '{':
			x = random.randint(self.min_braces,self.max_braces) 
			self.tree[i] = "{" + str(x) + "}"
			self.tree[(i*2) + 1] = 0
			self.tree[(i*2) + 2] = None # self.add_extra_leaf(extra) # or None
		
		elif self.tree[i] == 'cat':
			self.tree[(i*2) + 1] = random.choice(LAST)
			self.tree[(i*2) + 2] = random.choice(LAST)

		elif self.tree[i] == '|':
			self.tree[(i*2) + 1] = random.choice(LAST)
			self.tree[(i*2) + 2] = random.choice(LAST)
		
	def add_operator_arity2(self, i):
		'''
		Add a new operator with an arity of 2. That means the node i is in an intermediate layer

		i: The current node in the tree
		'''
		return random.choice(['cat', '|'])

	def build_full_tree(self, extra=False):
		'''
		Build the tree shape with FULL method
		extra: activate or not the function to extend the regex

		return the tree in list format
		'''

		self.define_leaves_in_tree()

		# Goes through all the nodes
		for node in range(self.max_nodes):
			if node == 0 :  # Root
				self.tree[node] = random.choice(['cat','|'])

			else:
				if self.tree[node] == 0: # Empty node
					
					if node in self.leaves:      # last layer
						self.tree[node] = self.add_leaf(node)

					elif node in self.mid_layer: # before last layer
						self.add_operator_arity1(node, extra)

					else:
						self.tree[node] = self.add_operator_arity2(node) # other layers

		# print(self.tree)
		return self.tree

			###########################################
			###                                     ###
			###             GROW INDI               ###
			###                                     ###
			###########################################

class Grow:

	def __init__(self, depth, min_braces, max_braces):
		self.depth = depth # Depth of the tree
		self.max_nodes = (2**depth)-1
		self.tree = [0] * self.max_nodes # General structure of the tree
		self.min_braces = min_braces
		self.max_braces = max_braces

		# print(self.depth)


	def new_node(self, i, spe_case=False):
		"""
		Add an operator from the OPCODE list, in the current node i, that is never see by 
		the algorithm

		i: the current node
		return a value in OPCODE list
		"""

		# Root node
		if i == 0:
			return random.choice(["cat","|"])
		# Other node
		else:
			# Avoid the conflict between {} and +
			if spe_case:
				tmp_list = []
				for i in LAST:
					tmp_list.append(i)

				tmp_list.append('cat')
				tmp_list.append('^')
				tmp_list.append('[')

				return random.choice(tmp_list)


			else:
				return random.choice(list(OPCODE.keys()))

	def add_leaf(self, i):
		"""
		Add a value for specific operator like: [], {}, + . The node i take the value ([], {}, +)
		and the child node (i*2) +1 take an other value 
		operator: [], {}, +

		i: index/position of the element in the tree

		return the node value
		"""

		operator = self.tree[i] # The current parent with arrity 1

		# CASE []
		if operator == '[':	
			# Return a k length (min=1) list of unique elements chosen from the ALPHABET. 
			a = random.randint(1, MAX_IN_SQUARE)
			self.tree[(i*2) + 1] = random.sample(ALPHABET, a)
			return "[]" # add current nodes
		
		# CASE [^]
		elif operator == '^':
			self.tree[(i*2) + 1] = ['^']
			self.tree[(i*2) + 1]+=random.sample(ALPHABET, random.randint(MAX_IN_SQUARE, len(ALPHABET)-1))
			return "[^]" # add current nodes

		# CASE {}
		elif operator == "{":
			# Random toss between min and max values
			x = random.randint(self.min_braces,self.max_braces) 
			self.tree[(i*2) + 1] = 1 # add an element as child node

			return "{" + str(x) + "}" # add current node

		# # CASE * (Matches the preceding element zero or more times)
		# elif operator == "*":
		# 	self.tree[(i*2) + 1] = 0 #random.choice(elements)  # add a new child node

		# 	return "*" # add current node
		
		# CASE + (Matches the preceding element one or more times)
		elif operator == "+":
			self.tree[(i*2) + 1] = 1 #random.choice(elements)  # add a new child node

			return "+" # add current node

	def add_child(self, i, arity):
		"""
		Add X child(s), depending of the number of ther arity of parent node

		i: the current node in index i
		arity: the number of child for the current node i
		"""

		# Node with 0 child
		if arity == 0: 
			try:
				self.tree[(i*2) + 1] = None
				self.tree[(i*2) + 2] = None
			except IndexError:
				pass

		# Node with 1 child. 4 types of leaves with 1 child is possible: [], {}, + and *
		elif arity == 1: 
			try:
				self.tree[i] = self.add_leaf(i)
				self.tree[(i*2) + 2] = None
			except IndexError:
				pass

		# Node with 2 children.
		elif arity == 2: 
			try:
				self.tree[(i*2) + 1] = 0
				self.tree[(i*2) + 2] = 0
			except IndexError:
				pass

	def check_last_layer(self):
		"""
		Check if leaves (nodes in the last layer) are composed by only LAST values,
		such as: Amino acids, nucleic acids, or joker (.). Operators like {}, | etc
		can't be used in leaf node
		"""
		# Extract all leaves node (in the last layer)
		for i in range(2**(self.depth-1)):
			x = random.choice([0,1]) # avec cette ligne on exploite le bug qui cree des regex interessante
			x=0
			if x == 1:
				# If a leaf is a list, such as: [A, G] it is a good leaf
				if isinstance(self.tree[-(i+1)], list):
					pass
				# If the leaf value is not in LAST list, then change the leaf value
				# Si la profondeur =3 ou 4 par ex, alors ca genere des motifs interessant
				# ex: (GTC.C.G)(T+)
				elif self.tree[-(i+1)] not in LAST or self.tree[-(i+1)] != None: # BUG !!!!!!!!!!!!!!!!!
					self.tree[-(i+1)] = random.choice(LAST)
			else:
				# If a leaf is a list, such as: [A, G] it is a good leaf
				if isinstance(self.tree[-(i+1)], list):
					pass
				elif self.tree[-(i+1)] not in LAST and self.tree[-(i+1)] != None:
					self.tree[-(i+1)] = random.choice(LAST)
	
	def build_grow_tree(self):
		"""
		Builds a tree representing an individual. The construction method of the tree is GROW
		The structure is as follows:
		node 0 = root 
		node i = parent
		node (i*2)+1 = child_1 of i node
		node (i*2)+2 = child_2 of i node

		Ex: parent = node 2, child_1 = (2*2)+1=5, child_2 = (2*2)+2 = 6
		Ex: [root, p1, p2, c1p1,c2p1, c1p2, c2p2 ]

		return a list representing the individual in a tree shape
		"""

		# nbr_node = ((2**self.depth) -1) # Nombre de noeud dans la couche
		# print(nbr_node, self.max_nodes)

		# Goes through all the nodes
		for i in range(self.max_nodes):

			# New node case (value 0=never traversed by the algorithm)
			if self.tree[i] == 0:
				# Modify the current i node with an Operator
				self.tree[i] = self.new_node(i) 
				# Modify the child nodes, depending of the value of the Operator

				self.add_child(i, OPCODE[self.tree[i]]) 
				# print(self.tree) # temporary tree
			elif self.tree[i] == 1:
				self.tree[i] = self.new_node(i, spe_case=True) 
				self.add_child(i, OPCODE[self.tree[i]]) 


			# Node with None value, means that parent node is a leaf (without children) 
			elif self.tree[i] == None:
				try:
					self.tree[(i*2) + 1] = None
					self.tree[(i*2) + 2] = None
				except IndexError:
					pass

			# All other cases
			else:
				self.add_child(i, 0)

		# Check if the last layer (leaf) has good value
		self.check_last_layer()

		return self.tree

	
# test()
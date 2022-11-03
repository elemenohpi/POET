print('import Regex ok')
import random
import numpy as np
'''
def create_regex(size_regex, min_braces, max_braces):
	indi = RegexPattern(size_regex, min_braces, max_braces)
	indi.build_grow_tree()
	print(indi.tree)

	return indi.tree2regex()
'''	
def indi_grow(size_regex, min_braces, max_braces):
	indi = Grow(size_regex, min_braces, max_braces)
	indi.build_grow_tree()
	print(indi.tree)

	return indi.tree2regex()


def indi_full(size_regex, min_braces, max_braces):
	indi = Full(size_regex, min_braces, max_braces)
	indi.build_full_tree()
	print(indi.tree)

	return indi.tree2regex()

def indi_half(size_regex, min_braces, max_braces):
	x = np.random.choice([0,1], p=[0.5, 0.5])
	if x == 1:
		print('Type Grow')
		return indi_grow(size_regex, min_braces, max_braces)
	else:
		print('Type Full')
		return indi_full(size_regex, min_braces, max_braces)

# GROW: brnaches avec profondeurs variables
# FULL: brnaches completes, on va jusau'a la profondeur max
# RAMPED HALF-AND-HALF: 



OPCODE = {'A':0, 'C':0, 'G':0, 'T':0, '.':0,
		   '[':1 ,"{":1, "*":1,'^':1, #"+":1
		  '|':2, 'cat':2}

ALPHABET = ["A", "C", "G", "T"]
LAST = ['A', 'C', 'G', 'T', '.']
OPERATOR = ['+', '{', '[', '*', '|', 'cat', '+', '^']


class Full:

	def __init__(self, depth, min_braces, max_braces):
		self.depth = depth # Depth of the tree
		self.max_nodes = (2**depth)-1
		self.tree = [0] * self.max_nodes # General structure of the tree
		self.min_braces = min_braces
		self.max_braces = max_braces
		self.leaves = []
		self.mid_layer = []

	def define_leaves_in_tree(self):
		# Nodes in last layer
		for i in range(2**(self.depth-1)):
			self.leaves.append(self.max_nodes - i)
		
		# Before last layer nodes
		for i in range(2**(self.depth-2)):
			self.mid_layer.append((self.max_nodes - 2**(self.depth-1)) - i)

	def add_extra_leaf(self, proba=[0.3, 0.7]):
		'''
		Probability of (default) 30% to replace a None leaf by a leaf in LAST list
		'''
		x = np.random.choice([0,1], p=proba)
		if x == 0:
			return random.choice(LAST)
		else:
			return None

	def add_leaf(self, i):
		return random.choice(LAST)

	def add_operator_arity1(self, i):
		self.tree[i] = random.choice(['[]',"*",'[^]', 'cat', '|', '{'])
		
		if self.tree[i] == '[]':
			self.tree[(i*2) + 1] = random.sample(ALPHABET, random.randint(1, len(ALPHABET)-1))
			self.tree[(i*2) + 2] = self.add_extra_leaf() # or None

		elif self.tree[i] == '*':
			self.tree[(i*2) + 1] = random.choice(LAST)
			self.tree[(i*2) + 2] = self.add_extra_leaf() # or None

		elif self.tree[i] == '[^]':
			self.tree[(i*2) + 1] = ['^']
			self.tree[(i*2) + 1]+=random.sample(ALPHABET, random.randint(1, len(ALPHABET)-1))
			self.tree[(i*2) + 2] = self.add_extra_leaf() # or None

		elif self.tree[i] == 'cat':
			self.tree[(i*2) + 1] = random.choice(LAST)
			self.tree[(i*2) + 2] = random.choice(LAST)
		elif self.tree[i] == '|':
			self.tree[(i*2) + 1] = random.choice(LAST)
			self.tree[(i*2) + 2] = random.choice(LAST)

		elif self.tree[i] == '{':
			x = random.randint(self.min_braces,self.max_braces) 
			self.tree[i] = "{" + str(x) + "}"
			self.tree[(i*2) + 1] = 0
			self.tree[(i*2) + 2] = self.add_extra_leaf() # or None

		

	def add_operator_arity2(self, i):
		return random.choice(['cat', '|'])

	def build_full_tree(self):
		self.define_leaves_in_tree()

		# Goes through all the nodes
		for node in range(self.max_nodes):
			if node == 0 :  # Root
				self.tree[node] = random.choice(['cat','|'])

			else:
				if self.tree[node] == 0:
					if node+1 in self.leaves:
						self.tree[node] = self.add_leaf(node)
		
					elif node+1 in self.mid_layer:
						self.add_operator_arity1(node)
	
					else:
						# add operators	
						self.tree[node] = self.add_operator_arity2(node)

			'''
			# New node case (value 0=never traversed by the algorithm)
			if self.tree[i] == 0:
				# Modify the current i node with an Operator
				self.tree[i] = self.new_node(i) 
				# Modify the child nodes, depending of the value of the Operator
				self.add_child(i, OPCODE[self.tree[i]]) 
				# print(self.tree) # temporary tree
	
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
			'''
		print(self.tree)
		return self.tree


	def check_last_layer(self):
		"""
		Check if leaves (nodes in the last layer) are composed by only LAST values,
		such as: Amino acids, nucleic acids, or joker (.). Operators like {}, | etc
		can't be used in leaf node
		"""

		# Extract all leaves node (in the last layer)
		for i in range(2**(self.depth-1)):
			# If a leaf is a list, such as: [A, G] it is a good leaf
			if isinstance(self.tree[-(i+1)], list):
				pass
			# If the leaf value is not in LAST list, then change the leaf value

			# elif self.tree[-(i+1)] not in LAST or self.tree[-(i+1)] != None: # BUG !!
				# print('ICI', self.tree[-(i+1)])
			elif self.tree[-(i+1)] not in LAST and self.tree[-(i+1)] != None:
				self.tree[-(i+1)] = random.choice(LAST)


	def build_regex(self, arr):
		for i, node in enumerate(arr):
			if isinstance(node,list):
				node = ''.join(node)
				arr[i] = "[" + node + "]"
		es = ''.join(arr[1:-1]).replace("cat","")
		return es

	def tree2regex(self):
		array=[]
		self.explore_tree(self.tree, array, index_node=0)

		regex = self.build_regex(array)
	
		return regex

	'''
	def IsOperator(self, s):
		"""
		Vérifie si le caractère est un opérateur.
		"""
		if '+' in s:
			return True
		if '{' in s:
			return True
		if '[' in s:
			return True
		if '*' in s:
			return True
		if '|' in s:
			return True
		if 'cat' in s:
			return True
		if '+' in s:
			return True
		return False
	'''

	def explore_tree(self, tree, arr, index_node=0):

		# The node contain a value
		if tree[index_node] != None:

			if tree[index_node] in OPERATOR:
				arr.append('(')
			# if self.IsOperator(tree[index_node]): # Si operator on ouvre les (
			# 	arr.append('(')
	
			try:
				self.explore_tree(tree, arr, index_node=(index_node*2)+1)
			except:
				pass
	
			# if tree[index_node] == "[]":
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
				self.explore_tree(tree, arr, index_node=(index_node*2)+2)
			except:
				pass

			if tree[index_node] in OPERATOR:
			# if self.IsOperator(tree[index_node]):
				arr.append(')')





class Grow:

	def __init__(self, depth, min_braces, max_braces):
		self.depth = depth # Depth of the tree
		self.max_nodes = (2**depth)-1
		self.tree = [0] * self.max_nodes # General structure of the tree
		self.min_braces = min_braces
		self.max_braces = max_braces


	def new_node(self, i):
		"""
		Add an operator from the OPCODE list, in the current node i, that is never see by 
		the algorithm

		i: the current node
		"""
		# If the node is the root, add the concatenate node (+)
		if i == 0:
			return random.choice(["cat","|"])
		# Other node
		else:
			return random.choice(list(OPCODE.keys()))


	def add_leaf(self, i):
		"""
		Add a value for specific operator like: [], {}, + or *. The node i take the value ([], {}, + or *)
		and the child node (i*2) +1 take an other value 
		operator: [], {}, + or *
		i: index/position of the element in the tree
		"""

		operator = self.tree[i] # The current parent with arrity 1

		# CASE [] and [^]
		if operator == '[':	
			# Return a k length (min=1) list of unique elements chosen from the ALPHABET. 
			self.tree[(i*2) + 1] = random.sample(ALPHABET, random.randint(1, len(ALPHABET)-1))
			return "[]" # add current nodes

		elif operator == '^':
			self.tree[(i*2) + 1] = ['^']

			# Return a k length (min=1) list of unique elements chosen from the ALPHABET. 
			self.tree[(i*2) + 1]+=random.sample(ALPHABET, random.randint(1, len(ALPHABET)-1))
			return "[^]" # add current nodes


		# CASE {}
		elif operator == "{":
			# Random toss between min and max values
			x = random.randint(self.min_braces,self.max_braces) 
			self.tree[(i*2) + 1] = 0 # add an element as child node

			return "{" + str(x) + "}" # add current node

		# CASE * (Matches the preceding element zero or more times)
		elif operator == "*":
			self.tree[(i*2) + 1] = 0 #random.choice(elements)  # add a new child node

			return "*" # add current node
		
		# CASE + (Matches the preceding element one or more times)
		elif operator == "+":
			self.tree[(i*2) + 1] = 0 #random.choice(elements)  # add a new child node

			return "+" # add current node

	def add_child(self, i, arity):
		"""
		Add X child(s), depending of the number of child (arity)

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
			# If a leaf is a list, such as: [A, G] it is a good leaf
			if isinstance(self.tree[-(i+1)], list):
				pass
			# If the leaf value is not in LAST list, then change the leaf value

			# elif self.tree[-(i+1)] not in LAST or self.tree[-(i+1)] != None: # BUG !!!!!!!!!!
				# print('ICI', self.tree[-(i+1)])
			elif self.tree[-(i+1)] not in LAST and self.tree[-(i+1)] != None:
				self.tree[-(i+1)] = random.choice(LAST)

	


	def build_grow_tree(self):
		"""
		Builds a tree representing an individual. The type of the tree is GROW
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

	def build_regex(self, arr):
		for i, node in enumerate(arr):
			if isinstance(node,list):
				node = ''.join(node)
				arr[i] = "[" + node + "]"
		es = ''.join(arr[1:-1]).replace("cat","")
		return es

	def tree2regex(self):
		array=[]
		self.explore_tree(self.tree, array, index_node=0)

		regex = self.build_regex(array)
	
		return regex

	'''
	def IsOperator(self, s):
		"""
		Vérifie si le caractère est un opérateur.
		"""
		if '+' in s:
			return True
		if '{' in s:
			return True
		if '[' in s:
			return True
		if '*' in s:
			return True
		if '|' in s:
			return True
		if 'cat' in s:
			return True
		if '+' in s:
			return True
		return False
	'''

	def explore_tree(self, tree, arr, index_node=0):

		# The node contain a value
		if tree[index_node] != None:

			if tree[index_node] in OPERATOR:
				arr.append('(')
			# if self.IsOperator(tree[index_node]): # Si operator on ouvre les (
			# 	arr.append('(')
	
			try:
				self.explore_tree(tree, arr, index_node=(index_node*2)+1)
			except:
				pass
	
			# if tree[index_node] == "[]":
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
				self.explore_tree(tree, arr, index_node=(index_node*2)+2)
			except:
				pass

			if tree[index_node] in OPERATOR:
			# if self.IsOperator(tree[index_node]):
				arr.append(')')






# if __name__ == '__main__':

# 	# i = RegexPattern(4,1,5)
# 	# i.tree = ['{', 'cat', '|', '^', 'cat', '|', 'cat', ['^', 'G'], None, 'T', 'T', 'G', '.', 'A', '.']
# 	# print()
# 	# print(i.tree2regex())
# 	# exit()
# 	a = Full(5,1,3)
# 	a.build_full_tree()
# 	print(a.tree2regex())


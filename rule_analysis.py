import pandas as pd
import sys

model_path = sys.argv[1]
model_number = sys.argv[2]

filename = model_path+"/model_"+model_number+".csv"
print(filename)
model=pd.read_csv(filename, index_col=0).dropna(axis = 0)
#print(model)

pattern = model['pattern']
unique_chars = [p.lower() for p in list(set(''.join(pattern)))]
pattern = [" ".join(list(p.lower())) for p in pattern]
print(pattern)
print(unique_chars)

from sklearn.feature_extraction.text import TfidfVectorizer
vectorizer = TfidfVectorizer(token_pattern=r'[^\s+]')
#vectorizer.fit(unique_chars)
X = vectorizer.fit_transform(pattern)
vocab = vectorizer.vocabulary_

#print(X.toarray()[0])
print(vocab)
#print(pattern[0])

from sklearn.decomposition import PCA

pca = PCA(n_components = 2)
tr = pca.fit_transform(X.toarray())
xy = pca.components_[0:2, :]

from matplotlib import pyplot as plt
print(pca.n_components_)
print(pca.n_samples_)
print(pca.explained_variance_ratio_)
print(xy)
print(tr[0])
markers = sorted(vocab)
#print(markers)
fig, ax = plt.subplots()
ax.scatter(xy[0, :], xy[1, :], marker='none')
ax.set_title("PCA Components")
ax.set_xlabel("Component 1")
ax.set_ylabel("Component 2")
for i in range(xy.shape[1]):
	x = xy[0, i]
	y = xy[1, i]
	ax.annotate(markers[i], (x,y), color='black', fontsize='large', horizontalalignment='center', verticalalignment='center', weight='heavy')
	
plt.savefig(f'{model_path}/pca.png')

import numpy as np
fig, ax = plt.subplots(figsize=(10, 10))
rules = list(model['pattern'])
weight = np.array(model['weight']).reshape(-1, 1)
color = ['red' if p < 0 else 'black' for p in weight]
#print(color)
#print(list(rules))
#print(weight)
from sklearn.preprocessing import StandardScaler, Normalizer, MinMaxScaler
scaler = MinMaxScaler(feature_range=(0, 1))

normalized_weight = scaler.fit_transform(np.abs(weight))*25
#print(normalized_weight)
ax.scatter(tr[:, 0], tr[:, 1], marker = 'none')
for i in range(tr.shape[0]):
	x = tr[i, 0]
	y = tr[i, 1]
	#print(x, y)
	ax.annotate(rules[i], (x,y), color=color[i], fontsize=normalized_weight[i], horizontalalignment='center', verticalalignment='center', weight = 'heavy')
ax.set_title("Pattern Clustering")
ax.set_xlabel("Component 1")
ax.set_ylabel("Component 2")

plt.savefig(f'{model_path}/pca_transform.png')

#print(xy)


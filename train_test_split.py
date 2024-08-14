import pandas as pd
from sklearn.model_selection import train_test_split

name = 'EPG_POET'
full = pd.read_csv(f'data/{name}.csv')
train,test = train_test_split(full, test_size = 0.3, random_state = 110)
train.to_csv('data/EPG_POET_Train.csv', index = False)
test.to_csv('data/EPG_POET_Test.csv', index = False)

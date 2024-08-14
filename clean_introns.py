import pandas as pd
import sys

model_path = sys.argv[1]
model_number = sys.argv[2]

from pathlib import Path
Path((model_path+"/cleaned")).mkdir(parents=True, exist_ok=True)

filename = model_path+"/model_"+model_number+".csv"
print(filename)
model=pd.read_csv(filename, index_col=0)
#print(model[model['status')
cleaned_model = model[model['status'] == 1].drop_duplicates()
cleaned_model = cleaned_model.round(2)
cleaned_model.to_csv((model_path+"/cleaned/model_"+model_number+".csv"))

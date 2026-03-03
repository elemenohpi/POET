import os
import eletility
import fitness as F
import individual as I
import copy
import pandas as pd
from predictor import Sequence


configparser = eletility.ConfigParser()
config = configparser.read("config.ini")

model_path = "output/Mar-02-2026-upar_epoch_0/models/cleaned"
objF = F.Fitness(config)
files = [f for f in os.listdir(model_path) if os.path.isfile(os.path.join(model_path, f))]
ensemble = []
for model in files:
    if model.split(".")[-1] != "csv":
        continue
    individual = I.Individual(config)
    individual.makeFromFile(os.path.join(model_path, model))
    # individual.remove_unexpressed()
    ensemble.append(individual)


learn_data = pd.read_csv(config["learn_data"])
data_set_rules = learn_data["sequence"].tolist()

pop = []
for seq in data_set_rules:
    fitness = 0
    for model in ensemble:
        _, tempF = objF.eval(seq, 0, model, True)
        fitness += tempF
    fitness = fitness / len(ensemble)
    pop.append(Sequence(seq, fitness))

print("Sequence,Fitness")
for i, seq in enumerate(pop):
    print(f'{seq.pattern}, {seq.fitness}')


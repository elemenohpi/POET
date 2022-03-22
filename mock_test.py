import fitness as F
import individual as I
import pandas as pd
import settings
import matplotlib.pyplot as plt
import numpy as np

# 10 (learn8) and 35 (mock) and 31 (overall)
model = "best_epoch_models\E8.csv"
# model = "..\..\..\poetnewres\myfile\POET1\output\model\model_15.csv"
dataset = "data/mock.csv"

Indv = I.Individual()
Indv.makeFromFile(model)

Evaluator = F.Fitness()

df = pd.read_csv(dataset)
df = df.sort_values("fitness")

predictions = []

for seq in df["sequence"]:
    # error, prediction = Evaluator.eval(seq, 0, Indv, True)
    error, prediction = Evaluator.eval(seq, 0, Indv, True)

    couple = [seq, prediction]
    predictions.append(couple)

prediction_df = pd.DataFrame (predictions, columns=["sequence", "prediction"])
prediction_df = prediction_df.sort_values("prediction")
# print(prediction_df)
# print(df)
# exit()

predicted_ranks = []
actual_ranks = []

score = 0
print("Predicted Rank \t Actual Rank \t Sequence \t Error")
for index, data in enumerate(prediction_df["sequence"]):
    predicted_rank = str(index+1)
    sequence = data
    actual_rank = df[df["sequence"]==sequence].index.values[0] + 1
    predicted_ranks.append(int(predicted_rank))
    actual_ranks.append(actual_rank)
    if index < 10 and actual_rank <= 10:
        score += 1

    error = abs(int(predicted_rank) - actual_rank)
    print(predicted_rank, "\t\t", actual_rank, "\t\t", sequence, "\t\t", error)

print("Model score is: ", score)

pre_series = pd.Series(predicted_ranks)
act_series = pd.Series(actual_ranks)

corr = pre_series.corr(act_series)

print (corr)

# plt.style.use('fast')
#
# # plot
# fig, ax = plt.subplots()
#
# plt.xlabel("Protein ID")
# plt.ylabel("Predicted Rank")
# plt.title("Ranking of the Proteins of the Mock Set Using the Best Model of Epoch 8")
#
# x = list(range(43))
#
# # ax.fill_between(gens, q25s, q75s, alpha=.5, linewidth=0)
# # ax.plot(gens, medians, linewidth=2.0)
# ax.bar(x, predicted_ranks, width=1, edgecolor="white", linewidth=0.7, alpha=0.5)
# ax.bar(x, actual_ranks, width=1, edgecolor="white", linewidth=0.7, alpha=0.5)
# plt.show()
# # plt.savefig("./gen_evo/fig.png", dpi=300)


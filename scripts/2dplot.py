import math
import statistics
from os import listdir
from os.path import isfile, join

from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import pandas as pd


def polygon_under_graph(x, y):
    """
    Construct the vertex list which defines the polygon filling the space under
    the (x, y) line graph. This assumes x is in ascending order.
    """
    return [(x[0], 0.), *zip(x, y), (x[-1], 0.)]

gens = 10000

epochs_means = []

for i in range(1, 9):
    path = "D:\\POETPaperPeerJRes\\poet" + str(i) + "res/output/evolution/"
    # path = "D:\\POET_Correlation_April5\\" + str(i) + "/output/evolution/"
    files = [f for f in listdir(path) if isfile(join(path, f))]
    dataframes = []
    for file in files:
        file_tokens = file.split(".")
        if file_tokens[len(file_tokens) - 1] != "log":
            continue
        file = path + file
        # print(file)
        df = pd.read_csv(file, index_col=0)
        if df.size < 80000:
            continue
        dataframes.append(df["best fitness"].values)
    # exit()
    medians = []
    # print(len(dataframes))
    # exit()

    x = range(gens)

    for k in range(gens):
        median_list = []
        for j in range(len(dataframes)):
            median_list.append(dataframes[j][k])
        median = statistics.median(median_list)
        # median = statistics.mean(median_list)

        # if i == 1:
        #     median /= 42
        # elif i == 2:
        #     median /= 51
        # elif i == 3:
        #     median /= 61
        # elif i == 4:
        #     median /= 71
        # elif i == 5:
        #     median /= 82
        # elif i == 6:
        #     median /= 92
        # elif i == 7:
        #     median /= 102
        # elif i == 8:
        #     median /= 112

        # print(median)

        # medians.append(abs(math.log(median)))
        medians.append(median)

    epoch_means = polygon_under_graph(x, medians)
    epochs_means.append(epoch_means)

verts = epochs_means

facecolors = plt.colormaps['viridis_r'](np.linspace(0, 1, len(verts)))

ax = plt.figure().add_subplot(projection='3d')
lambdas = range(1, 9)

poly = PolyCollection(verts, facecolors=facecolors, alpha=.7)
ax.add_collection3d(poly, zs=lambdas, zdir='y')

ax.set(xlim=(0, gens), ylim=(1, 8), zlim=(0, 12),
       xlabel='Generations', ylabel=r'Epoch', zlabel='Error (Training RMSE)')

plt.show()

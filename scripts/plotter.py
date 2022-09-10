import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np


class Plotter:
	# Takes directories, returns a plot of average with q75 and q25. Can plot multiple directories to compare.
	def plot_experiment(self, paths, gens, title, xlabel, ylabel, line_labels, output=None):
		plt.style.use('fast')

		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.title(title)

		all_directory_based_data = []
		for index, path in enumerate(paths):
			q25s = []
			q75s = []
			medians = []
			dataframes = []
			files = [f for f in listdir(path) if isfile(join(path, f))]
			for file in files:
				extension = file.split(".")[-1]
				if extension != "log":
					continue
				file = join(path, file)
				df = pd.read_csv(file, index_col=0)
				if len(df.index) >= gens:
					dataframes.append(df)
				else:
					print("data frame {} has {} rows and thus is disregarded".format(file, len(df.index)))
			# print(len(df.index), file)
			# exit()
			for i in range(gens):
				best_fitness_values_at_gen_i = []
				for df in dataframes:
					# print(df)
					# exit()
					try:
						best_fitness_values_at_gen_i.append(df["best fitness"][i])
					except:
						print(df)
						exit()

				best_fitness_values_at_gen_i_df = pd.DataFrame(best_fitness_values_at_gen_i)
				quantiles = best_fitness_values_at_gen_i_df.quantile([0.25, 0.75])
				medians.append(best_fitness_values_at_gen_i_df.median().values[0])
				q25s.append(quantiles[0].iloc[0])
				q75s.append(quantiles[0].iloc[1])

			directory_data = [line_labels[index], q25s, q75s, medians]
			all_directory_based_data.append(directory_data)

		x_axis_data = range(0, gens)

		# colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
		colors = ["#0570b0", "#3690c0", "#74a9cf", "#cc4c02", "#fe9929", "#dd3497", "#ae017e", "#7a0177"]
		for data_index, line_data in enumerate(all_directory_based_data):
			plt.fill_between(x_axis_data, line_data[1], line_data[2], alpha=.1, linewidth=0)
			plt.plot(x_axis_data, line_data[3], linewidth=2.0, label=line_data[0], color=colors[data_index])
		plt.legend(loc='upper right')
		plt.show()
# plt.savefig(output, dpi=300)


plotter = Plotter()

plotter.plot_experiment(
	["D:\\POETPaperPeerJRes\\poet1res/output/evolution/", "D:\\POETPaperPeerJRes\\poet2res/output/evolution/",
	 "D:\\POETPaperPeerJRes\\poet3res/output/evolution/", "D:\\POETPaperPeerJRes\\poet4res/output/evolution/",
	 "D:\\POETPaperPeerJRes\\poet5res/output/evolution/", "D:\\POETPaperPeerJRes\\poet6res/output/evolution/",
	 "D:\\POETPaperPeerJRes\\poet7res/output/evolution/", "D:\\POETPaperPeerJRes\\poet8res/output/evolution/"],
	10000, "Evolution of POET Models in Each Epoch", "Generations", "Error (Training RMSE)", ["E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"])

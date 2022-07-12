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
				if extension != "csv":
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
						best_fitness_values_at_gen_i.append(df[" fitness"][i])
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

		for line_data in all_directory_based_data:
			plt.fill_between(x_axis_data, line_data[1], line_data[2], alpha=.1, linewidth=0)
			plt.plot(x_axis_data, line_data[3], linewidth=2.0, label=line_data[0])
		plt.legend(loc='lower right')
		plt.show()
	# plt.savefig(output, dpi=300)




plotter = Plotter()




# ================================================================== Static Loop (solving a^20) ========================
# plotter.plot_experiment(
# 	["Output/29apr-statemprmse/Evo", "Output/29apr-staprogrmse/Evo"],
# 	1000, "x^20 Problem with RMSE", "Generations", "Fitness", ["rmse_p", "rmse_s"])
#
# plotter.plot_experiment(
# 	["Output/29apr-staprog/Evo", "Output/29apr-statemp/Evo"],
# 	50, "x^20 Problem with Correlation", "Generations", "Fitness", ["correl_p", "correl_s"])

# ================================================================== Dyna Loop (solving 2^a) ===========================

# plotter.plot_experiment(
# 	["Output/29apr-dynaprog/Evo", "Output/29apr-dynatemp/Evo"],
# 	50, "2^x Problem with Correlation", "Generations", "Fitness", ["correl_p", "correl_s"])
#
# plotter.plot_experiment(
# 	["Output/29apr-dynaprogrmse/Evo", "Output/29apr-dynatemprmse/Evo"],
# 	1000, "2^x Problem with RMSE", "Generations", "Fitness", ["rmse_p", "rmse_s"])

# ================================================================== Dyna Loop (solving a^b) ===========================
#
# plotter.plot_experiment(
# 	["Output/29apr-powerprog/Evo", "Output/29apr-powertemp/Evo"],
# 	1000, "x^y Problem with Correlation", "Generations", "Fitness", ["correl_p", "correl_s"])

# plotter.plot_experiment(
# 	["Output/29apr-powerprogrmse/Evo", "Output/29apr-powertemprmse/Evo"],
# 	1000, "x^y Problem with RMSE", "Generations", "Fitness", ["rmse_p", "rmse_s"])

# ================================================================== RL + Loop (Ant Problem) ===========================

# plotter.plot_experiment(
# 	["Output/29apr-ant-prog/Evo", "Output/29apr-ant-prog-const/Evo", "Output/29apr-ant-prog-const-reg/Evo"],
# 	999, "Comparing Different QGP Settings for Solving the Ant Problem", "Generations", "Fitness", ["prog", "prog-const", "prog-const-reg"])

# plotter.plot_experiment(
# 	["Output/29apr-ant-spatio/Evo", "Output/29apr-ant-spatio-const/Evo", "Output/29apr-ant-spatio-const-reg/Evo"],
# 	999, "Comparing Different QGP Settings for Solving the Ant Problem", "Generations", "Fitness", ["temporo", "temporo-const", "temporo-const-reg"])

# plotter.plot_experiment(
# 	["Output/29apr-ant-prog/Evo", "Output/29apr-ant-prog-const/Evo", "Output/29apr-ant-prog-const-reg/Evo", "Output/29apr-ant-spatio/Evo", "Output/29apr-ant-spatio-const/Evo", "Output/29apr-ant-spatio-const-reg/Evo"],
# 	999, "Comparing Different Settings for Solving the Ant Problem", "Generations", "Fitness", ["prog", "p_const", "p_const_math", "spatial", "s_const", "s_const_math"])

# ================================================================== State/Decision Making (TicTacToe) ===========================

plotter.plot_experiment(
	["Output/TicTacToe/Evo"],
	500, "Evolution of Models for Solving the TicTacToe Problem", "Generations", "Fitness", ["prog"])




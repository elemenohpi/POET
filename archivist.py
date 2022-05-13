import os
import pandas as pd
import datetime
import shutil
import sys
import eletility

class Archivist:
	def __init__(self, config):
		self.config = config
		pass

	def saveCSV(self, df, path, filename):
		try:
			if not os.path.exists(path):
				os.makedirs(path)
			if os.path.exists(path + "/" + filename + ".csv"):
				os.remove(path + "/" + filename + ".csv")
		except:
			print("P-Predictor: Could not initialize the folders correctly")
		try:
			df.to_csv(path + "/" + filename + ".csv")
		except:
			print("P-Predictor: Could not save the file: " + path + "/" + filename + ".csv")

	def save(self, filename, string):
		file = open(filename, "a")

		file.write(string + "\n")

		file.close()

	def saveEvo(self, string):
		path = self.config["output_evo"]
		file = open(path, "a")
		file.write("{}\n".format(string))
		file.close()

	def saveModel(self, df):
		path = self.config["output_model"]
		df.to_csv(path)

	def setup(self, archive=False):
		# Ensure the appropriate directories exist. Remove any existing log file and archive them using appropriate time/date
		outputDir = "./output"
		archiveDir = "./archive"

		# Check if everything is in place

		if not os.path.exists(outputDir):
			os.makedirs(outputDir)
		if not os.path.exists(archiveDir):
			os.makedirs(archiveDir)

		F = eletility.Files()
		F.writeTruncate(self.config["output_evo"], "")
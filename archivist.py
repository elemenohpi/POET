import os
import pandas as pd
import datetime
import shutil
import settings
import sys

class Archivist:
	def __init__(self):
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
			df.to_csv(path+"/"+filename+".csv")
		except:
			print("P-Predictor: Could not save the file: "+path+"/"+filename+".csv")

	def save(self, filename, string):
		file = open(filename, "a")

		file.write(string + "\n")

		file.close()

	def saveEvo(self, string):
		path = "./output/evolution/evo_log_{}.log".format(settings.output_file_name)
		file = open(path, "a")
		file.write("{}\n".format(string))
		file.close()

	def saveModel(self, df):
		path = "./output/model/model_{}.csv".format(settings.output_file_name)
		df.to_csv(path)


	def setup(self, archive=False):
		# Ensure the appropriate directories exist. Remove any existing log file and archive them using appropriate time/date
		outputDir = "./output"
		evoLogDir = "./output/evolution"
		modelDir = "./output/model"
		archiveDir = "./archive"
		subDir = "./subs"
		batchDir = "./batch"

		# Check if everything is in place
		try:
			if not os.path.exists(outputDir):
				os.makedirs(outputDir)
			if not os.path.exists(modelDir):
				os.makedirs(modelDir)
			if not os.path.exists(evoLogDir):
				os.makedirs(evoLogDir)
			if not os.path.exists(archiveDir):
				os.makedirs(archiveDir)
			if not os.path.exists(subDir):
				os.makedirs(subDir)
			if not os.path.exists(batchDir):
				os.makedirs(batchDir)
		except:
			print(sys.exc_info())
			exit(1)

		# if archive, run this, else, clear the existing files and return.
		if not archive:
			# reset the output file
			open("{}/{}{}.log".format(evoLogDir, "evo_log_", settings.output_file_name), 'w').close()
			return

		evoLogs = os.listdir(evoLogDir)
		models = os.listdir(modelDir)
		
		archivePath = "{}/{}".format(archiveDir, datetime.date.today())

		try:
			if not os.path.exists(archivePath):
				os.makedirs(archivePath)
				os.makedirs("{}/{}".format(archivePath, "model"))
				os.makedirs("{}/{}".format(archivePath, "evolution"))
		except:
			raise("Couldn't create the archive folder")

		for file in models:
			path = "{}/{}".format(modelDir, file)
			pathTo = "{}/{}/{}".format(archivePath, "model", file)
			# Make sure the file already doesn't exists. 
			# ToDo:: Make this overwrite protection a bit smarter maybe? So it doesn't add the numbers to the end.
			tempPath = pathTo
			counter = 1
			while os.path.exists(tempPath):
				tempPath = pathTo + repr(counter)
				counter += 1
			pathTo = tempPath
			shutil.move(path, pathTo)

		for file in evoLogs:
			path = "{}/{}".format(evoLogDir, file)
			pathTo = "{}/{}/{}".format(archivePath, "evolution", file)
			# Make sure the file already doesn't exists. 
			# ToDo:: Make this overwrite protection a bit smarter maybe? So it doesn't add the numbers to the end.
			tempPath = pathTo
			counter = 1
			while os.path.exists(tempPath):
				tempPath = pathTo + repr(counter)
				counter += 1
			pathTo = tempPath
			shutil.move(path, pathTo)

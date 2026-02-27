import os
import pandas as pd
import eletility


class Archivist:
    def __init__(self, config):
        self.config = config

    def saveEvo(self, string):
        with open(self.config["output_evo"], "a") as f:
            f.write("{}\n".format(string))

    def saveModel(self, df):
        path = self.config["output_model"]
        df.to_csv(path)

    def setup(self, archive=False):
        os.makedirs("./output", exist_ok=True)
        os.makedirs("./archive", exist_ok=True)
        eletility.Files().writeTruncate(self.config["output_evo"], "")

import os
import pandas as pd
import eletility


class Archivist:
    def __init__(self, config):
        self.config = config

    @staticmethod
    def _ensure_parent_dir(path):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    def saveEvo(self, string):
        self._ensure_parent_dir(self.config["output_evo"])
        with open(self.config["output_evo"], "a") as f:
            f.write("{}\n".format(string))

    def saveModel(self, df):
        path = self.config["output_model"]
        self._ensure_parent_dir(path)
        df.to_csv(path, index=False)

    def setup(self, archive=False):
        os.makedirs("./output", exist_ok=True)
        os.makedirs("./archive", exist_ok=True)
        self._ensure_parent_dir(self.config["output_evo"])
        self._ensure_parent_dir(self.config["output_model"])
        eletility.Files().writeTruncate(self.config["output_evo"], "")

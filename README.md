# POET Software Installation and Usage Guide

POET is a tool used for replacing directed evolution to find fitter proteins with respect to a protein function.

## Prerequisites

- Ensure you have Python installed on your system.
- It's recommended to have PyCharm installed for easier setup and execution. You can use the community edition or obtain a student license for free.

## Installation

### Setting up a Virtual Environment (venv)

A `venv` or virtual environment is a Python tool that allows for the installation of libraries or modules in an isolated environment instead of system-wide. This helps in managing dependencies specific to different projects.

If you're using PyCharm:
1. Open the POET folder as a Python project: `File -> Open`.
2. If prompted to select an interpreter (which sets up a venv for you), locate "<No interpreter>" on the bottom right of PyCharm and create/add a new interpreter.
3. If the "Base interpreter" does not locate your Python installation, you have to manually tell it where Python is installed.

### Installing Required Libraries

The required libraries for POET are `pandas` and `scipy`.

#### Using PyCharm (Recommended):

1. Click on the Python interpreter button which should now be on the bottom right of the PyCharm app and select Interpreter Settings.
2. Here, click on the `+` sign, search for `pandas` and `scipy`, and click on the "install package" button for each.

#### Using Terminal:

1. Run the following commands:
   ```bash
   pip install pandas
   pip install scipy
   ```

**Note:** In some cases, installing via pip might require troubleshooting or additional steps. Make sure to check for meaningful error messages that might guide the resolution.

#### Running POET through PyCharm:

1- Clone the project code to your desired directory. This can by directly downloading the project files from the repository, cloning the project through Github desktop application or by cloning the project using a terminal like the following:

```
git clone https://github.com/elemenohpi/POET.git
``` 

1- Navigate through the project files using the internal file browser (usually on the left side) and find poet.py.

2- Right-click on poet.py and select Run 'poet'.

3- For run-time options, locate 'poet' on the top right section of PyCharm next to the play/run button. Click on it and select Edit configuration. In the Parameters textbox, add the run-time options (e.g., -predict).


#### Running POET through Terminal:

You can run POET through a terminal using the following command in the main directory of the project:

```
python run.py
```

additional options and runtime commands can be found by adding `-h` option to the command above to access the help instructions:

```
python run.py -h
```

If there are further issues related to Python not being recognized, ensure Python is added to your system's PATH. [Here's a guide on how to do that.](https://realpython.com/add-python-to-path/)

Make sure all the datasets (`available in data/`) and configuration files (`config.ini`) are updated according to your goal before running POET experiments. 
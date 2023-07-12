import argparse
import os
import strutils
import strformat

# poet.py -h
# 
# -config CONFIG configures the application
# -o O output filename, but usually determined in config file itself
# -r R number of GP iterations, but usually determined in config file itself
# -predict
#   count [int]
#   sequence size [int]
#   generations [int]
#   model [file_path]
#   config [file_path]
# -f F [file_path] returns fitness of model according to a dataset
# -md MD [file_path] returns table of sequence, actual_value, predicted_value
# -c [dir_path] evaluates all models according to dataset, shows all, but also shows best
# -learn LEARN [file_path] path to the dataset
# -hpcc Runs on the HPCC, asks interactive questions and overrides config
#   exper title [string]
#   walltime [int] (hr)
#   repeats [int] (replicates or jobs)
#   generations [int]
#   random seed [int]
#   configuration file [path]:
#   confirmation (y/n)
# 
#   dir created: May-03-2023-{expr-title}/
#     config.ini errors/ logs/ models/ slurms/ subs/

var parser = newParser:
  option("-config", help="Takes the config file to configure the application")
  option("-o", help="Output file name. Include extension")
  option("-r", help="Number of GP iterations to find a model")
  option("-seed", help="The random seed")
  option("-mo", help="Specifies a path to the model output. Include extension")
  flag("-predict", help="Number of potential protein sequences you want the program to predict")
  # count [int]
  # sequence size [int]
  # generations [int]
  # model [file_path]
  # config [file_path]
  option("-f", help="Gets path to a model as it's input and returns the fitness of it")
  option("-md", help="Expects a model to be given. Returns a table of predictions, actual values and RMSE for that model")
  option("-c", help="evaluates all models according to dataset, shows all, but also shows best")
  flag("-hpcc", help="Runs on the HPCC, interactive prompting and overrides config")
  option("-learn", help="Path to the learn data (format: csv)", default=some("data/learnall.csv"))
  # exper title [string]
  # walltime [int] (hr)
  # repeats [int] (replicates or jobs)
  # generations [int]
  # random seed [int]
  # configuration file [path]:
  # confirmation (y/n)

type RunMode = enum
  ModelFitness
  MeasureDatasetAgainstModels
  CompareModels
  HPCC
  Predict
type FitnessAlg = enum
  RMSE
  Correlation

proc parse[T](data: string): T =
    try:
      when T is int:
        result = parseInt data
      elif T is float:
        result = parseFloat data
      elif T is bool:
        result = parseBool data
      elif T is string:
        result = strip data
      elif T is RunMode:
        result = parseEnum[T](data)
    except ValueError:
      echo "Bad data: '", data, "'"

proc userget[T](msg: string): T =
  var input: string
  var input_needed = true
  while input_needed:
    stdout.write msg
    try:
      input = stdin.readLine
      result = parse[T](input)
      input_needed = false
    except ValueError:
      echo "Bad input: '", input, "'"

type PredictParams = object
  count, sequenceSize, iterations: int
  model, config: string
type HPCCParams = object
  title, config: string
  hours, reps: int
  confirmation: bool
type Config = object
  runmode: RunMode
  path: string = ""
  # general
  seed = 110
  # translation
  translationTable = "amino_to_amino.csv"
  diversitySelection = true
  diversityWeight = 0.01
  # Hyper-parameters
  populationSize = 100
  maximumRuleSize = 6
  maximumRuleCount = 100
  ruleWeightMin = 0.0
  ruleWeightMax = 10.0
  generations = 10005
  crossoverUnusedSelectionChance = 0.2
  fitnessAlg: FitnessAlg
  tournamentSize = 5
  # output
  popLogInterval = -1
  outputEvo = "output/evo.csv"
  outputModel = "data/unseen.csv"
  learnData = "data/learnall.csv"
  # Mutation - on model
  mutAddRule = 0.2
  mutRemoveRule = 0.2
  # Mutation - on rule
  mutChangeWeight = 0.2
  mutAddToPattern = 0.1
  mutRemoveFromPattern = 0.1
  # Predictor
  predictionModel = "output/mode/model19.csv"

import parsecfg
import std/[strutils, streams]
import tables

proc readConfig(name: string): Config =
  # 1) first load k,v strings from file
  # 2) convert to typed Config object
  if not fileExists name: return
  var settings: Table[string, string]
  var file = newFileStream(name, fmRead)
  assert file != nil, "cannot open " & name
  var cfgParser: CfgParser
  open(cfgParser, file, name)
  while true:
    var e = next(cfgParser)
    case e.kind
    # ignore these first 3
    of cfgEof: break
    of cfgSectionStart: break
    of cfgOption: break
    # key = value
    of cfgKeyValuePair:
      settings[e.key.nimIdentNormalize] = e.value
    of cfgError:
      echo e.msg
  close(cfgParser)
  # 2) convert to typed Config object (result)
  for k,v in fieldPairs result:
    if settings.hasKey k.nimIdentNormalize:
      v = parse[typeof v](settings[k.nimIdentNormalize])

proc writeConfig(config: Config) =
  var contents = ""
  for k,v in fieldPairs config:
    contents &= "$1 = $2\n".format(k,v)
  writeFile("config.ini", contents)

# try parsing command line params
try:
  var opts = parser.parse(commandLineParams())
except ShortCircuit as err:
  if err.flag == "argparse_help":
    echo err.help
    quit(1)
except UsageError:
  stderr.writeLine getCurrentExceptionMsg()
  quit(1)

var config = readConfig "config.ini"
var opts = parser.parse(commandLineParams())

if opts.seed_opt.isSome: config.seed = parse[int] opts.seed
if opts.r_opt.isSome: config.generations = parse[int] opts.r
if opts.o_opt.isSome: config.outputEvo = opts.o
if opts.mo_opt.isSome: config.outputModel = opts.mo
if opts.learn_opt.isSome: config.learnData = opts.learn

if opts.f_opt.isSome:
  config.runMode = RunMode.ModelFitness
  config.path = opts.f
if opts.c_opt.isSome:
  config.runMode = RunMode.CompareModels
  config.path = opts.c
if opts.hpcc:
  config.runMode = RunMode.HPCC
if opts.predict:
  config.runMode = RunMode.Predict
if opts.md_opt.isSome:
  config.runMode = RunMode.MeasureDatasetAgainstModels
  config.path = opts.md

echo config


# Write config file
#var config = Config()
#writeConfig config

#var config = readConfig "config.ini"
#echo config

# -predict
#   count [int]
#   sequence size [int]
#   generations [int]
#   model [file_path]
#   config [file_path]

#echo """Predicting proteins:
#Population pool is set to be 10000
#================================"""
#p.predict.count = userget[int] "Enter prediction count: "
#p.predict.sequenceSize = userget[int] "Enter protein sequences size: "
#p.generations = userget[int] "Enter number of evolutionary iterations (Larger values results in more confident and " &
#                             "yet similar predictions.\n Lower values makes room for novelty but the prediction might" &
#                             " not be as accurate): "
#p.predict.model = userget[string] "Enter the path to the predictor model(s) directory: "
#p.predict.config = userget[string] "Enter the path to the config file used to run the experiment: "

# do HPCC
#echo """
#POET is running on the HPCC experiment mode. Please enter the following information
#======================================"""
#p.hpcc.title  = userget[string] "Experiment Title: "
#p.hpcc.hours  = userget[int] "Requested Hours: "
#p.hpcc.reps = userget[int] "Number of Experiment Repeats: "
#p.generations = userget[int] "Evolutionary Generations: "
#p.seed        = userget[int] "Starting Random Seed: "
#p.hpcc.config = userget[string] "Experiment configuration file: "
#echo """
##======================================
##Experiment Title: $#
##Hours: $#
##Repeats: $#
##Generations: $#
##Seed: $#
##Config: $#""".format(p.hpcc.title, p.hpcc.hours, p.hpcc.reps, p.generations, p.seed, p.hpcc.config)
#echo "(press ENTER to continue if okay, CTRL-C otherwise)"
#discard stdin.readLine

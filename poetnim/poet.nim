import tables, sequtils, algorithm, strformat, strutils, sets
import hashes
from random import rand, randomize, sample
from math import round, copySign, sqrt, sum
from sugar import collect
from algorithm import sort, SortOrder
from std/decls import byaddr
import times # for benchmarking
import os

## community libraries
import seqmath # math ops defined for lists
import manu # matrix math
import flatty # serialization


const
  GENS = 2_000
  POPULATION_SIZE = 100
  RULE_WEIGHT_MIN = 0.0
  RULE_WEIGHT_MAX = 10.0
  RULE_SIZE_MAX = 6
  RULE_COUNT_MAX = 100

## non-const so they can be
## modified during optional simulated annealing
var
  MUT_MODEL_ADD_RULE = 0.2
  MUT_MODEL_REMOVE_RULE = 0.2
  MUT_RULE_CHANGE_WEIGHT = 0.2
  MUT_RULE_ADD_AA = 0.1
  MUT_RULE_REMOVE_AA = 0.1
  MUT_RULE_CHANGE_AA = 0.1
  UNUSED_SELECTION_CHANCE = 0.2


## optional compile flags

## general debug mode
## provides extra checks, at runtime cost
when defined(debug):
  echo "-d:debug (on)"
else:
  echo "-d:debug (off)"

## simulated annealing
## for a variety of properties
when defined(annealing):
  echo "-d:annealing (on)"
else:
  echo "-d:annealing (off)"

## deduplication of rules
## removes duplicate patterns
when defined(dedup):
  echo "-d:dedup (on)"
else:
  echo "-d:dedup (off)"

## enables the intended form
## of model shrinking
## 1) delete short-to-long unused rules
## 2) delete short-to-long used rules
when defined(shrinkend):
  echo "-d:shrinkend (on)"
else:
  echo "-d:shrinkend (off)"


## enables evaluating all rules, regardless
## if they have been found in the target sequence
## already at that position
when defined(evalallrules):
  echo "-d:evalallrules (on)"
  const EVAL_ALL_RULES = true
else:
  echo "-d:evalallrules (off)"
  const EVAL_ALL_RULES = false


iterator mitemsReverse*[T](a: var openArray[T]): var T {.inline.} =
  ## like the build-in mitems, but in reverse
  ## iterates over each item of `a` so that you can modify the yielded value.
  var i = a.high
  while i >= 0:
    yield a[i]
    dec i


iterator mpairsReverse*[T](a: var openArray[T]): tuple[key: int, val: var T] {.inline.} =
  ## like the build-in mpairs, but in reverse
  ## useful if modifying/deleting elements in the list
  ## Iterates over each item of `a`. Yields `(index, a[index])` pairs.
  ## `a[index]` can be modified.
  var i = a.high
  while i >= 0:
    yield (i, a[i])
    dec i


## simulated annealing
template lerp(a,b,weight:float):float =
  if a<b: (b-a)*weight+a
  else:   a-(a-b)*weight
template interpolateProps(gen,maxgen: int) =
  when defined(annealing):
    let weight = gen.float / maxgen.float
    MUT_MODEL_ADD_RULE = lerp(0.4, 0.2, weight)
    MUT_MODEL_REMOVE_RULE = lerp(0.4, 0.2, weight)
    MUT_RULE_CHANGE_WEIGHT = lerp(0.4, 0.2, weight)
    MUT_RULE_ADD_AA = lerp(0.3, 0.1, weight)
    MUT_RULE_REMOVE_AA = lerp(0.5, 0.1, weight)
    MUT_RULE_CHANGE_AA = lerp(0.3, 0.1, weight)

type
  Datum = object
    sequence: string
    score: float
  Data = seq[Datum]
  Rule = object
    pattern: string
    weight: float
    status: bool
  Model = object
    rules: seq[Rule]
    fitness: float


## comparisons defined for sorting to work
func `<`(a,b:Model):bool = a.fitness < b.fitness
func `>`(a,b:Model):bool = a.fitness < b.fitness

template showmem(force=false) =
  ## debug, showing memory consumption
  when defined(debug) or force:
    echo getOccupiedMem()


const AASYMBOLS = "ACDEFGHIKLMNPQRSTVWY"
proc randAASequence(length: int): string =
  ## helper to create a random
  ## amino acid sequence of length
  result.setLen length
  for i in 0..<length:
    result[i] = sample AASYMBOLS


## hydrophobicity stability map
const HSmap =
 {'A': 0.17,
  'C': -0.24,
  'D': 1.23,
  'E': 2.02,
  'F': -1.13,
  'G': 0.01,
  'H': 0.96,
  'I': -0.31,
  'K': 0.99,
  'L': -0.56,
  'M': -0.23,
  'N': 0.42,
  'P': 0.45,
  'Q': 0.58,
  'R': 0.81,
  'S': 0.13,
  'T': 0.14,
  'V': 0.07,
  'W': -1.85,
  'Y': -0.94}.toTable()
type Feature = enum
  Hydrophobicity
proc features(sequence: string, features:static[set[Feature]]): TableRef[string, float] =
  ## extendable way to calculate many different features for a protein
  ## and return the results in an annotated dictionary
  ## features example: "ACCA".features({Hydrophobicity})
  ## return example: {"hydrophobicity": 1.25,...}
  new result
  # state variables
  when Hydrophobicity in features: (var hydro_phobicity = 0.0)
  # amino acid-level analysis
  for aa in sequence:
    when Hydrophobicity in features: hydro_phobicity += HSmap[aa]
  # whole sequence-level analysis
  when Hydrophobicity in features: result[$Hydrophobicity] = hydro_phobicity


func reversed*(sv:string):string =
  ## string reversal
  var s = newString(sv.len)
  for i in 0..<sv.len: s[i] = sv[sv.high - i]
  s


type OP = enum
  lt, gt, eq
proc makeData(epochid: static[int], epochop: static[OP]): Data =
  ## loads training/testing data
  ## allows to subset the data
  ## epochid: an epoch id
  ## epochop: keep all data less than, greater than, or equal to the <epoch id>
  const original = staticRead("sources.csv")
  result = newSeq[Datum]()
  for line in original.strip.splitLines:
    let fields = line.split(',')
    let (sequence, score, epoch) = (fields[0], fields[1].parseFloat, fields[2].parseInt)
    when epochop == OP.lt:
      if epoch < epochid:
        result.add Datum(sequence: sequence, score: score)
    when epochop == OP.gt:
      if epoch > epochid:
        result.add Datum(sequence: sequence, score: score)
    when epochop == OP.eq:
      if epoch == epochid:
        result.add Datum(sequence: sequence, score: score)


proc add(mdl: var Model, rule: sink Rule) = mdl.rules.add rule
## when adding a rule to a model, use move for performance


## convenience printing of models
func `$`(mdl: Model):string =
  result = &"pattern,weight,status ({mdl.fitness.formatFloat(ffDecimal, 2)})\n"
  for rule in mdl.rules:
    result.add &"{rule.pattern},{rule.weight.formatFloat(ffDecimal, 2)},{rule.status.int}\n"

func `$`(pop: seq[Model]):string =
  for model_i,model in pop:
    if model_i > 0: result.add "\n"
    result.add &"{model.rules.len} ({model.fitness.formatFloat(ffDecimal, 2)}): {model}"


proc sort(mdl: var Model) =
  ## sorts all rules in model longest to shortest
  ## defines a closure as the comparison function
  ## depending on compile flags
  let lengthComparator = proc(x,y:Rule):int =
    when defined(dedup):
      let strlen = system.cmp(x.pattern.len, y.pattern.len)
      let strcmp = system.cmp(x.pattern, y.pattern) * -1
      if   strlen < 0: return -1
      elif strlen > 0: return 1
      elif strlen == 0 and strcmp < 0: return -1
      elif strlen == 0 and strcmp > 0: return 1
      elif strlen == 0 and strcmp == 0: return 0
    else:
      system.cmp(x.pattern.len, y.pattern.len)
  mdl.rules.sort(cmp=lengthComparator, order=SortOrder.Descending)


proc sort(models: var seq[Model]) =
  ## convenience function
  ## applies sort to all models
  # sorts all rules in all models longest to shortest
  for model in models.mitems:
    sort model

var globalCounter = 1
proc predict(mdl: var Model, target: string): float =
  ## returns score, and sets status bit of model
  ## Note: assumes mdl rules are sorted longest to shortest rule
  when defined(debug):
    ## in debug mode, do some extra safety checking for misuse
    var sortableModel = mdl
    let beforeSorting = mdl
    sort sortableModel
    let afterSorting = sortableModel
    if beforeSorting != afterSorting:
      raise (ref ValueError)(msg:"Predict: Model has unsorted rules, detected in debug mode. Predict() requires sorted models.")
  const NO_POS = -1
  var prevFound: set[0..127] # assume sequences will never be longer than 128
  var pos: int # index of found position
  var poslist = newSeqOfCap[int](12)
  var timesFound: int
  for rulei,rule in mdl.rules.mpairs:
    let pattern = rule.pattern.`$`
    let patternR = rule.pattern.`$`.reversed
    timesFound = 0
    # find forward pattern
    pos = target.find(pattern, start=0)
    while pos != NO_POS:
      #prevFound.incl pos
      #inc forwardFound
      if (pos notin prevFound) or EVAL_ALL_RULES:
        prevFound.incl pos
        inc timesFound
      pos = target.find(pattern, start=pos+1)
    # find backward battern
    pos = target.find(patternR, start=0)
    while pos != NO_POS:
      #prevFound.incl pos
      #inc backwardFound
      if (pos notin prevFound) or EVAL_ALL_RULES:
        prevFound.incl pos
        inc timesFound
      pos = target.find(patternR, start=pos+1)
    result += rule.weight * timesFound.float
    if timesFound.bool and rule.status == false:
      rule.status = true
    #mdl.rules[rulei].status = seen


proc mutate(mdl: var Model) =

  ## add rule
  if rand(1.0) <= MUT_MODEL_ADD_RULE and mdl.rules.len < RULE_COUNT_MAX:
    let length = rand 1..RULE_SIZE_MAX
    let newsequence = randAASequence(length)
    mdl.add Rule(pattern: newsequence, weight: rand(RULE_WEIGHT_MIN .. RULE_WEIGHT_MAX).round(2), status: false)

  ## remove rule
  if rand(1.0) <= MUT_MODEL_REMOVE_RULE and bool(mdl.rules.len):
    mdl.rules.delete rand(mdl.rules.high)

  ## over all rules
  for rule_i,rule in mdl.rules.mpairsReverse:
    let pattern {.byaddr.} = rule.pattern

    ## change weight
    if rand(1.0) <= MUT_RULE_CHANGE_WEIGHT:
      rule.weight += rand(-1.0..1.0)

    ## add amino acid
    if rand(1.0) <= MUT_RULE_ADD_AA and pattern.len < RULE_SIZE_MAX:
      pattern.insert($sample(AASYMBOLS), rand pattern.len)

    ## remove amino acid
    if rand(1.0) <= MUT_RULE_REMOVE_AA and bool(pattern.len):
      let pos = rand pattern.high
      pattern.delete(pos..pos)
      if pattern.len == 0:
        mdl.rules.delete rule_i
        continue

    ## point-mutation amino acid (not in original py-POET)
    #if rand(1.0) < MUT_RULE_CHANGE_AA:
    #  let pos = rand thisPattern.high
    #  thisPattern[pos] = sample AASYMBOLS


proc mutate(models: var seq[Model]) =
  ## convenience function
  ## allowing mutating each model
  for model in models.mitems:
    mutate model


## element-wise ops
func `-`[T](x,y:seq[T]):seq[T] = eSub(x,y)
func `-`[T](x:seq[T],y:T):seq[T] = eSub(x,y)
func `+`[T](x,y:seq[T]):seq[T] = eAdd(x,y)
func `+`[T](x:seq[T],y:T):seq[T] = eAdd(x,y)
func `*`[T](x,y:seq[T]):seq[T] = eMul(x,y)
func `*`[T](x:seq[T],y:T):seq[T] = eMul(x,y)
func `/`[T](x,y:seq[T]):seq[T] = eDiv(x,y)
func `/`[T](x:seq[T],y:T):seq[T] = eDiv(x,y)

#####################
## stats functions ##
#####################

proc mean[T](x:openArray[T]): float =
  ## same as in std/stats
  x.sum.float / x.len.float


proc linearRegression[T:SomeFloat](x,y:seq[T]):seq[T] =
  ## returns lowest power first (the y-intercept),
  ## which is opposite of how python does it
  var X = matrix[T](x.len, 2)
  let Y = matrix[T](y, y.len)
  X[0..<X.m, 1..1] = matrix[T](x, x.len)
  X[0..<X.m, 0..0] = ones[T](X.m, 1)
  let r = (X.transpose * X).inverse * (X.transpose * Y)
  result = r.getColumnPacked


proc rmse[T](x,y: seq[T]): float =
  ## standard rmse
  let n = len(x)
  when defined(debug):
    if n != len(y):
      raise newException(ValueError,"x and y must have the same length.")
    if n < 2:
      raise newException(ValueError,"x and y must have length at least 2.")
  for i in 0..<n:
    result += (y[i]-x[i])*(y[i]-x[i])
  result = sqrt(result / n.float)


proc pearsonr[T](x,y: seq[T]): float =
  ## standard pearsonr
  # ported from SciPy: https://github.com/scipy/scipy/blob/v1.9.3/scipy/stats/_stats_py.py#L4225-L4492
  let n = len(x)
  when defined(debug):
    if n != len(y):
      raise newException(ValueError,"x and y must have the same length.")
    if n < 2:
      raise newException(ValueError,"x and y must have length at least 2.")
  when defined(debug):
    if x.allIt(it == x[0]) or y.allIt(it == x[0]):
      echo "Warnings: An input array is constant; the correlation coefficient is not defined. returning NaN."
      return NaN
  if n == 2:
    let r = copySign(1.0.T, x[1] - x[0])*copySign(1.0.T, y[1] - y[0])
    return r
  let
    xm = x - x.mean
    ym = y - y.mean
    normxm = sqrt(sum(xm*xm))
    normym = sqrt(sum(ym*ym))
  when defined(debug):
    const threshold = 1e-13.T
    if normxm < threshold*abs(x.mean) or normym < threshold*abs(y.mean):
      echo "Warning: An input array is nearly constant; the computed correlation coefficient may be inaccurate."
  let r = sum ((xm/normxm) * (ym/normym))
  return max(min(r, 1.0.T), -1.0.T)


var ydata:Table[int, seq[float]]
proc evaluate(mdl: var Model, data: Data) =
  ## determines fitnesses according to data
  ## caches data depending on dataframe length
  ## so epoch 1 and epoch 2 etc. are unique and cached
  var x = newSeq[float](data.len)
  if data.len notin ydata:
    ydata[data.len] = newSeq[float](data.len)
    for datum_i, datum in data:
      ydata[data.len][datum_i] = datum.score
  # reset all rule status
  for rule in mdl.rules.mitems:
    rule.status = false
  for datum_i, datum in data:
    x[datum_i] = mdl.predict datum.sequence
  let r = pearsonr(x,ydata[data.len])
  mdl.fitness = r*r


proc evaluate(models: var seq[Model], data: Data) =
  ## convenience function to evaluate all models
  for model in models.mitems:
    evaluate(model, data)


proc predict(mdl: var Model, data: Data): seq[float] =
  ## prediction with alignment (need better names btwn evaluate,predict,predict...)
  var x = newSeq[float](data.len)
  if data.len notin ydata:
    ydata[data.len] = newSeq[float](data.len)
    for datum_i, datum in data:
      ydata[data.len][datum_i] = datum.score
  for datum_i, datum in data:
    x[datum_i] = predict(mdl, datum.sequence)
  let alignment = linearRegression(x,ydata[data.len]) # [y-int, slope]
  result = x * alignment[1] + alignment[0]


proc rmse(mdl: var Model, data: Data):float =
  ## convenience function for getting rmse
  ## of a model relative to data
  var x = newSeq[float](data.len)
  if data.len notin ydata:
    ydata[data.len] = newSeq[float](data.len)
    for datum_i, datum in data:
      ydata[data.len][datum_i] = datum.score
  x = predict(mdl,data)
  let rmse = rmse(x,ydata[data.len])
  return rmse


proc makePopulation(size = POPULATION_SIZE; maxrules = RULE_COUNT_MAX div 3; rulewmin = RULE_WEIGHT_MIN; rulewmax = RULE_WEIGHT_MAX):seq[Model] =
  ## samples from source to create an initial
  ## population of models with random weights
  var models = collect:
    for _ in 1..size:
      let rules = collect:
        let numrules = rand 1..maxrules
        for _ in 1..numrules:
          let length = rand 1..RULE_SIZE_MAX
          let newsequence = randAASequence(length)
          Rule(pattern: newsequence, weight: rand(rulewmin .. rulewmax).round(2), status: false)
      Model(rules: rules)
  models


proc trim(mdl: var Model) =
  ## Note! assumes sorted rules
  when defined(dedup):
    ## remove dups
    var index = 0
    if mdl.rules.len == 0: return
    var pattern = ""
    while index < mdl.rules.len:
      if pattern != mdl.rules[index].pattern:
        pattern = mdl.rules[index].pattern
        inc index
      else:
        mdl.rules.delete index
  when defined(shrinkend):
    ## this is the conceptually intented shrinking method
    ## remove 0-status then 1-status short rules first
    var status1RulesToRemove = 0
    for i in countDown(mdl.rules.high,RULE_COUNT_MAX):
      if not mdl.rules[i].status: mdl.rules.delete i
      else: inc status1RulesToRemove
      if mdl.rules.len == RULE_COUNT_MAX: break
    if status1RulesToRemove >= mdl.rules.len:
      mdl.rules.setLen(0)
    elif status1RulesToRemove > 0:
      for i in countDown(mdl.rules.high,mdl.rules.high-status1RulesToRemove):
        mdl.rules.delete i
  else:
    ## this is the default py-POET style shrinking that
    ## begins deleting from toward the top of the list
    ## thereby preserving most length-1 patterns
    var countGreens = 0
    while len(mdl.rules) > RULE_COUNT_MAX:
      countGreens = 0
      for index in countDown(mdl.rules.high, 0):
        if (countGreens >= index):
          delete(mdl.rules,index)
          break
        else:
          if (mdl.rules[index].status == false):
            delete(mdl.rules,index)
            break
          else:
            inc countGreens

proc crossover(parentA,parentB:Model):Model =
  ## parentA,parentB: the parents
  ## result: 1 child
  let model {.byaddr.} = result
  model.rules = newSeqOfCap[Rule](RULE_COUNT_MAX) # preallocate up to max for speed
  let lenA = parentA.rules.len
  let lenB = parentB.rules.len
  let maxLen = max(lenA,lenB)
  for j in 0..<maxLen:
    if j < lenA and (parentA.rules[j].status or rand(1.0) < UNUSED_SELECTION_CHANCE):
      model.rules.add parentA.rules[j]
    if j < lenB and (parentB.rules[j].status or rand(1.0) < UNUSED_SELECTION_CHANCE):
      model.rules.add parentB.rules[j]

proc makeTournamentIndex(size,popsize:int):seq[int] =
  ## size: how many ids to generate (tournament size)
  ## popsize: 0..<popsize possible ids to sample from
  ## result: list of ids for 1 tournament
  result = collect:
    for _ in 1..size:
      rand(popsize-1)
  sort(result,order=SortOrder.Descending)

## wrappers for sorting by index
## maybe a better way to do this using std lib stuff
from std/heapqueue import HeapQueue, `[]`, pop, push, pushpop
type IDXWrapper[T] = object
  idx: int
  item: T
func `<`[T](a,b:IDXWrapper[T]):bool = a.item < b.item
func `<`[T](a:IDXWrapper[T],b:Model):bool = a.item < b.fitness
proc findBestModelsIndex(a: openArray[Model], n: Positive, sub: openArray[int]): seq[int] =
  ## a: list of models
  ## n: how many best to find
  ## sub: subset of indices to filter consideration
  ## result: list of indices of the best n models from a
  let N = min(min(n, a.len),sub.len)
  var hq: HeapQueue[IDXWrapper[type(Model.fitness)]]
  for i in 0..<N: hq.push IDXWrapper[type(Model.fitness)](idx: sub[i], item: a[sub[i]].fitness)
  var small = hq[0]
  for i in N..<sub.len:
    if a[sub[i]] > small:
      discard hq.pushpop IDXWrapper[type(Model.fitness)](idx: sub[i], item: a[sub[i]].fitness)
      small = hq[0]
  result.setLen(N)
  for i in 0..<N: result[N-1-i] = hq[i].idx

## vvv serialization - should remake without the flatty lib
proc serialize(mdl: Model): string =
  let fitness = mdl.fitness.toFlatty
  let rulesPatterns = mdl.rules.mapIt(it.pattern.toFlatty).toFlatty
  let rulesWeights = mdl.rules.mapIt(it.weight.toFlatty).toFlatty
  let rulesStati = mdl.rules.mapIt(it.status.toFlatty).toFlatty
  result = @[fitness, rulesPatterns, rulesWeights, rulesStati].toFlatty

proc deserialize(s: string, _: typedesc[Model]):Model =
  let guts = s.fromFlatty(seq[string])
  result.fitness = guts[0].fromFlatty(float)
  let rulesPatterns = guts[1].fromFlatty(seq[string]).mapIt(it.fromFlatty(string))
  let rulesWeights = guts[2].fromFlatty(seq[string]).mapIt(it.fromFlatty(float))
  let rulesStati = guts[3].fromFlatty(seq[string]).mapIt(it.fromFlatty(bool))
  for i in 0..<rulesPatterns.len:
    result.rules.add Rule(weight: rulesWeights[i], pattern: rulesPatterns[i], status: rulesStati[i])

proc save(mdl: Model, filename: string) =
  writeFile(filename=filename, content=mdl.serialize)

proc load(filename: string): Model =
  filename.readFile.deserialize(Model)
## ^^^ serialization - should remake without the flatty lib

func unique(rules: seq[Rule]): HashSet[string] =
  ## returns all unique rule patterns
  for rule in rules:
    result.incl rule.pattern

func hist(rules: seq[Rule]): seq[int] =
  ## returns histogram counts of unique rules, unsorted
  var counts: Table[string, int]
  for rule in rules:
    inc counts.mgetOrPut(rule.pattern, 0)
  result = counts.values.toSeq

proc loadPyModel(filename: string): Model =
  ## specifically reads in py-POET style saved models
  ## that look like this:
  #,pattern,weight,status
  #0,KKYNKNQN,5.120000000000001,0
  #1,GLNWFND,7.380000000000003,0
  let model {.byaddr.} = result # readability convenience
  var line: string
  let file = open filename
  defer: close file
  # skip header
  discard file.readLine(line)
  # read file
  while file.readLine(line):
    let fld = line.split(',')
    model.add Rule(pattern: fld[1], weight: fld[2].parseFloat, status: fld[3].parseInt.bool)

proc loadPyData(filename: string): Data =
  ## specifically reads in py-POET style training data
  ## that look like this:
  #sequence,fitness
  #KKKKKKKKKKKK,12.5
  #KSKSKSKSKSKS,17
  let data {.byaddr.} = result # readability convenience
  var line: string
  let file = open filename
  defer: close file
  # skip header
  discard file.readLine(line)
  # read file
  while file.readLine(line):
    let fld = line.split(',')
    data.add Datum(sequence: fld[0], score: fld[1].parseFloat)

proc main =
  echo &"population size: {POPULATION_SIZE}"
  echo &"rule count max: {RULE_COUNT_MAX}"

  ## load data
  let data = makeData(epochid = 9, epochop = OP.lt)
  let testData = makeData(epochid = 9, epochop = OP.eq)

  var elite, best: Model

  ## used to manually load a model and test it
  #block:
  #  var model = loadPyModel("models/model/model_2.csv") # example loading pyPOET formatted data
  #  let data = loadPyData("data/learn8.csv") # example loading pyPOET formatted model
  #  #var model = load("best.mdl") # example other kind of loading/saving
  #  model.evaluate(data)
  #  echo model
  #  echo &"fitness: {model.fitness}"
  #  echo &"rmse (train): {model.rmse(data)}"
  #  echo &"rmse ( test): {model.rmse(testData)}"
  #  echo ""
  #  #let predictions = model.predict(testData)
  #  #let r = pearsonr(predictions, testData.mapIt(it.score))
  #  #echo r*r
  #  quit(0)

  randomize(getCurrentProcessId()) # randomize the seed, provide int arg to set seed

  var models = makePopulation()
  var newModels = newSeq[Model](models.len)

  ## sort
  sort models

  echo &"generations: {GENS}"
  let startTime = cpuTime()

  echo "gen  fitness  max-rmse-train  max-rmse-test  max-#-rules  memory(bytes)"

  for gen in 0..<GENS:
    ## evaluate all models
    evaluate(models,data)
    ## identify best
    best = models.max
    ## make new population & seed with best
    newModels = newSeqOfCap[Model](POPULATION_SIZE) # preallocates capacity for speed
    newModels.add best
    ## tournament selection
    for _ in 2..POPULATION_SIZE:
      let group = makeTournamentIndex(size=5, popsize=models.len) # return <size> pop indices
      let best = models.findBestModelsIndex(n=2, sub=group) # return <n> pop indices
      newModels.add crossover(models[best[0]], models[best[1]])
      sort newModels[^1] # [^1] like [-1] in py
      trim newModels[^1] # only place that models are shrunk!
    ## population <- new population
    swap(models, newModels)
    ## copy elite
    elite = models[0]
    ### mutation & sort
    mutate models
    sort models
    ## elite greedy hill-climbing
    evaluate(models[0], data)
    if models[0].fitness < elite.fitness:
      models[0] = elite

    ## periodic stuff, including progress output
    if gen mod 10 == 0:
      # annealing, works if compiled with flag 'annealing'
      interpolateProps(gen, GENS)
      best = models[0]
      stdout.write gen,' ',best.fitness.formatFloat(ffDecimal,4),' ',(best.rmse(data)).formatFloat(ffDecimal,4),' ',(best.rmse(testData)).formatFloat(ffDecimal,4),' ',best.rules.len,' '
      showmem(force=true)

  let stopTime = cpuTime()
  echo &"Elapsed Time: {(stopTime-startTime).round(2)}s"

  ## update status & fitness for displaying correct information
  evaluate(models,data)

  ## num unique patterns vs total rules (number of repeats)
  ## show the number of unique rules / total # of rules
  ## The last 2 numbers at the bottom are the
  ## average number of length-1 rules, and average number of unique rules
  ## this "signature" was used to compare to python and validate the algorithm
  var stati, first, ones, uniques: seq[int]
  for model in models:
    let hist = model.rules.hist.sorted(order=SortOrder.Descending)
    #echo &"{model.rules.unique.len} / {model.rules.len}: {hist}"
    ones.add hist.count(1)
    uniques.add model.rules.unique.len
    if hist.len >= 1:
      first.add hist[0]
  stati.add models[0].rules.mapIt(it.status).count(true)
  echo &"{ones.mean.formatFloat(ffDecimal,2)} {uniques.mean.formatFloat(ffDecimal,2)} {first.mean.formatFloat(ffDecimal,2)} {stati[0]}"

  ## save all models
  for i,model in models:
    let f = open(&"model-{i}.txt",fmWrite)
    f.writeLine &"model {i}"
    f.writeLine $model
    close f

  #models.max.save("best.mdl")
  #echo "wrote 'best.mdl'"

main()

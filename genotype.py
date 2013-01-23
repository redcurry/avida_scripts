# Enable true division (1 / 2 = 0.5, not 0)
from __future__ import division

import os
import string
import random
import itertools

import avida
import stats


DATA_PATH = 'data'
TASKS = [] # Is initialized upon calling determineTasks()


### Interfaces #################################################################


# Interface for creating hybrids
class Hybridizer:
  
  def getHybrids(genotype_1, genotype_2):
    pass


# Interface for writing genotypes to a file
class GenotypeWriter:

  def write(self, path, comps):
    pass


# Interface for reading genotypes from a file
class GenotypeReader:

  def read(self, path):
    pass



### Classes ####################################################################

class TaskPerformer:

  def __init__(self):
    self.tasks = []

  # Set tasks from task string (e.g., '101100101')
  def setTasks(self, tasks):
    self.tasks = []
    for i in range(len(tasks)):
      if tasks[i] == '1':
        self.tasks.append(TASKS[i])


class Genotype(TaskPerformer):

  def __init__(self, sequence, ancestor = None):
    TaskPerformer.__init__(self)
    if isinstance(sequence, Genotype):
      self.ancestor = sequence.ancestor
      self.sequence = sequence.sequence
      self.fitness = sequence.fitness
      if ancestor != None:
        self.ancestor = ancestor
    else:
      self.ancestor = ancestor
      self.sequence = sequence
      self.fitness = -1

  def __eq__(self, other):
    if other != None:
      return self.sequence == other.sequence
    else:
      return False

  def __getitem__(self, key):
    return self.sequence[key]

  def __len__(self):
    return len(self.sequence)

  def __str__(self):
    if self.ancestor == None:
      return self.sequence
    else:
      insts = list(self.sequence)
      for i in self.mutations():
        insts[i] = self.emphasize(insts[i])
      return string.join(insts, '')

  def __hash__(self):
    return hash(self.sequence)

  def emphasize(self, c):
#    return '\033[36m' + c + '\033[0m'  # Light blue
    return '\033[31m' + c + '\033[0m'   # Red

  def relativeFitness(self):
    if self.fitness == -1:
      raise Exception('Fitness not set')
    if self.ancestor == None:
      return 1.0
    else:
      return self.fitness / self.ancestor.fitness

  def mutations(self):
    if self.ancestor == None:
      return []
    else:
      return [i for i in range(len(self.sequence)) \
        if self.sequence[i] != self.ancestor.sequence[i]]

  def mutate(self, inst, pos):
    self.sequence = \
      self.sequence[0:pos] + inst + self.sequence[(pos + 1):len(self.sequence)]


class CompensatedGenotype(Genotype):

  def __init__(self, sequence, mutant = None):
    Genotype.__init__(self, sequence)
    if isinstance(sequence, CompensatedGenotype):
      self.mutant = sequence.mutant
    else:
      self.ancestor = mutant.ancestor
      self.mutant = mutant

  def __str__(self):
    insts = list(self.sequence)
    for i in self.mutations():
      if i in self.mutant.mutations():
        insts[i] = self.mutant.emphasize(insts[i])
      else:
        insts[i] = self.emphasize(insts[i])
    return string.join(insts, '')

  def emphasize(self, c):
    return '\033[34m' + c + '\033[0m'


class HybridGenotype(Genotype):

  def __init__(self, sequence, parents = None, loci = None):
    Genotype.__init__(self, sequence)
    if isinstance(sequence, HybridGenotype):
      self.parents = list(sequence.parents)
      self.loci = list(sequence.loci)
    else:
      if parents[0].ancestor != None:
        self.ancestor = parents[0].ancestor
      else:
        self.ancestor = parents[1].ancestor
      self.parents = parents
      self.loci = loci

  def __str__(self):
    if self.loci == None:
      return Genotype.__str__(self)
    insts = list(self.sequence)
    for i in range(len(insts)):
      if i in self.loci:
        insts[i] = self.emphasize(insts[i])
    return string.join(insts, '')

  def emphasize(self, c):
    return '\033[33m' + c + '\033[0m'


class UnlinkedHybridGenotype(Genotype):

  def __init__(self, sequence, parent_1, parent_2):
    Genotype.__init__(self, sequence)
    if isinstance(sequence, UnlinkedHybridGenotype):
      self.mutations_1 = list(sequence.mutations_1)
      self.mutations_2 = list(sequence.mutations_2)
    else:
      self.ancestor = parent_1.ancestor
      self.mutations_1 = \
        [pos for pos in self.mutations() if self[pos] == parent_1[pos]]
      self.mutations_2 = \
        [pos for pos in self.mutations() if pos not in self.mutations_1]

  def __str__(self):
    return ''.join(['\033[34m' + self[i] + '\033[0m' if i in self.mutations_1 \
               else '\033[31m' + self[i] + '\033[0m' if i in self.mutations_2 \
               else self[i] for i in range(len(self))])


class OneModuleHybridizer(Hybridizer):

  def getHybrids(self, genotype_1, genotype_2):

    seq_1 = genotype_1.sequence
    seq_2 = genotype_2.sequence

    hybrids = []
    for startPos in range(0, len(seq_1)):
      for endPos in range(startPos + 1, len(seq_1)):
        hybridSeqs = self.getCrossedHybrids(seq_1, seq_2, startPos, endPos)
        hybrids.append(HybridGenotype(hybridSeqs[0], [genotype_1, genotype_2], \
          range(startPos, endPos)))
        hybrids.append(HybridGenotype(hybridSeqs[1], [genotype_1, genotype_2], \
          range(0, startPos) + range(endPos, len(seq_1))))

    return hybrids


  def getRandomHybrids(self, genotype_1, genotype_2, n):

    seq_1 = genotype_1.sequence
    seq_2 = genotype_2.sequence

    hybrids = []
    for m in range(n):
      startPos = random.randrange(0, len(seq_1))
      if startPos == len(seq_1) - 1:  # Handle special case
        endPos = startPos
      else:
        endPos = random.randrange(startPos + 1, len(seq_1))
      hybridSeqs = self.getCrossedHybrids(seq_1, seq_2, startPos, endPos)
      hybrid_1 = HybridGenotype(hybridSeqs[0], [genotype_1, genotype_2], \
        range(startPos, endPos))
      hybrid_2 = HybridGenotype(hybridSeqs[1], [genotype_1, genotype_2], \
        range(0, startPos) + range(endPos, len(seq_1)))
      hybrids.append(random.choice([hybrid_1, hybrid_2]))

    return hybrids


  def getCrossedHybrids(self, seq_1, seq_2, startPos, endPos):

    hybrid_1 = seq_1[0:startPos] + seq_2[startPos:endPos] + seq_1[endPos:]
    hybrid_2 = seq_2[0:startPos] + seq_1[startPos:endPos] + seq_2[endPos:]

    return [hybrid_1, hybrid_2]


  def getMeanHybridRelativeFitness(self, genotype_1, genotype_2):

    hybrids = self.getHybrids(genotype_1, genotype_2)
    measureFitness(hybrids)
    return stats.mean([h.relativeFitness() for h in hybrids])


  def getMeanExpectedHybridRelativeFitness(self, genotype_1, genotype_2):

    ancestor_1 = genotype_1.ancestor
    ancestor_2 = genotype_2.ancestor

    if ancestor_1 == None:
      ancestor_1 = genotype_1

    if ancestor_2 == None:
      ancestor_2 = genotype_2

    hybrids_1 = self.getHybrids(genotype_1, ancestor_1)
    hybrids_2 = self.getHybrids(ancestor_2, genotype_2)

    measureFitness(hybrids_1 + hybrids_2)

    expectedHybridFitness = [h1.relativeFitness() * h2.relativeFitness() \
      for h1, h2 in zip(hybrids_1, hybrids_2)]

    return stats.mean(expectedHybridFitness)


class UnlinkedHybridizer:

  def getRandomHybrids(self, genotype_1, genotype_2, n):
    return [self.getRandomHybrid(genotype_1, genotype_2) for i in range(n)]

  def getRandomHybrid(self, genotype_1, genotype_2):

    seq_1 = genotype_1.sequence
    seq_2 = genotype_2.sequence

    hybrid_seq_list = \
      [c_1 if random.random() < 0.5 else c_2 for c_1, c_2 in zip(seq_1, seq_2)]
    hybrid_seq = ''.join(hybrid_seq_list)

    return UnlinkedHybridGenotype(hybrid_seq, genotype_1, genotype_2)


# Recombines by treating mutations (not every instruction) as unlinked
# (i.e., returns hybrids with every combination of mutations), but has the
# option to do random combinations if the number of mutations is high;
# a limit of 0 does every combination; otherwise, it does n random combinations
class UnlinkedMutationsHybridizer:

  def __init__(self, n, genotypes_to_skip = []):

    self.n = n
    self.genotypes_to_skip = genotypes_to_skip

  def getHybrids(self, genotype_1, genotype_2, include_1 = [], include_2 = []):

    include_1 = set(include_1)
    include_2 = set(include_2)
    include = include_1 | include_2

    # (Genotypes should share their ancestor)
    ancestor = genotype_1.ancestor

    mutations_1 = set(genotype_1.mutations())
    mutations_2 = set(genotype_2.mutations())

    # No mutations to recombine, return a parent (parents should be identical),
    # unless hybrids that look like parents should not be produced
    if len(mutations_1) == 0 and len(mutations_2) == 0:
      parent = [Genotype(genotype_1, ancestor)]
      if parent in self.genotypes_to_skip:
        return []
      else:
        return parent

    seq_1 = genotype_1.sequence
    seq_2 = genotype_2.sequence

    if len(mutations_1 & mutations_2) > 0:
      print 'Some mutations overlap between parents'

    # Combine mutations (set union)
    mutations = mutations_1 | mutations_2

    hybrids = []

    samples_1 = []
    samples_2 = []

    for i in range(self.n):

      sample_1 = self.sample_mutations(mutations_1)
      sample_2 = self.sample_mutations(mutations_2)

      # Ensure at least one sample contains a mutation
      while len(sample_1) == 0 and len(sample_2) == 0:
        sample_1 = self.sample_mutations(mutations_1)
        sample_2 = self.sample_mutations(mutations_2)

      # Ensure that sample_1 or sample_2 contains at least one
      # mutation from their respective include_1 or include_2
      if len(include_1) > 0:
        if len(include_2) > 0:
          while len(include_1 & sample_1) == 0 and \
                len(include_2 & sample_2) == 0:
            sample_1 = self.sample_mutations(mutations_1)
            sample_2 = self.sample_mutations(mutations_2)
        else:
          while len(include_1 & sample_1) == 0:
            sample_1 = self.sample_mutations(mutations_1)
      else:
        if len(include_2) > 0:
          while len(include_2 & sample_2) == 0:
            sample_2 = self.sample_mutations(mutations_2)

      samples_1.append(sample_1)
      samples_2.append(sample_2)

    # Create hybrids with each combination of mutations,
    # ensuring that mutational overlaps are resolved by giving
    # priority to the parent whose mutation is in an 'include'
    for sample_1, sample_2 in zip(samples_1, samples_2):
      hybrid = Genotype(ancestor, ancestor)

      combination = sample_1 | sample_2
      for mutation in combination:
        if mutation in sample_1:
          if mutation in sample_2:
            if mutation in include_1:
              if mutation in include_2:
                # If the mutation is in both samples and in both includes,
                # then pick a random parent to use
                hybrid.mutate(random.choice([seq_1, seq_2])[mutation], \
                  mutation)
              else:
                # Mutation is in both samples, but only in include_1
                hybrid.mutate(seq_1[mutation], mutation)
            else:
              if mutation in include_2:
                # Mutation is both samples, but only in include_2
                hybrid.mutate(seq_2[mutation], mutation)
              else:
                # Mutation is in both samples, but in neither include,
                # then pick a random parent to use
                hybrid.mutate(random.choice([seq_1, seq_2])[mutation], \
                  mutation)
          else:
            # Mutation is only in sample_1
            hybrid.mutate(seq_1[mutation], mutation)
        else:
          # Mutation is only in sample_2
          hybrid.mutate(seq_2[mutation], mutation)

      # Skip certain hybrids
      if hybrid in self.genotypes_to_skip:
        continue

      hybrids.append(hybrid)

    return hybrids

  def sample_mutations(self, mutations):
    if len(mutations) > 0:
      return set(random.sample(mutations, random.randint(0, len(mutations))))
    else:
      return set()


class SharedNeutralReader(GenotypeReader):

  def read(self, path):

    # Add neutral genotypes to this list
    neutrals = []

    # Open the file for reading
    file = open(path, 'r')

    # Get number of ancestors
    ancestorCount = int(file.readline())

    for ancestorNum in range(ancestorCount):

      # Read ancestor as a genotype
      ancestor = Genotype(file.readline().strip())

      # Get number of neutral genotypes for current ancestor
      neutralCount = int(file.readline())
        
      for neutralNum in range(neutralCount):

        # Read compensated genotype
        neutral = Genotype(file.readline().strip(), ancestor)
        neutrals.append(neutral)

    file.close()

    return neutrals


### Functions ##################################################################


def measureFitness(genotypes, environmentFileName = None):

  # If genotypes is a single genotype, convert to a list
  if not isinstance(genotypes, list):
    genotypes = [genotypes]

  # Copy the sequences of each genotype into a new sequence string list
  # because avida.calc_info takes in a list of sequence strings
  seqs = [genotype.sequence for genotype in genotypes]

  # Get fitnesses and viability
  avida.set_environment(environmentFileName)
  seqs_info = avida.get_info(seqs, ['fitness', 'viable'])
  fitnesses = seqs_info['fitness']
  viability = seqs_info['viable']
  avida.reset_environment()

  # Update each genotype, ensuring unviable sequences have 0 fitness
  for genotype, fitness, viable in zip(genotypes, fitnesses, viability):
    genotype.fitness = fitness * viable


# Returns dictionary {update: genotype} for the specified run
# Note: Does not calculate fitnesses
def readRunGenotypes(runPath, auto_ancestor = False, ancestor = None):

  # Read the run file to get update and sequence lists
  runData = avida.read_data(runPath, ['update', 'dom_sequence'])
  updates = runData['update']
  sequences = runData['dom_sequence']

  # Convert data to dictionary {update: genotype}
  genotypes = {}
  for update, sequence in zip(updates, sequences):
    genotypes[update] = Genotype(sequence, ancestor)

  # Set the first genotype (min update) as the ancestor
  if auto_ancestor:
    ancestor = genotypes[min(updates)]
    for update in updates:
      genotypes[update].ancestor = ancestor

  return genotypes


# Same as readRunGenotypes(), but reads from a simple file
# without a header in the format: update sequence
def readPlainRunGenotypes(run_path, auto_ancestor = False):

  run = {}
  for line in open(run_path):
    line_parts = line.split();
    update = int(line_parts[0])
    sequence = line_parts[1]
    run[update] = Genotype(sequence)

  if auto_ancestor:
    ancestor = run[min(run)]
    for update in run:
      run[update].ancestor = ancestor

  return run


# Returns dictionary {update: mean_fitness} for the specified run
def readRunMeanFitness(runPath):

  # Read the run file to get update and mean fitness lists
  runData = avida.read_data(runPath, ['update', 'ave_fitness'])
  updates = runData['update']
  fitnesses = runData['ave_fitness']

  # Convert data to dictionary {update: fitness}
  fitnessData = {}
  for update, fitness in zip(updates, fitnesses):
    fitnessData[update] = fitness

  return fitnessData


# Returns list of all genotypes for the specified population file
# Note: Does not calculate fitness
# Note: If z == True, assumes compressed file: 'detail-10000.pop.gz'
def readPopulation(runPath, z = False):

  # Read population file
  popData = avida.read_data(runPath, ['num_cpus', 'sequence'], z)
  numbers = popData['num_cpus']
  sequences = popData['sequence']

  # Create genotypes based on the sequences and their number
  genotypes = []
  for number, sequence in zip(numbers, sequences):
    for n in range(int(number)):
      genotypes.append(Genotype(sequence))

  return genotypes


# Returns the consensus genotype of a list of genotypes
# Note: Does not calculate its fitness
def getConsensusGenotype(genotypes):

  genome_length = len(genotypes[0].sequence)
  consensus_seq_list = [i for i in range(genome_length)]
  for i in range(genome_length):
    inst_counts = {}
    for genotype in genotypes:
      inst = genotype.sequence[i]
      if inst in inst_counts:
        inst_counts[inst] += 1
      else:
        inst_counts[inst] = 1

    max_inst = inst_counts.keys()[0]
    for inst in inst_counts:
      if inst_counts[inst] > inst_counts[max_inst]:
        max_inst = inst
    consensus_seq_list[i] = max_inst

  return Genotype(''.join(consensus_seq_list))


# Returns the consensus genotypes of a run as {update: genotype}
# Note: Does not calculate fitness
# Note: Population file must be of the form, e.g., 'detail-10000.pop'
# Note: If z == True, assumes compressed file: 'detail-10000.pop.gz'
def getRunConsensusGenotypes \
  (runDir, z = False, ancestor = None, start = -1, end = -1):

  run = {}

  runDataDir = os.path.join(runDir, DATA_PATH)
  for popFileName in os.listdir(runDataDir):

    # Must be a population file
    if not (popFileName.endswith('.pop') or popFileName.endswith('.pop.gz')):
      continue

    # Must have a '-' in the file name
    if popFileName.find('-') < 0:
      continue

    # Get the update from the file name
    update = int(popFileName.replace('-', '.').split('.')[1])

    if start != -1:
      if update < start:
        continue

    if end != -1:
      if update > end:
        continue

    # Get the consensus genotype from the population file
    popPath = os.path.join(runDataDir, popFileName)
    cg = getConsensusGenotype(readPopulation(popPath, z))
    cg.ancestor = ancestor

    # Update dictionary
    run[update] = cg

  return run


def getRunRandomGenotypes(runDir):

  global DATA_PATH

  run = {}

  runDataDir = os.path.join(runDir, DATA_PATH)
  for popFileName in os.listdir(runDataDir):

    if not popFileName.endswith('.pop'):
      continue
    if popFileName.find('-') < 0:
      continue

    update = int(popFileName.replace('-', '.').split('.')[1])

    popPath = os.path.join(runDataDir, popFileName)
    pop = readPopulation(popPath)
    g = pop[random.randint(0, len(pop) - 1)]

    run[update] = g

  return run


# Determine which tasks specified genotypes can perform
def determineTasks(genotypes):

  global TASKS

  # Initialize task list
  if TASKS == []:
    TASKS = avida.get_tasks()

  # If genotypes is a single genotype, convert to a list
  try:
    iter(genotypes)
  except:
    genotypes = [genotypes]

  # Copy the sequences of each genotype into a new sequence string list
  # because avida.calc_info takes in a list of sequence strings
  seqs = [genotype.sequence for genotype in genotypes]

  # Get task list
  seqs_info = avida.get_info(seqs, ['task_list'])
  task_list = seqs_info['task_list']

  # Update each genotype with its tasks
  for genotype, tasks in zip(genotypes, task_list):
    genotype.setTasks(tasks)


# Return list of named instructions for specified genotype
def getInstructions(genotype):

  INSTS = { 'a': 'nop-A', 'b': 'nop-B', 'c': 'nop-C', 'd': 'if-n-equ', \
    'e': 'if-less', 'f': 'pop', 'g': 'push', 'h': 'swap-stk', 'i': 'swap', \
    'j': 'shift-r', 'k': 'shift-l', 'l': 'inc', 'm': 'dec', 'n': 'add', \
    'o': 'sub', 'p': 'nand', 'q': 'IO', 'r': 'h-alloc', 's': 'divide-sex', \
    't': 'h-copy', 'u': 'h-search', 'v': 'mov-head', 'w': 'jmp-head', \
    'x': 'get-head', 'y': 'if-label', 'z': 'set-flow' }

  insts = []
  for inst in genotype.sequence:
    insts.append(INSTS[inst])

  return insts


# Returns all single mutants
def getSingleMutants(genotype):

  muts = []

  seq = genotype.sequence
  for pos in range(len(seq)):
    for inst in 'abcdefghijklmnopqrstuvwxyz':
      if seq[pos] != inst:
        mutSeq = avida.mutate(seq, pos, inst)
        muts.append(Genotype(mutSeq, genotype))

  return muts


# Returns m mutants with n random mutations
def getMutants(genotype, m, n):

  mutants = []

  for i in range(m):

    # Initialize mutant to-be
    mutant = Genotype(genotype, genotype)

    for j in range(n):
    
      # Get a random position and instruction (except on previous mutation)
      randPos = avida.rand_pos_except(mutant.sequence, mutant.mutations())
      randInst = avida.rand_inst_except(mutant.sequence[randPos])

      # Mutate the genotype
      mutant.sequence = avida.mutate(mutant.sequence, randPos, randInst)

    mutants.append(mutant)

  return mutants


# Returns the gene diversity (proportion of polymorphic loci across genome)
def getGeneDiversity(pop):

  polymorphicLoci = 0
  for c in range(len(pop[0].sequence)):
    inst = pop[0].sequence[c]
    for g in pop:
      if inst != g.sequence[c]:
        polymorphicLoci += 1
        break

  return polymorphicLoci / len(pop[0].sequence)


# Returns the heterozygosity (mean probability of forming a
# heterozygote--hypothetically since these are haploid--per locus)
def getHeterozygosity(pop):
  pass


# Returns the mean number of alleles per locus
def getMeanAllelesPerLocus(pop):
  pass

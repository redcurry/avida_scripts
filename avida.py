import os
import random
import string
import gzip


AVIDA_DIR = '/mnt/home/carlosja/config/snowball'
AVIDA_CONFIG_PATH = AVIDA_DIR + '/avida.cfg'
ANALYZE_CONFIG_PATH = 'analyze.cfg'
ENVIRONMENT_CONFIG_PATH = AVIDA_DIR + '/environment.cfg'
INSTSET_PATH = AVIDA_DIR + '/instset-heads-sex.cfg'
START_CREATURE_PATH = AVIDA_DIR + '/default-heads-sex.org'
EVENT_PATH = AVIDA_DIR + '/events.cfg'


def set_environment(environment_file_name):
  global ENVIRONMENT_CONFIG_PATH
  if environment_file_name != None:
    ENVIRONMENT_CONFIG_PATH = environment_file_name


def reset_environment():
  global ENVIRONMENT_CONFIG_PATH
  ENVIRONMENT_CONFIG_PATH = AVIDA_DIR + '/environment.cfg'


def set_instset(instset_filename):
  global INSTSET_PATH
  INSTSET_PATH = instset_filename


def reset_instset():
  global INSTSET_PATH
  INSTSET_PATH = 'instset-heads.cfg'


# Note: if z == True, assumes compressed file
def read_data(path, detail_names, z = False):

  # Open specified file for reading
  if z:
    file = gzip.open(path, "r")
  else:
    file = open(path, "r")

  # Create a dictionary with each specified detail
  # name as a key, and initialize each entry with
  # an empty list, where the data read will be stored
  data = {}
  for detail_name in detail_names:
    data[detail_name] = []

  for line in file:
    # Read the format (detail) list of the file
    if line.startswith("#format"):

      format_list = line.split()
      format_list.remove("#format")

    # Ignore comments
    elif not line.startswith("#"):

      data_row = line.split()

      # Ignore blank lines
      if len(data_row) > 0:

        for detail_name in detail_names:

          # Get the index within data_row of the desired detail,
          # and add the appropriate data to the data dictionary
          data_row_index = format_list.index(detail_name)
          data[detail_name].append(data_row[data_row_index])

  # Because the data were read as strings, fix their types
  for detail_name in detail_names:
    for i in range(len(data[detail_name])):
      data[detail_name][i] = correct_type(data[detail_name][i], detail_name)

  # Close the file
  file.close()

  return data


# Converts the type of each item in the data list
# according to its detail name (e.g., fitness)
def correct_type(item, detail_name):

  if detail_name == "update":
    return int(item)

  elif detail_name == "fitness":
    return float(item)

  elif detail_name == "viable":
    return int(item)

  else:
    return item


# Returns a list with the requested info (info - list of desired info)
# for each sequence specified (seqs - list of sequences)
def get_info(seqs, info):

  # Create a file to write sequences in
  seq_file = open("seqs.dat", "w")

  # Write header information for Avida analyze mode
  print >> seq_file, "#filetype genotype_data"
  print >> seq_file, "#format sequence"

  # Write sequences to file
  for seq in seqs:
    print >> seq_file, seq

  seq_file.close()

  # Move file to Avida directory
#  os.system('mv seqs.dat ' + AVIDA_DIR)

  # Create analyze file (and move to Avida directory)
  analyze_file = open(ANALYZE_CONFIG_PATH, "w")
  print >> analyze_file, "load seqs.dat"
  print >> analyze_file, "recalc"
  print >> analyze_file, "detail seqs_info.dat sequence " + string.join(info)
  analyze_file.close()
#  os.system('mv ' + ANALYZE_CONFIG_PATH + ' ' + AVIDA_DIR)

  # Run Avida in analyze mode
#  curdir = os.getcwd()
#  os.chdir(AVIDA_DIR)
  os.system('avida -c ' + AVIDA_CONFIG_PATH + \
    ' -set ANALYZE_MODE 1 -set ANALYZE_FILE ' + ANALYZE_CONFIG_PATH +
    ' -set ENVIRONMENT_FILE ' + ENVIRONMENT_CONFIG_PATH +
    ' -set INST_SET ' + INSTSET_PATH +
    ' -set START_CREATURE ' + START_CREATURE_PATH +
    ' -set EVENT_FILE ' + EVENT_PATH +
    ' > /dev/null')

  # Read info (including the sequence itself) as a list
  seqs_info_list = read_data("data/seqs_info.dat", ["sequence"] + info)

  # Delete files
  os.system("rm data/seqs_info.dat")
  os.system("rm seqs.dat")
  os.system("rm analyze.cfg")

  # Go back to working directory
#  os.chdir(curdir)

  return seqs_info_list


# Return a list of tasks based on environment file
def get_tasks():

  env_file = open(os.path.join(AVIDA_DIR, ENVIRONMENT_CONFIG_PATH))

  tasks = []
  for line in env_file:
    if line.startswith('REACTION'):
      tasks.append(line.split()[2])

  return tasks


# Return a random instruction
def rand_inst():
  return chr(random.choice(range(ord('a'), ord('z') + 1)))


# Return a random instruction that is not in the specified list
# Assumption: the specified list doesn't exhaust all instructions
def rand_inst_except(inst_list):

  # Convert to list if single instruction is given
  try:
    iter(inst_list)
  except:
    inst_list = [inst_list]

  r_inst = rand_inst()
  while(r_inst in inst_list):
    r_inst = rand_inst()

  return r_inst


# Return a random position (index) from the specified sequence
def rand_pos(seq):
  return random.randint(0, len(seq) - 1)


# Return a random position that is not in the specified list
# Assumption: the specified list doesn't exhaust all positions
def rand_pos_except(seq, pos_list):

  r_pos = rand_pos(seq)
  while(r_pos in pos_list):
    r_pos = rand_pos(seq)

  return r_pos


# Mutate the specified sequence at the specified
# position with the specified instruction
def mutate(seq, pos, inst):
  return seq[0:pos] + inst + seq[(pos + 1):len(seq)]

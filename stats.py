import random


# Returns a random subset of size n from alist (without replacement)
# Previous function name: randomSubset
def sample(alist, n):

  subsetIndeces = []

  while n > 0:
    index = random.randint(0, len(alist) - 1)
    while index in subsetIndeces:
      index = random.randint(0, len(alist) - 1)
    subsetIndeces.append(index)
    n -= 1

  subset = [alist[i] for i in subsetIndeces]

  return subset


def mean(alist):

  mean = 0.0
  for item in alist:
    mean += item

  if len(alist) > 0:
    mean /= len(alist)

  return mean


def var(alist):

  average = mean(alist)

  sum = 0.0
  for item in alist:
    sum += (item - average) * (item - average)
  sum /= len(alist) - 1

  return sum


def stdev(alist):
  return pow(var(alist), 0.5)


def cov(alist, blist):

  average1 = mean(alist)
  average2 = mean(blist)

  sum = 0.0
  for i in range(len(alist)):
    sum += (alist[i] - average1) * (blist[i] - average2)
  sum /= len(alist) - 1

  return sum

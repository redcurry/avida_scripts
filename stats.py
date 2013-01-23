import random
import numpy


# Returns a random subset of size n from alist (without replacement)
# Previous function name: randomSubset
def sample(alist, n):
  return random.sample(alist, n)


def mean(alist):
  return numpy.mean(alist)


def var(alist):
  return numpy.var(alist)


def stdev(alist):
  return numpy.std(alist)


def cov(alist, blist):
  amean = mean(alist)
  bmean = mean(blist)

  return numpy.sum([(alist[i] - amean) * (blist[i] - bmean)
    for i in range(len(alist))]) / float((len(alist) - 1))

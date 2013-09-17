
import sys
import math
from scipy import stats

# FBAT - ISB implementation.
# david l gibbs
# dgibbs@systemsbiology.org
# June 20, 2013

# Goals:
# single marker fbat test
# rare-variant region-based test
#
 
# Input:
# tped and tfam files as specified by plink.  No recoding needed.
# for now: trios only.
 
# Output:
# tab separated table of test results, p-values

# the test statistic is defined as:
# U = sum_k U_k
# where U_k = sum_ij (Y_k - mu)(X_ijk - E(X_ijk) = T_ij * X_ijk - E(X_ijk) * T_ij
# k indexes the pedigrees
# j is the offspring of the ith nuclear family
#
# Var(U) = sum_i (U_i)^2 for one nuclear family per pedigree
# where U_i is the i-th families contribution to U

class SingleTest:
    """ performs a single test"""
    useidx = [] # an index of the samples that pass validation
      
    def __init__(self, familyList, worker):
        self.Z = 0
        self.familyList = familyList
        self.worker = worker
             
    def test(self):
        """ perform the single marker fbat test """
        # for each family:
        for k in self.familyList:
            # compute the individual stats
            self.worker.singleStats(self.familyList, k)
            self.familyList[k].printfam()



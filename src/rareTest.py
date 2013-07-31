
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
# and additive coding for X .. i.e. the number of A alleles  
#    also expected to be biallelic
 
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

"""
From the Madsen, Browning paper:
The weight, w_i, is the estimated standard deviation of the total number of mutations in the sample (including affected and unaffected individuals), under the null hypothesis of no frequency differences between affected and unaffected individuals. It is used to down-weight mutation counts in constructing the weighted-sum score;

We estimate qi according to the mutation-frequency in the unaffected individuals only, rather than the frequency in the combined population of affected and unaffected individuals. We use this approach so that a true signal from an excess of mutations in the affected individuals is not deflated by using the total number of mutations in both affected and unaffected individuals. By using a permutation-based test, we account for using only the unaffected individuals when scaling the mutation frequency, and we are hence able to increase the power of detecting very rare disease-associated mutations. The drawback of this approach is a higher variance of the scaled mutation-frequency, and hence a loss of power when the frequency of the mutation is high. Adding one to the numerator and two to the denominator of the frequency estimate, qi, avoids zero estimates which would lead to numerical problems in the genetic score used below, and is based on the Bayesian posterior-mean estimate of a binomial proportion when using a uniform prior.
"""

import numpy as np
from math import sqrt
from singleTest import SingleTest

class RareTest:
    """ performs a single test"""
    def __init__(self, region, markers, chrms, freqs, listofgs, ys, famidx, childidx, paridx, weighted):
        """ any initial tasks """
        self.region = region # the name of this region
        self.gs = listofgs # for each marker .. the set of genotypes
        self.ys = ys # the phenotypes for this marker, adjusted for offset == T_ij
        self.chrms = chrms
        self.famidx = famidx # the index into families.
        self.childidx = childidx # index of children within families
        self.paridx = paridx
        self.markers = markers # maps to each of the genotype sets
        self.freqs = freqs # the allele freqs for each marker
        self.markerfreq = []
        self.freqdict = dict()
        self.weighted = weighted
        self.weights = [] # the weights to use ... mapping to each marker... from allele freq of each marker.
        self.M = len(markers)
        self.Um =[] # the list of statistics from each marker in the region
        self.W =[] # the weighted sum of the Us
        self.D = np.array(0) # the diagonal matrix with elements of Um, square rooted
        self.Ve = np.array(0) # the covariance matrix between all Um's
        self.Va = np.array(0) # the 
        self.VarW = 0
        self.Z = 0
        self.pvalue = 1
        self.qvalue = 1

    def buildFreqDict(self):
        for x in self.freqs:
            y = x.split("\t")
            self.freqdict[y[0]] = float(y[1])

    def kaviarFilter(self):
        # filter the data for only markers in kaviar.
        1==1
        
    def computeWeights(self):
        """ the weighting methods """
        # Madsen and Browning method #
        if self.weighted == "madsen":
            w = []
            for m in self.markers:   # for each marker
                i = self.markers.index(m)
                # get number of people with X > 0
                unaffected = []
                for j in range(len(self.ys)):
                    if self.ys[j] <= 0:
                        unaffected.append(self.gs[i][2*j-1])
                        unaffected.append(self.gs[i][2*j])
                mark = self.Um[i].markers[len(self.Um[i].markers)-1]
                mu = sum([1.0 for k in unaffected if k == mark])/2
                nu = len(unaffected)/2
                ni = len(self.ys)
                q = (mu+1.0)/(2.0*nu+2.0)
                # what if mu and nu are 0?
                w.append(1/(sqrt(ni*q*(1-q))))
            self.weights = map(lambda x: x/sum(w), w)
        elif self.weighted == "kaviar": # here we assume the data has been filtered
            w = []                      # so that only markers used are in Kaviar
            for m in self.markers:      # for each marker
                a = self.freqdict[m]    # get the allele freq
                if a < 1:
                    w.append(1/sqrt(len(self.ys)*a*(1-a)))
                else:
                    w.append(0)
            self.weights = w
        else:
            self.weights = [1 for i in range(len(self.Um))]

    def computeSingleStats(self):
        """ a list of single tests """
        self.Um = [SingleTest(self.markers[i], self.chrms[i], self.gs[i],
                              self.ys, self.famidx, self.childidx,
                              self.paridx, False, 0.0) for i in range(len(self.gs))]
        for i in range(len(self.Um)):
            self.Um[i].test(False)
            self.markerfreq.append(self.Um[i].allelefreq)
        self.markers = [self.markers[i] for i in range(len(self.Um)) if self.Um[i].U != -1]
        self.gs = [self.gs[i] for i in range(len(self.Um)) if self.Um[i].U != -1]
        self.chrms = [self.chrms[i] for i in range(len(self.Um)) if self.Um[i].U != -1]
        self.Um = [x for x in self.Um if x.U != -1] # remove any that are empty
        self.M = len(self.Um)
        
    def computeW(self):
        """ the rare test statistic """
        usum = map(lambda x: sum(x), [x.U for x in self.Um])
        self.W = map(lambda x, y: x*y, self.weights, usum)  # W =  w_m * U_m

    def computeD(self):
        """ diagonal matrix with entries sqrt(Var(U_m)) """
        self.D = np.array([0 for i in range(pow(self.M,2))], dtype='f').reshape(self.M, self.M)
        for i in range(self.M):
            x = sum(self.Um[i].V)
            if x > 0:
                self.D[i][i] = sqrt(x)
        
    def computeVe(self):
        """ MxM matrix """
        self.Ve = np.array([0 for i in range(pow(self.M,2))], dtype='f').reshape(self.M, self.M)
        cidx = map(lambda x,y: x+y[0], self.Um[0].famidx, self.Um[0].childidx)
        for p in range(self.M):
            for q in range(self.M): # for each pair of markers #
                e = 0
                wp = self.weights[p]
                wq = self.weights[q]
                for k in range(len(self.Um[0].X)): # for each offspring (family)
                    t = pow(self.Um[0].ts[cidx[k]],2) # all have same T_ij .. get kids phenotype T=Y-mu
                    #a = wp * t * self.Um[p].X[k] - self.Um[p].EofX[k]
                    #b = wq * t * self.Um[q].X[k] - self.Um[q].EofX[k]
                    a = self.Um[p].X[k] - self.Um[p].EofX[k]
                    b = self.Um[q].X[k] - self.Um[q].EofX[k]
                    e += t*a*b
                self.Ve[p][q] = e

    def diagonalizeRoot(self, A):
        """ take a matrix, return a matrix
        with only diagonal entries, 1/square rooted """
        B = np.copy(A)
        for i in range(B.shape[0]):
            for j in range(B.shape[1]):
                if i == j and B[i][j] > 0:
                    B[i][j] = 1/sqrt(B[i][j])
                else:
                    B[i][j] = 0
        return(B)
        
    def computeVa(self):
        diagVe = self.diagonalizeRoot(self.Ve)
        x = diagVe.dot(self.Ve).dot(diagVe)
        self.Va = self.D.dot(x).dot(self.D)

    def computeVarW(self):
        oners = np.array([1 for i in range(self.Va.shape[0])], dtype='f')
        self.VarW = oners.transpose().dot(self.Va).dot(oners)
                
    def computeZ(self):
        """ the main statistic """
        self.Z = sum(self.W) / sqrt(self.VarW)

    def computePvalues(self):
        """ compute either pvalues using either ChiSq or Z """
        ''' Z: abs(Z), lower.tail, *2 for two sided test.'''
        f = stats.norm()
        self.pvalue = (1-f.cdf(abs(self.Z))) * 2

    def printTest(self):
        #for i in range(len(self.Um)):
        #    self.Um[i].printTest()
        print(self.chrms[0] + "\t" + self.region + "\t" + str(sum(self.W)) + "\t"
              + str(self.VarW) + "\t" + str(self.Z)
              + "\t" +  str(self.pvalue) + "\t" +
              str(len(self.Um)) + "\t" + str(self.weights))
        #print("weights: " + str(self.weights))
        #print("W: " + str(self.W))
        #print("W: " + str(sum(self.W)))
        #print("Ve: " + str(self.Ve))
        #print("D: " + str(self.D))
        #print("VarW: " + str(self.VarW))
        #print("Z: " + str(self.Z))
        #print("pvalue: " + str(self.pvalue))

    def test(self):
        """ perform the single marker fbat test """
        if self.freqs != []:
            self.buildFreqDict()
        self.computeSingleStats()
        if len(self.Um) > 1:
            self.computeWeights()
            self.computeW()
            self.computeD()
            self.computeVe()
            self.computeVa()
            self.computeVarW()
            self.computeZ()
            self.computePvalues()
            self.printTest()
            

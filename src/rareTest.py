
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

import numpy as np
from math import sqrt
from singleTest import SingleTest

class RareTest:
    """ performs a single test"""
    def __init__(self, region, markers, freqs, listofgs, ys, famidx, childidx, paridx, weighted):
        """ any initial tasks """
        self.region = region # the name of this region
        self.gs = listofgs # for each marker .. the set of genotypes
        self.ys = ys # the phenotypes for this marker, adjusted for offset == T_ij
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
        if self.weighted == "data":
            w = []
            for m in self.markers:   # for each marker
                i = self.markers.index(m)
                # get number of people with X > 0
                unaffected = []
                for j in range(len(self.ys)):
                    if self.ys[i] <= 0:
                        unaffected.append(self.gs[i][2*j-1])
                        unaffected.append(self.gs[i][2*j])
                mark = self.Um[i].markers[len(self.Um[i].markers)-1]
                mu = sum([1.0 for j in unaffected if j == mark])
                nu = len(unaffected)
                q = (mu+1.0)/(2.0*nu+2.0)
                w.append(1/sqrt(len(self.ys)*q*(1-q)))
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
            self.weights = [1 for i in range(len(self.markers))]

    def computeSingleStats(self):
        """ a list of single tests """
        self.Um = [SingleTest(self.markers[i], self.gs[i],
                              self.ys, self.famidx, self.childidx,
                              self.paridx) for i in range(len(self.gs))]
        for i in range(len(self.Um)):
            self.Um[i].test(False)
            self.markerfreq.append(self.Um[i].allelefreq)
        
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
            else:
                self.D[i][i] = 0
        
    def computeVe(self):
        """ MxM matrix """
        self.Ve = np.array([0 for i in range(pow(self.M,2))], dtype='f').reshape(self.M, self.M)
        cidx = map(lambda x,y: x+y[0], self.Um[0].famidx, self.Um[0].childidx)
        for p in range(self.M):
            for q in range(self.M): # for each pair of markers
                e = 0
                wp = self.weights[p]
                wq = self.weights[q]
                for k in range(len(self.Um[0].X)): # for each offspring (family)
                    t = pow(self.Um[0].ts[cidx[k]],2) # all have same T_ij .. get kids phenotype T=Y-mu
                    #a = wp * t * self.Um[p].X[k] - self.Um[p].EofX[k]
                    #b = wq * t * self.Um[q].X[k] - self.Um[q].EofX[k]
                    a = t * self.Um[p].X[k] - self.Um[p].EofX[k]
                    b = t * self.Um[q].X[k] - self.Um[q].EofX[k]
                    e += a*b
                self.Ve[p][q] = e

    def diagonalizeRoot(self, A):
        """ take a matrix, return a matrix
        with only diagonal entries, square rooted """
        B = np.copy(A)
        for i in range(B.shape[0]):
            for j in range(B.shape[1]):
                if i == j:
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
        x = stats.norm()
        self.pvalue = x.pdf(self.Z)

    def printTest(self):
        #for i in range(len(self.Um)):
        #    self.Um[i].printTest()
        print(self.region + "\t" + str(sum(self.W)) + "\t"
              + str(self.VarW) + "\t" + str(self.Z)
              + "\t" +  str(self.pvalue) + "\t" +
              str(self.markerfreq) + "\t" + str(self.weights))
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
        self.buildFreqDict()
        self.computeSingleStats()
        self.computeWeights()
        self.computeW()
        self.computeD()
        self.computeVe()
        self.computeVa()
        self.computeVarW()
        self.computeZ()
        self.computePvalues()
        self.printTest()
            

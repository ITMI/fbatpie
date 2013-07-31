
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

class SingleTest:
    """ performs a single test"""
    useidx = [] # an index of the samples that pass validation
      
    def __init__(self, marker, chrm, gs, ys, famidx, childidx, paridx, silent, freqcutoff):
        """ any initial tasks """		      
        self.gs = gs # the genotypes for this marker
        self.ts = ys # the phenotypes for this marker, adjusted for offset == T_ij
        self.famidx = famidx # the index into families.
        self.childidx = childidx # index of children within families
        self.paridx = paridx
        self.silent = silent # do we print the mendelian errors?
        self.freqcutoff = freqcutoff
        self.chrm = chrm
        self.markers = list(set(gs))
        self.markers.sort() # markers will go [0,1,2] or [0,A,B] etc.
        self.thismarker = marker
        self.markerCount = [0,0,0]
        self.allelefreq = 0
        self.X = []
        self.S = []     # sum_j of X_ij * T_ij
        self.EofX = []  # E(X_ij) 
        self.E = []     # E(S_i) = E(X_ij) * T_ij
        self.U = []     # U = S_i - E(S_i)
        self.VofX = []  # Var(X_ij)
        self.V = []     # V(S_i) = sum_j Var(X_ij) * T_ij
        self.fp = FbatProb()
        self.pvalue = 1
        self.qvalue = 1
        self.Z = 0
            
        # gs ['a', 'b', 'a', 'a', 'b', 'b', ... ]
        # idx  0    1    2    3    4    5   
        # ys [ 1, 2, 1, ... ]
        # idx  0  1  2
        # famidx [0, 3, .. ]
        # cidx = [ 1 ]
        # pidx = [[0,2], ... ]

    def markerSet(self):
        if len(self.markers) == 1:
            self.markers = ['0', self.markers[0]]
        elif len(self.markers) < 3 and self.markers[0] != '0':
            self.markers = ['0', self.markers[0], self.markers[1]]
        if len(self.markers) > 2:
            self.markerCount[0] = self.gs.count(self.markers[0])
            self.markerCount[1] = self.gs.count(self.markers[1])
            self.markerCount[2] = self.gs.count(self.markers[2])

    def validate (self):
        """ look for problems in inheritance - log errors """
        for i in range(len(self.famidx)):             # for each family
            idx  = self.famidx[i]                     # family idx
            cidx = map(lambda x:x+idx, self.childidx[i])     # get the child index e.g. [1]
            pidx = map(lambda x:x+idx, self.paridx[i])       # get the parents index e.g. [0,2]
            for j in cidx: # for each child in this family
                cg = [self.gs[2*j], self.gs[2*j+1]] # child genotypes  ### EXPECTING ONE CHILD ###
                pg1 = [self.gs[2*pidx[0]], self.gs[2*pidx[0]+1]]
                pg2 = [self.gs[2*pidx[1]], self.gs[2*pidx[1]+1]]
                if not(((cg[0] == pg1[0] or cg[0] == pg1[1]) and 
                    (cg[1] == pg2[0] or cg[1] == pg2[1])) or
                    ((cg[0] == pg2[0] or cg[0] == pg2[1]) and 
                    (cg[1] == pg1[0] or cg[1] == pg1[1]))):
                    # not valid!
                    if self.silent == True:
                        sys.stderr.write("warning: family: " + str(i) + " has a mendelian error")
                    cg[0] = '0' # zero it out
                    cg[1] = '0'

    def Xfun(self, g):
        return(self.markers.index(g)-1)
        
    def Xfunp(self, gs):
        return([(self.markers.index(gs[0])-1),(self.markers.index(gs[1])-1)])

    def computeAlleleFreq(self):
        # estimated from parents 
        mark = self.markers[len(self.markers)-1] # check the last marker
        a = 0.0 # the number of alleles in parents
        b = 0.0 # total number minus missing
        for i in self.famidx: # i is the index into the ts
            for j in self.paridx[self.famidx.index(i)]: # j is the offset from i
                for k in range(2):
                    if self.gs[(2*(j+i))+k] == mark:
                        a += 1.0
                    if self.gs[2*(i+j)+k] != '0':
                        b += 1.0
        self.allelefreq = float(a)/float(b)
        if self.allelefreq > 0.5 and len(self.markers) > 2: # then the "alt" genotype is the major allele
            m = [self.markers[0],self.markers[2],self.markers[1]]
            self.markers = m
            self.allelefreq = 1-self.allelefreq
        
    def computeS(self):
        """ S_i = sum over children j in family i of X_ij * T_ij """
        for i in self.famidx: # i is the index into the ts
            for j in self.childidx[self.famidx.index(i)]: # j is the offset from i
                g1 = self.gs[2*(j+i)]
                g2 = self.gs[2*(j+i)+1]
                if g1 != '0' and g2 != '0':
                    self.X.append((self.Xfun(g1) + self.Xfun(g2)))
                    self.S.append((self.Xfun(g1) + self.Xfun(g2))* self.ts[i+j])
                else:
                    self.X.append(0)
                    self.S.append(0)
        
    def computeEofXandV(self):
        """ with parent genotypes ... compute E(X_ij) """
        for i in range(len(self.famidx)):             # for each family
            idx  = self.famidx[i]
            cidx = map(lambda x:x+idx, self.childidx[i])     # get the child index [1]
            pidx = map(lambda x:x+idx, self.paridx[i])       # get the parents index [0,2]
            for j in cidx: # for each child in this family
                cg = [self.gs[2*j], self.gs[2*j+1]] # child genotypes  ### EXPECTING ONE CHILD ###
                pg1 = [self.gs[2*pidx[0]], self.gs[2*pidx[0]+1]]
                pg2 = [self.gs[2*pidx[1]], self.gs[2*pidx[1]+1]]
                self.EofX.append(self.fp.expLookup(pg1, pg2, cg, self.markers))
                self.VofX.append(self.fp.varLookup(pg1,pg2,cg,self.markers) - pow(self.fp.expLookup(pg1,pg2,cg,self.markers),2))
                        
    def computeUandV(self):
        """ U = sum over families i of S_i - E(S_i) 
        E(S_i) = sum over children in family i of E(X_ij) * T_ij
        E(X_ij) = sum over genotypes possible in family i of X(g) * P(g)
        where X(g) == number of A alleles,
        and P(g) is the probability given the parental genotypes"""
        childT = map(lambda x,y: x+y[0], self.famidx, self.childidx)
        self.E = map(lambda x,y:x*y, self.EofX, [self.ts[i] for i in childT])  # E(S_i) =  E(X_ij) * T_ij
        self.U = map(lambda x,y:x-y, self.S, self.E)                           # U = sum_i S_i - E(S_i)
        self.V = map(lambda x,y:x*y, self.VofX, [pow(self.ts[i],2) for i in childT] )  # V(S_i) = Var(X_ij) * T_ij^2
        self.famN = len(filter(lambda x: x != 0, self.U))                      # number of families with non-zero Us
        
        
    def computePvalue(self):
        """ compute either pvalues using either ChiSq or Z """
        ''' Z: abs(Z), lower.tail, *2 for two sided test.'''
        f = stats.norm()
        self.pvalue = (1-f.cdf(abs(self.Z))) * 2
             
    def test(self, printit):
        """ perform the single marker fbat test """
        self.markerSet()
        self.validate()
        self.computeAlleleFreq()
        self.computeS()
        self.computeEofXandV()
        self.computeUandV()
        if len(self.markers) > 2 and self.allelefreq > self.freqcutoff and sum(self.V) > 0:
            if sum(self.V) == 0 or sum(self.U) == 0:
                self.Z = 0
            else:
                self.Z = sum(self.U) / math.sqrt(sum(self.V) + 0.0000001)
            self.computePvalue()
        else:
            self.E = -1
            self.V = -1
            self.U = -1
            self.famN = -1
        if(printit == True and self.E != -1):
            self.printTest()

    def printTest(self):
        childT = map(lambda x,y: x+y[0], self.famidx, self.childidx)
        print(self.chrm + "\t" + self.thismarker + "\t" + str(self.markers) + "\t" + str(self.markerCount) + "\t" 
              + str(self.allelefreq) + "\t" + str(self.famN) + "\t"
              + str(sum(self.U)) + '0' + "\t" + str(sum(self.V)) + "\t" +  str(self.Z) + "\t" + str(self.pvalue))
        #print("T:    " + str([self.ts[i] for i in map(lambda x,y: x+y[0], self.famidx, self.childidx)]))
        #print("X:    " + str(self.X))
        #print("EofX: " + str(self.EofX))
        #print("V:    " + str(self.V))
        #print()
        
class FbatProb:
    """
    used for finding the expectation
    
    As example:
    Parental genotypes AB AB
    A  B
    A  0  1
    B  1  2
    
    P(AA) = 1/4  P(AB) = 1/2  P(BB) = 1/4
    
    E(X|P) = P(AA) * 0 + P(AB) * 1 + P(BB) + 2 = 1
    
    Var(X|P) = E(X^2) - (E(x))^2
    = [P(AA) * 0^2 + P(AB) * 1^2 + P(BB) * 2^2] - (1)^2 = 0.5 
    
    To use the tables of probabilities for missing parental information,
    have to have at least one complete parent.
    
    When we have just a trio, with missing parental information,
    the rule is the observed child genotype has probability 1
    Say the child has genotype AB.
    then E(X_ij) = P(AA) * X(AA) + P(AB) * X(AB) + P(BB) * X(BB) = 1
                 =   0   *  0    +  1    *  1    +  0    *   2   = 1
    
    With more than one child, we can make inferences on the missing parental
    genotypes. These come from the tables found in fbat_technical_guide
    
    """
    
    def __init__(self):
        """ assuming perfect trios, with parental genotypes in a single list
        coded as 0,1,2 ... where the 2 genotype is counted.
        These are the expected number of 2 genotypes ... """		
        self.exp = dict()
        self.exp[ str([0,0,0,0])] = 0
        self.exp[ str([1,1,1,1])] = 0
        self.exp[ str([1,1,1,2])] = 0.5
        self.exp[ str([1,1,2,1])] = 0.5
        self.exp[ str([2,1,1,1])] = 0.5
        self.exp[ str([1,2,1,1])] = 0.5
        self.exp[ str([1,1,2,2])] = 1
        self.exp[ str([2,2,1,1])] = 1
        self.exp[ str([1,2,1,2])] = 1
        self.exp[ str([1,2,2,1])] = 1
        self.exp[ str([2,1,1,2])] = 1
        self.exp[ str([2,1,2,1])] = 1
        self.exp[ str([1,2,2,2])] = 1.5
        self.exp[ str([2,1,2,2])] = 1.5
        self.exp[ str([2,2,2,1])] = 1.5
        self.exp[ str([2,2,1,2])] = 1.5
        self.exp[ str([2,2,2,2])] = 2

        # The variance of the allele count given the parents
        self.var = dict()
        self.var[ str([0,0,0,0])] = 0
        self.var[ str([1,1,1,1])] = 0
        self.var[ str([1,1,1,2])] = 0.5   #   1  1
        self.var[ str([1,1,2,1])] = 0.5   # 1 0  0  => 0.5
        self.var[ str([2,1,1,1])] = 0.5   # 2 1  1
        self.var[ str([1,2,1,1])] = 0.5
        self.var[ str([1,1,2,2])] = 1.0 
        self.var[ str([2,2,1,1])] = 1.0 
        self.var[ str([1,2,2,1])] = 1.5 
        self.var[ str([2,1,1,2])] = 1.5 
        self.var[ str([1,2,1,2])] = 1.5   
        self.var[ str([2,1,2,1])] = 1.5
        self.var[ str([1,2,2,2])] = 2.5
        self.var[ str([2,1,2,2])] = 2.5
        self.var[ str([2,2,2,1])] = 2.5
        self.var[ str([2,2,1,2])] = 2.5        
        self.var[ str([2,2,2,2])] = 4

        # the expectation and variance when
        # parent information is missing
        self.expmissing = dict()
        self.expmissing[ str([1,1]) ] = 0
        self.expmissing[ str([1,2]) ] = 1
        self.expmissing[ str([2,1]) ] = 1
        self.expmissing[ str([2,2]) ] = 2

        self.varmissing = dict()
        self.varmissing[ str([1,1]) ] = 0
        self.varmissing[ str([1,2]) ] = 1
        self.varmissing[ str([2,1]) ] = 1
        self.varmissing[ str([2,2]) ] = 4
        
    def toidx(self, a):
        return(map(lambda x: self.markers.index(x), a))
    
    def toints(self, a):
        return(map(lambda x: int(x), a))
            
    def expLookup(self,p1,p2,cg,markers):
        self.markers = markers  
        if cg[0] == '0' or cg[1] == '0':  # if the child is missing
            return(0) 
        if (p1[0] == '0' or p1[1] == '0' or
            p2[0] == '0' or p2[1] == '0'):  # if a parent is missing
            return(self.expmissing[str(self.toidx(cg))])
        return(self.exp[str(self.toidx(p1+p2))])
      
    def varLookup(self,p1,p2,cg,markers):
        if cg[0] == '0' or cg[1] == '0':  # if the child is missing
            return(0) 
        if (p1[0] == '0' or p1[1] == '0' or
            p2[0] == '0' or p2[1] == '0'):  # if a parent is missing
            return(self.varmissing[str(self.toidx(cg))])
        return(self.var[str(self.toidx(p1+p2))]) #otherwise normal ops



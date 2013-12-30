from math import pow
from math import sqrt
from itertools import product
from scipy import stats
import numpy as np
import sys
    
class Computer:

    def __init__(self, mode, offset, afreq):
        self.mode = mode
        self.testtype = "single"
        self.offset = float(offset)
        self.afreq = float(afreq)
        self.markers = []
        self.U = 0
        self.VarU = 0
        self.Z = 0
        self.M = 0
        self.W = 0
        self.VarW = 0

    def update(self, markers):
        # the marker_name, chrm, and genotypes
        self.markers = markers
        self.U = 0
        self.VarU = 0
        self.Z = 0
        
    def adjPhenotypes(self,pi):
        # fbat takes ped format: 2 == affected, 1 == unaffected
        # and codes them as 1 and 0 respectively.
        if pi == '2':
            return('1')
        elif pi == '1':
            return('0')
        else:
            return('-1') # NA
        
    def computeT(self, family):
        pheno = self.adjPhenotypes(family.phenotypes['NB'])
        # then apply offset
        if pheno != '-1':
            return(float(pheno)-float(self.offset))
        else:
            return(-9)
            
    def markerIndex(self, g):
        return(self.markers.index(g)-1.0)

    def pow2(self, s):
        return(pow(s,2))

    def computeX(self, gs):
        if self.mode == "additive":
            return(sum([self.markerIndex(gs[0]),self.markerIndex(gs[1])]))
        elif self.mode == "recessive":
            # if both are the minor allele
            if self.markerIndex(gs[0]) == 1 and self.markerIndex(gs[1]) == 1:
                return(1)
            else:
                return(0)
        elif self.mode == "dominant":
            # if either are the minor allele
            if self.markerIndex(gs[0]) == 1 or self.markerIndex(gs[1]) == 1:
                return(1)
            else:
                return(0)
        else:
            sys.stderr.write("ERROR!! incorrect mode specification [additive, recessive, dominant]\n")
            sys.exit(1)

    def computeEs(self, family):
        if '0' in family.genotypes['M'] or '0' in family.genotypes['F']:
            # then one of the parents has missing data, so use NB for E and E^2
            return([self.computeX(family.genotypes['NB']),
                    pow(self.computeX(family.genotypes['NB']),2)])
        else:
            # else return the expected number of alleles from parents
            # which is a count on the combinations of genotypes possible
            genotypes = product(family.genotypes['M'], family.genotypes['F'])
            scores = [self.computeX(g) for g in genotypes]                
            return([(sum(scores)/4.0), (sum(map(self.pow2,scores))/4)]) 
                
    def computeEsX(self, family):
        gender = family.genders['NB']
        if gender == '1': # it's a boy
            if '0' in family.genotypes['M']:
                # then the mother has missing data so use NB for E and E^2
                return([self.computeX(family.genotypes['NB']),
                        pow(self.computeX(family.genotypes['NB']),2)])
            else:
                if self.mode == "additive":
                    # else return the expected number of alleles from parents
                    # which is a count on the combinations of genotypes possible
                    score = self.computeX(family.genotypes['M'])
                    return([(score/2.0), (score/2.0)])
                else: # it's dominant...
                    g1 = [family.genotypes['M'][0], family.genotypes['M'][0]]
                    g2 = [family.genotypes['M'][1], family.genotypes['M'][1]]
                    genotypes = [g1,g2]
                    scores = [self.computeX(g) for g in genotypes]                
                    return([(sum(scores)/2.0), (sum(map(self.pow2,scores))/2.0)])
                        
        else: # it's a girl
            if ('0' in family.genotypes['M'] or '0' in family.genotypes['F']):
                # then the mother has missing data so use NB for E and E^2
                return([self.computeX(family.genotypes['NB']),
                        pow(self.computeX(family.genotypes['NB']),2)])
            else:
                # else return the expected number of alleles from parents
                # which is a count on the combinations of genotypes possible
                genotypes = product(family.genotypes['M'], [family.genotypes['F'][0]])
                scores = [self.computeX(g) for g in genotypes]        
                return([(sum(scores)/2.0), (sum(map(self.pow2,scores))/2.0)])
            
    def computeU(self, family):
        return(family.stat['T'] * (family.stat['X'] - family.stat['E']))
            
    def singleStats(self, familyList, k, marker):
        T = 0; X = 0; E = 0; E2 = 0; V = 0

        # first T
        T = self.computeT(familyList[k])
        familyList[k].updateStat('T', T)
        chrm = familyList[k].curChr
        if T != -9 and '0' not in familyList[k].genotypes['NB']:

            # then X
            X = self.computeX(familyList[k].genotypes['NB'])
            if 'X' in chrm and familyList[k].genders['NB'] == '1' and self.mode == "additive":
                X = 0.5*X
            familyList[k].updateStat('X', X)

            # then E(X|P)
            if 'X' in chrm:
                E,E2 = self.computeEsX(familyList[k])
                familyList[k].updateStat('E', E)
                familyList[k].updateStat('E2', E2)
            else:
                E,E2 = self.computeEs(familyList[k])
                familyList[k].updateStat('E', E)
                familyList[k].updateStat('E2', E2)
                
            # then Var = V = E(X^2) - [E(X)]^2
            V = E2-pow(E,2)
            familyList[k].updateStat('V', V)

            # then U = T * (X - E(x|P))
            U = self.computeU(familyList[k])
            familyList[k].updateStat('U', U)

            # and Var(U) = T^2 * V
            VarU = pow(T,2)*V
            familyList[k].updateStat('VarU',VarU)

            # we accumulate across families
            self.U += U
            self.VarU += VarU
        else:
            #familyList[k].updateStat('T', 0)
            familyList[k].updateStat('X', 0)
            familyList[k].updateStat('E', 0)
            familyList[k].updateStat('V', 0)
            familyList[k].updateStat('U', 0)
            familyList[k].updateStat('VarU',0)
            familyList[k].hasHet = 0
        if self.testtype == "region":
            self.updateMarkerDict(marker,k,T,X,E)

    def computeZ(self):
        if self.VarU > 0:
            self.Z = self.U/sqrt(self.VarU)
        else:
            self.Z = 0

    def computePValue(self):
        """ compute either pvalues using either ChiSq or Z """
        ''' Z: abs(Z), lower.tail, *2 for two sided test.'''
        f = stats.norm()
        return((1-f.cdf(abs(self.Z))) * 2.0)

    def cumStats(self):
        self.computeZ()
        pvalue = self.computePValue()
        return([self.U, self.VarU, self.Z, pvalue])

    #                     ^
    # single marker tests |
############################################################################
    # region test         |
    #                    \/ 

    
    def setupRegion(self):
        self.testtype = "region"
        self.markerdict = dict()

    def updateRegionM(self, M):
        self.M = int(M)

    def computeD(self, vs):
        #diagonal matrix with entries sqrt(Var(U_m))
        self.D = np.array([0 for i in range(int(pow(self.M,2)))], dtype='f').reshape(self.M, self.M)
        for i in range(self.M):
            if vs[i] > 0:
                self.D[i][i] = sqrt(vs[i])

    def updateMarkerDict(self,marker,k,T,X,E):
        if marker in self.markerdict:
            tdict, xdict, edict = self.markerdict[marker]
            tdict[k] = T
            xdict[k] = X
            edict[k] = E
            self.markerdict[marker] = [tdict,xdict,edict] 
        else:
            tdict = dict()
            xdict = dict()
            edict = dict()
            tdict[k] = T # family k's stats
            xdict[k] = X
            edict[k] = E
            self.markerdict[marker] = [tdict,xdict,edict] 

    def computeVeEntry(self, markerp, markerq):
        tp, xp, ep = self.markerdict[markerp]
        tq, xq, eq = self.markerdict[markerq]
        res0 = 0
        up = 0
        uq = 0
        for k in tp:
            res0 += tp[k]*tq[k]*(xp[k]-ep[k])*(xq[k]-eq[k])
        return(res0)

    def computeVe(self, dat, familyList, markers):
        # MxM matrix taking the U and VarU statistics from each marker
        self.Ve = np.array([0 for i in range(int(pow(self.M,2)))], dtype='f').reshape(self.M, self.M)
        for p in range(self.M):
            for q in range(self.M): # for each pair of markers #
                self.Ve[p][q] = self.computeVeEntry(markers[p], markers[q])

    def diagonalizeRoot(self, A):
        # take a matrix, return a matrix
        #with only diagonal entries, 1/square rooted
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
                
    def computeRegionZ(self):
        # the main statistic
        if self.VarW > 0:
            self.Z = self.W / sqrt(self.VarW)
        else:
            self.Z = 0
            
    def computeRegionStats(self, idx, markers, us, vs, dat, familyList):
        # idx: the index to the markers that are used here
        # Us: the U stats for each marker
        # Vs: the Var stat for each marker
        # dat: the data for each marker
        self.W = sum(us)
        self.computeD(vs)
        self.computeVe(dat, familyList, markers)
        self.computeVa()
        self.computeVarW()
        self.computeRegionZ()
        pvalue = self.computePValue()
        return([self.W,self.VarW,self.Z,pvalue])

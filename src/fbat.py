
# FBAT - ISB implementation.
# david l gibbs
# dgibbs@systemsbiology.org
# June 20, 2013

# Goals:
# single marker fbat test
# rare-variant region-based test
# other variants .. how about log. reg.?
 
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
from itertools import product

import sys
import math
import fileloader
import computer
import family
import validator
import numpy as np
import gzip
from math import sqrt

class FBAT:
    """ the fbat test object """
    
    def __init__(self, analysis, tfam, tped, mode, offset, afreq, nfams, verbose):
        """ any initial tasks """
        self.analysis = analysis
        self.loader = fileloader.Loader(tfam, tped)
        self.worker = computer.Computer(mode, offset, afreq)
        self.familyList = self.loader.buildFamilyList()
        self.validator = validator.Validator(verbose)
        self.afreq = afreq
        self.nfams = nfams
        #for key in self.familyList:
        #    self.familyList[key].printfam()
        self.printFBAT(mode)

    def printFBAT(self, m):
        """ print a little statement about the object """
        print("#FBAT")
        print("#files: " + self.loader.tfamfile + "  " + self.loader.tpedfile)
        print("#Number of families: " + str(len(self.familyList)))
        print("#Offset: " + str(self.worker.offset))
        print("#Mode: " + m)
        

#####################################################################################
    
    def single(self):
        """ perform the single marker test """
        print("AnalysisName\tChr\tMarker\tAllele\tFamilies\tAlleleFreq\tU\tVar(U)\tZ\tP")
        x = self.loader.getMarker()
        while x != '':
            # update family list
            self.loader.updateFamilyList(x, self.familyList)
            # validate the marker
            self.validator.validateMarker(self.familyList)
            # order the markers [missing, major, minor]
            markers, markerCounts = self.validator.markerCheck(x[2], self.familyList)
            # count the number of "informative families"
            info = self.validator.countInfoFams(markers, self.familyList)
            if len(markers) > 2:
                if (info[3] > self.afreq) and (int(info[0]) >= int(self.nfams)) :
                    # update the worker with alleles etc.
                    self.worker.update(markers)
                    # then perform the single marker fbat test #
                    U = 0; VarU = 0;
                    for k in self.familyList:
                        # compute the individual stats
                        self.worker.singleStats(self.familyList, k, '')
                    cumulativeStats = self.worker.cumStats()                    
                    print("\t".join(map(str, [self.analysis,x[1],x[0],markers,str(info[0])+":"
                                              +str(info[2]),info[3]] + cumulativeStats)))

            # AND GET THE NEXT MARKER TO TEST # 
            x = self.loader.getMarker()
            

###########################################################################################

    def regionsingle(self, x, markername):
        # update family list x = [markers, chrm, alleles]
        self.loader.updateFamilyList(x, self.familyList)
        # validate the marker
        self.validator.validateMarker(self.familyList)
        # order the markers [missing, major, minor]
        markers, markerCounts = self.validator.markerCheck(x[2], self.familyList)
        if len(markers) > 2:
            # thisfreq = float(markerCounts[2])/ float(sum(markerCounts))
            # update the worker with alleles etc.
            self.worker.update(markers)
            # then perform the single marker fbat test #
            for k in self.familyList:
                # compute the individual stats
                self.worker.singleStats(self.familyList, k, markername)
            cumulativeStats = self.worker.cumStats() #U\tVar(U)\tZ\tP
            # then we need the var across region
            return(cumulativeStats[0:2])
        else:
            return([0,0])

    def region(self, regionfile, indexfile):
        # first we need the index to the data #
        if indexfile == "none" or indexfile == '':
            indexfile = self.loader.writeIndex()
        tpedindex = self.loader.loadIndex(indexfile)
        # then we get the list of regions # 
        print("Analysis\tChr\tRegion\tW\tVarW\tZ\tpvalue")
        regions = open(regionfile,'r').read().strip().split("\n")
        regions = map(lambda x: x.strip(), regions)
        regionstats = []
        for r in regions:
            try:
                self.worker.setupRegion()
                vlist = []
                ulist = []
                rs = r.split("\t")
                regionname = rs[0]
                theseMarkers = rs[1:]
                # try to find them in the data index ...
                idx = [i for i in range(len(tpedindex[0])) if tpedindex[0][i] in theseMarkers]
                jdx = tpedindex[1][idx]
                markerset = [tpedindex[0][i] for i in range(len(tpedindex[0])) if tpedindex[0][i] in theseMarkers]
                # and if we found them .. which we should ...
                if len(jdx) > 0:
                    # then make the list of data using the file offsets .. call it thisData
                    markers, chrms, thisData = self.loader.buildDataSet(jdx)
                    # and for each marker in the list of data
                    for i in range(len(thisData)):
                        dat = [markers[i], chrms[i], thisData[i]]
                        # get the single marker test statistics
                        # need to retain the V from each test 
                        U,VarU = self.regionsingle(dat, markerset[i])
                        vlist.append(VarU)
                        ulist.append(U)
                    idx = [i for i in range(len(vlist)) if vlist[i] != 0]
                    ulistx = [ulist[i] for i in idx]  # should put these in computer.py
                    vlistx = [vlist[i] for i in idx]
                    datter = [thisData[i] for i in idx]
                    markernames = [markerset[i] for i in idx]
                    self.worker.updateRegionM(len(vlistx))
                    regionstats = self.worker.computeRegionStats(idx, markernames, 
                                                                 ulistx, vlistx,
                                                                 datter, self.familyList)
                    regionstring = map(str, regionstats)
                    chrm = ",".join(list(set(chrms)))
                    print("\t".join([self.analysis, chrm, regionname] + regionstring))
            except KeyboardInterrupt:
                sys.exit(0)
            except IOError:
                print("file io error")
                sys.exit(0)
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise
                #pass

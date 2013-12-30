# fileloader

import family
import sys
import numpy as np
import gzip

class Loader:

    def __init__(self, tfam_file, tped_file):
        """ load the data here """
        # ** THE RARE TEST REQUIRES UNZIPPED DATA! **
        if '.gz' in tped_file or '.gz' in tfam_file:
            sys.stderr.write("Using Zipped Data... affects region based tests.\n")
        self.tpedfile=tped_file
        self.tfamfile=tfam_file
        if '.gz' in tped_file:
            self.tped = gzip.open(self.tpedfile,'r')
        else:
            self.tped = open(self.tpedfile, 'r') # ready to read
        self.famOrder = []
        self.idOrder = []

    def buildFamilyList(self):
        tfam = open(self.tfamfile,'r').read().strip().split("\n")
        fams = dict()
        for f in tfam:
            fs = f.split("\t")
            if len(fs) != 6:
                print("Error: tfam has wrong number of columns")
                sys.exit(1)
            famid = fs[0]
            self.famOrder.append(fs[0]) # the order of markers in the tped file
            self.idOrder.append(fs[1])
            if famid not in fams:
                # then we need a new fam obj
                # order of fs: famid, id, motherID, fatherID, gender, phenotype)
                fams[famid] = family.Family(famid, fs[1], fs[2], fs[3], fs[4], fs[5])
            else:
                # then we need to add to a fam obj
                thisFam = fams[famid]
                thisFam.addMember(famid, fs[1], fs[2], fs[3], fs[4], fs[5])
        return(fams)

    def getMarker(self):
        varline = self.tped.readline()
        # if we still have variants ..
        if len(varline) > 0:
            gs = varline.strip().split("\t")
            marker = gs[1]
            chrm = gs[0]
            alleles = gs[4:]
            return([marker, chrm, alleles])
        else:
            return('')

    def updateFamilyList(self, markerData, familyList):
        # for each pair of genotypes .. update the family obj
        #print(len(self.idOrder))
        #print(len(familyList))
        #print(len(markerData[2])) # already just the alleles
        for i in range(len(self.idOrder)):
            genoidx = 2*i
            familyList[self.famOrder[i]].updateGenotype(markerData, genoidx, self.idOrder[i])

                                                   
#####################################################################################33

    def loadIndex (self, index_file):
        # this file indexes the current tped #
        return(np.load(self.tpedfile+".npy"))
        
    def writeIndex(self):
        # if there's no index file then we should write one
        sys.stderr.write("#Creating a data index ... ")
        finped = open(self.tpedfile, 'r')
        ped_line_offset = []
        ped_line_name = []
        offset = 0
        for line in finped:
            txt = line[0:128].split("\t")
            #ped_line_name.append("chr"+txt[0]+":"+txt[3])
            ped_line_name.append(txt[1])
            ped_line_offset.append(offset)
            offset += len(line)
        finped.close()
        idx = np.array([ped_line_name, ped_line_offset])
        np.save(self.tpedfile+".npy", idx)
        sys.stderr.write("done.\n")                         
        return(self.tpedfile+".npy")
        
    def getMarkersFromLine(self, gs):
        # taking a line first..
        if len(gs) > 0:
            marker = gs[1]
            chrm = gs[0]
            alleles = gs[4:]
            return(marker, chrm, alleles)
        else:
            return([])

    def buildDataSet(self, jdx):
        # this takes the jdx index to a tped, and builds a list of those vars
        # the index tells us where to seek #
        thisData = []
        chrms = []
        markers = []
        fin = open(self.tpedfile,'r')
        for j in jdx:
            fin.seek(int(j))
            l = fin.readline().strip().split("\t")
            a, b, c = self.getMarkersFromLine(l)
            markers.append(a)
            chrms.append(b)
            thisData.append(c)
        return(markers, chrms, thisData)


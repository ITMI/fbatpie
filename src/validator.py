
import sys

# validator #


class Validator:

    def __init__(self, verbose):
        # nothing to do now
        self.verbose = verbose

        
    def validateMarker (self, familyList):
        """ look for problems in inheritance - log errors - zero out problems"""
        for k in familyList:             # for each family
            if ('X' in familyList[k].curChr) and (familyList[k].genders['NB'] == '1'):
                #familyList[k].genotypes['NB'] = ['0','0']
                thisFam = familyList[k]
                cg =  thisFam.genotypes['NB'] # child genotypes  ### EXPECTING ONE CHILD ###
                mg = thisFam.genotypes['M']  # parent 1 genotypes
                fg = thisFam.genotypes['F']  # parent 1 genotypes
                if ('0' in mg) and ('0' in fg) and ('0' in cg):
                    if self.verbose != 'silent':
                        sys.stderr.write("warning: family: " + str(k) + 
                                         " has some missing data at marker " +familyList[k].curMarker + ".\n")
                elif not ((cg[0] == mg[0] and cg[1] == mg[0]) or (cg[0] == mg[1] and cg[1] == mg[1])):  
                    if self.verbose != 'silent':
                        sys.stderr.write("warning: family: " + str(k) + " has a mendelian error for marker " +
                                         familyList[k].curMarker + " .. zeroing out. \n")
                        print("        Boy:    " + str(familyList[k].genotypes))
                    familyList[k].genotypes['NB'] = ['0','0']
            else:
                thisFam = familyList[k]
                cg =  thisFam.genotypes['NB'] # child genotypes  ### EXPECTING ONE CHILD ###
                mg = thisFam.genotypes['M']  # parent 1 genotypes
                fg = thisFam.genotypes['F']  # parent 1 genotypes
                if '0' not in mg and '0' not in fg and '0' not in cg:
                    if not(((cg[0] == mg[0] or cg[0] == mg[1]) and 
                            (cg[1] == fg[0] or cg[1] == fg[1])) or
                        ((cg[0] == fg[0] or cg[0] == fg[1]) and   # must be possible for child genotype
                         (cg[1] == mg[0] or cg[1] == mg[1]))):    # from parents
                        # not valid!
                        if self.verbose != 'silent':
                            sys.stderr.write("warning: family: " + str(k) + " has a mendelian error for marker " +
                                             familyList[k].curMarker + " .. zeroing out. \n")
                        familyList[k].genotypes['NB'] = ['0','0']
                else:
                    if self.verbose != 'silent':
                        sys.stderr.write("warning: family: " + str(k) + 
                                         " has some missing data at marker " +familyList[k].curMarker + ".\n")


                        

    def markerCheck(self, gs, familyList):
        # determing the minor allele and making a count of alleles #
        # if only one allele, add in the missing allele too
        markers = list(set(gs))
        markers.sort()
        markerCount = [0]
        # to represent the frequencies post zeroing out things #
        familyalleles = []
        for k in familyList:
            familyalleles.append(familyList[k].genotypes['NB'][0])
            familyalleles.append(familyList[k].genotypes['F'][0])
            familyalleles.append(familyList[k].genotypes['M'][0])
            familyalleles.append(familyList[k].genotypes['NB'][1])
            familyalleles.append(familyList[k].genotypes['F'][1])
            familyalleles.append(familyList[k].genotypes['M'][1])
        gs = familyalleles
        if len(markers) == 1: 
            markerCount = [0,0]
            markers = ['0', markers[0]]
            markerCount[0] = gs.count(markers[0])
            markerCount[1] = gs.count(markers[1])
        # if only two alleles, but no missing, then add in the missing allele
        elif len(markers) < 3 and markers[0] != '0':
            markers = ['0', markers[0], markers[1]]
        if len(markers) == 2:
            # then we have a mono-allelic marker
            markerCount = [0,0]
            markerCount[0] = gs.count(markers[0])
            markerCount[1] = gs.count(markers[1])
        # with a biallelic marker
        if len(markers) > 2:
            markerCount = [0,0,0]
            markerCount[0] = gs.count(markers[0])
            markerCount[1] = gs.count(markers[1])
            markerCount[2] = gs.count(markers[2])
            # then we can change the order to reflect the minor allele in our data
            if markerCount[1] < markerCount[2]:
                z = markerCount[1]
                markerCount[1] = markerCount[2]
                markerCount[2] = z
                markers = ['0', markers[2], markers[1]]
        return([markers, markerCount])


    def countInfoFams(self, markers, familyList):
        caseInfo = 0
        controlInfo = 0
        totalInfo = 0
        altCount = 0.0
        totCount = 0.0
        freq = 0.0
        # One final validity check ...
        #print("---------------------------------------------------------")
        if len(markers) > 2: # don't run on no-variation markers
            alt = markers[2]
            for k in familyList:
                fg = familyList[k].genotypes['F']
                mg = familyList[k].genotypes['M']
                ng = familyList[k].genotypes['NB']
                if '0' not in mg: # checked this .. needs to be there
                    altCount += mg.count(alt) 
                    totCount += 2
                if '0' not in fg:
                    if familyList[k].curChr == 'X':  # this needs to be there too
                        altCount += fg.count(alt)/2.0
                        totCount += 1
                    else:
                        altCount += fg.count(alt)
                        totCount += 2.0
                if familyList[k].curChr == 'X' and familyList[k].genders['NB'] == '1':
                    if (  (mg.count(alt) == 1) and ('0' not in ng) and ('0' not in mg)  ): 
                        if '1' in familyList[k].phenotypes['NB']:
                            caseInfo += 1
                        else:
                            controlInfo += 1
                        totalInfo += 1
                        #print("par: " + k + str(familyList[k].genotypes) + "  " + str(totalInfo))

                elif (  ((fg.count(alt) == 1) or (mg.count(alt) == 1)) 
                      and ('0' not in ng) and ('0' not in mg) and ('0' not in fg)  ): 
                    if '1' in familyList[k].phenotypes['NB']:
                        caseInfo += 1
                    else:
                        controlInfo += 1
                    totalInfo += 1
                    #print("par: " + k + str(familyList[k].genotypes) + "  " + str(totalInfo))
        if totCount > 0:
            freq = altCount/totCount
        return([totalInfo, controlInfo, caseInfo, freq])


    def checkHet(self, gs):
        # need to check that the fathers are not Het.
        # if they are, then the family should be zeroed out.
        idx = 0
        fam = []
        for i in self.sampleIDs:
            # if its a father
            if 'F' in i:
                #check the genotype
                if gs[idx] != gs[idx+1]:
                    fam.append(i.split("-")[1])
            idx+=2
        if len(fam) > 0: # zero out the family
            idx = 0      # the genotype index
            for i in self.sampleIDs: 
                fi = i.split("-")[1]
                if fi in fam:
                    gs[idx] = '0'
                    gs[idx+1] = '0'
                idx += 2
        return(gs)

        

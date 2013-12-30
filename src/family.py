
# the family object #

# the tfam format
#     Family ID
#     Individual ID
#     Paternal ID
#     Maternal ID
#     Sex (1=male; 2=female; other=unknown)
#     Phenotype

#Affection status, by default, should be coded:
#    -9 missing 
#     0 missing
#     1 unaffected
#     2 affected


class Family:

        def __init__(self, famid, thisid, motherID, fatherID, gender, phenotype):
            self.famid = famid
            self.ids = dict()
            self.genders = dict()
            self.phenotypes = dict()
            self.genotypes = dict()
            self.roles = dict()
            self.stat = dict()
            self.curMarker = ''
            self.curChr = ''
            self.hasHet = 0
            self.addMember(famid, thisid, motherID, fatherID, gender, phenotype)

        def famRole(self, motherID, fatherID, gender):
            if motherID == '0' and fatherID == '0':
                # then its a founder
                if gender == '2':
                    # then it's a mother
                    return("mother")
                else:
                    # then it's a father
                    return("father")
            else:
                # it's a child (non-founder)
                return("child")

        def idToRole(self, thisid):
            if 'M' in thisid:
                return 'M'
            elif 'F' in thisid:
                return 'F'
            elif 'NB' in thisid:
                return 'NB'
            else:
                return 'id error'

        def addMember(self, famid, thisid, motherID, fatherID, gender, phenotype):
            role = self.famRole(motherID, fatherID, gender)
            if role == "mother":
                self.ids['M'] = thisid
                self.genders['M'] = '2'
                self.phenotypes['M'] = phenotype
                self.roles[thisid] = 'M'
            elif role == 'father':
                self.ids['F'] = thisid
                self.genders['F'] = '1'
                self.phenotypes['F'] = phenotype
                self.roles[thisid] = 'F'
            else:
                self.ids['NB'] = thisid
                self.genders['NB'] = gender
                self.phenotypes['NB'] = phenotype
                self.roles[thisid] = 'NB'

        def updateGenotype(self, markerData, genoidx, thisid):
            role = self.roles[thisid]
            self.curMarker = markerData[0]
            self.curChr = markerData[1]
            self.genotypes[role] = markerData[2][genoidx:(genoidx+2)]
            if 'X' in self.curChr:
                # if we're working with a marker on the X chromosome
                # and the father is het .. then it's an error and it should be zeroed out
                # if the NB is male and het, then it's an error. zero it out.
                if 'F' in thisid and self.genotypes['F'][0] != self.genotypes['F'][1]:
                    self.genotypes['F'] = ['0','0']
                if 'NB' in thisid and self.genders['NB'] == '1':
                    if self.genotypes['NB'][0] != self.genotypes['NB'][1]:
                        self.genotypes['NB'] = ['0','0']
	    if 'X' not in self.curChr:
		    if role == 'M' or role == 'F':
			    if self.genotypes[role][0] != self.genotypes[role][1]:
				    self.hasHet = 1
	    else:
		    if role == 'M':
			    if self.genotypes[role][0] != self.genotypes[role][1]:
				    self.hasHet = 1

        def updateStat(self, statName, stat):
            self.stat[statName] = stat

        def printfam(self):
            print("--------------------------------")
            print(self.famid)
            print(self.ids)
            print(self.genders)
            print(self.phenotypes)
            print(self.genotypes)
            print(self.roles)
            print(self.stat)

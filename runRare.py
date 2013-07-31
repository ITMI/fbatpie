
#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    tpedfile = ''     # tped file with genotypes .. constructed via vcftools
    tfamfile = ''     # tfam file showing pedigree and disease status
    offset=0.5        # the offset to the test, 0.5 balances the test.
    weights="madsen"  # the weighting scheme name
    freqfile='none'   # the file of allele freqs with two columns markername [tab] freq
    regionfile=''     # lists of transcripts and the markers associated with them
    indexfile = ''
    try:
        opts, args = getopt.getopt(argv,"h",["tped=","tfam=","offset=","weights=","regionfile=","freqfile=", "indexfile="])
    except getopt.GetoptError:
        print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset x --weights [data,kaviar,none] --regionfile --indexfile'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset x --weights [data,kaviar,none] --regionfile r'
            sys.exit()
        elif opt in ("--tped"):
            tpedfile = arg
        elif opt in ("--tfam"):
            tfamfile = arg
        elif opt in ("--offset"):
            offset = arg
        elif opt in ("--regionfile"):
            regionfile = arg
        elif opt in ("--indexfile"):
            indexfile = arg
        elif opt in ("--weights"):
            weights = arg
    if tfamfile == '' or tpedfile == '':
        print 'runSingle.py --tped <tped file> --tfam <tfam file>  --offset x --weights [madsen,kaviar,none] --regionfile r'
        sys.exit()
            
    f = fbat.FBAT()
    f.load(tfamfile, tpedfile)
    f.setOffset(offset)
    f.rare(regionfile,  # region
           indexfile, # index
           freqfile, # freq
           weights) # weighting

if __name__ == "__main__":
   main(sys.argv[1:])
   

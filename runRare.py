
#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    tpedfile = ''   # tped file with genotypes .. constructed via vcftools
    tfamfile = ''   # tfam file showing pedigree and disease status
    offset=0        # the offset to the test, 0.5 balances the test.
    weights="none"  # the weighting scheme name
    freqfile=''     # the file of allele freqs with two columns markername [tab] freq
    regionfile=''   # lists of transcripts and the markers associated with them
    try:
        opts, args = getopt.getopt(argv,"h",["tped=","tfam=","offset=","weights=","regionfile=","freqfile="])
    except getopt.GetoptError:
        print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset x --weights [data,kaviar,none]'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset x --weights [data,kaviar,none]'
            sys.exit()
        elif opt in ("--tped"):
            tpedfile = arg
        elif opt in ("--tfam"):
            tfamfile = arg
        elif opt in ("--offset"):
            offset = arg
    if tfamfile == '' or tpedfile == '':
        print 'runSingle.py --tped <tped file> --tfam <tfam file>  --offset x --weights [data,kaviar,none]'
        sys.exit()
            
    f = fbat.FBAT()
    f.load(tfamfile, tpedfile)
    f.setOffset(offset)
    f.rare(regionfile, freqfile, weights)

if __name__ == "__main__":
   main(sys.argv[1:])
   

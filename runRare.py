#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    analysis = 'test'
    tpedfile = ''
    tfamfile = ''
    offset=0
    afreq=0.0001
    nfams=1
    mode="additive"
    verbose="silent"
    regionfile=''     # lists of transcripts and the markers associated with them
    indexfile = ''

    try:
        opts, args = getopt.getopt(argv,"h",["analysis=","tped=","tfam=","mode=","offset=","afreq=","verbose=","regions=","regindex="])
    except getopt.GetoptError:
        print "ERR: runRare.py --analysis <a> --tped <tp> --tfam <tf> --offset <o> --afreq <cut> --nfams <f> --mode <m> --regions <r>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runRare.py --tped <tped file> --tfam <tfam file> --offset <o> --freq <cutoff> --nfams <f> --mode <mode> --regions <r>'
            sys.exit()
        elif opt in ("--analysis"):
            analysis = arg
        elif opt in ("--tped"):
            tpedfile = arg
        elif opt in ("--tfam"):
            tfamfile = arg
        elif opt in ("--offset"):
            offset = arg
        elif opt in ("--afreq"):
            afreq = arg
        elif opt in ("--nfams"):
            afreq = arg
        elif opt in ("--mode"):
            mode = arg
        elif opt in ("--verbose"):
            verbose = arg
        elif opt in ("--regions"):
            regionfile = arg
        elif opt in ("--regindex"):
            regindex = arg
    if tfamfile == '' or tpedfile == '':
        print "runRare.py --analysis <a> --tped <tp> --tfam <tf> --offset <o> --afreq <cut> --nfams <f> --mode <m> --regions<r>"
        sys.exit()

    f = fbat.FBAT(analysis, tfamfile, tpedfile, mode, offset, afreq, nfams, verbose)
    f.region(regionfile, regindex)
 
if __name__ == "__main__":
   main(sys.argv[1:])
   

#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    analysis = 'test'
    tpedfile = ''
    tfamfile = ''
    offset=0
    afreq=0.01
    nfams = 10
    mode="additive"
    verbose="silent"
    try:
        opts, args = getopt.getopt(argv,"h",["analysis=","tped=","tfam=","mode=",
                                             "offset=","afreq=","nfams=","verbose="])
    except getopt.GetoptError:
        print "runSingle.py --analysis <a> --tped <tp> --tfam <tf> --offset <o> --afreq <cut> --nfams <f> --mode <m> --verbose <v>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset <o> --afreq <cutoff> --nfams <f> --mode <m>'
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
            afreq = float(arg)
        elif opt in ("--mode"):
            mode = arg
        elif opt in ("--verbose"):
            verbose = arg
        elif opt in ("--nfams"):
            nfams = arg
    if tfamfile == '' or tpedfile == '':
        print "runSingle.py --analysis <a> --tped <tp> --tfam <tf> --offset <o> --afreq <cut> --nfams <f> --mode <m> --verbose <v>"
        sys.exit()
            
    f = fbat.FBAT(analysis, tfamfile, tpedfile, mode, offset, afreq, nfams, verbose)
    f.single()

if __name__ == "__main__":
   main(sys.argv[1:])
   

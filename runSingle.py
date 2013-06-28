
#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    tpedfile = ''
    tfamfile = ''
    offset=0
    try:
        opts, args = getopt.getopt(argv,"h",["tped=","tfam=","offset="])
    except getopt.GetoptError:
        print 'runSingle.py --tped <tped file> --tfam <tfam file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runSingle.py --tped <tped file> --tfam <tfam file>'
            sys.exit()
        elif opt in ("--tped"):
            tpedfile = arg
        elif opt in ("--tfam"):
            tfamfile = arg
        elif opt in ("--offset"):
            offset = arg
    if tfamfile == '' or tpedfile == '':
        print 'runSingle.py --tped <tped file> --tfam <tfam file>'
        sys.exit()
            
    f = fbat.FBAT()
    f.load(tfamfile, tpedfile)
    f.setOffset(offset)
    f.single()

if __name__ == "__main__":
   main(sys.argv[1:])
   

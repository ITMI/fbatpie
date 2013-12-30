
# Pipeline for producing FBAT results #

import sys
import gzip
import os
import time

if len(sys.argv) < 6:
    print("runFBATPipeline.py <analysis> <tped> <pheno> <offset> <mode> <regions> <regindex> <dir name>")
    sys.exit(1)

# runtime arguments #
analysis = sys.argv[1]
tpedfile = sys.argv[2]
tfamfile = sys.argv[3]
offset = sys.argv[4]
mode   = sys.argv[5]
regions = sys.argv[6]
regindex = sys.argv[7]
dirname = sys.argv[8]

# the directory and file prefix #
now = time.localtime()
thisID = str("_".join(map(lambda x: str(x), [dirname, "region", now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec])))

testinfo = open(thisID+".info", 'w')
testinfo.write("\n***************************\nWORKING ID: " + thisID + "\n***************************************\n")
testinfo.write(tpedfile+ "\n" +tfamfile+ "\n"+ offset +"\n"+ mode +"\n"+ regions +"\n")

# where the results go #
cmd1 = "mkdir " + thisID
testinfo.write(cmd1 + "\n")
os.system(cmd1)

# annotating the tfam #
#cmd4 = ("/tools/bin/python2.7 scripts/tfam_annotate.py df4.2.samples " +
#                " annotation/df4.annotated.tfam > " +thisID+"/"+thisID+".anno.tfam annotation/clinical_head.txt "+phenofile)
#testinfo.write(cmd4+"\n")
#os.system(cmd4)

# running the test
cmd5 = ("/tools/bin/python2.7 fbatpie2/runRare.py --analysis " + analysis
        + " --tfam " + tfamfile
        + " --tped " +tpedfile+ " --mode " + mode
        + " --offset " + offset + " --regions " + regions + " --regindex " + regindex
        + " > " +thisID+"/"+ "fbat.out")
testinfo.write(cmd5+"\n")
os.system(cmd5)


cmd7 = ("mv " + thisID + ".info " + thisID)
os.system(cmd7)

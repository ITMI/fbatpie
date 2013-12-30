# Pipeline for producing FBAT results #

import sys
import gzip
import os
import time

onisb = "no"

if len(sys.argv) < 7:
    print("runFBATPipeline.py <vcffile> <tped file> <sample file> <pheno file> <offset> <freq> <dirname> <isb,not> ")
    sys.exit(1)

# runtime arguments #
vcffile = sys.argv[1]
tpedfile = sys.argv[2]
samplefile = sys.argv[3]
phenofile = sys.argv[4]
offset = sys.argv[5]
freq = sys.argv[6]
dirname = sys.argv[7]
onisb = sys.argv[8]
mode = sys.argv[9]

if onisb == "isb":
    py = "/tools/bin/python2.7"
else:
    py = "python"


# the directory and file prefix #
now = time.localtime()
thisID = str("_".join(map(lambda x: str(x), [dirname, now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec])))

testinfo = open(thisID+".info", 'w')
testinfo.write("\n***************************\nWORKING ID: " + thisID + "\n***************************************\n")
testinfo.write(tpedfile +"\n"+ samplefile +"\n"+ phenofile +"\n"+ offset +"\n"+ freq  +"\n")

# where the results go #
cmd1 = "mkdir " + thisID
testinfo.write(cmd1 + "\n")
os.system(cmd1)

# annotating the tfam #
cmd4 = (py + " scripts/tfam_annotate.py " +samplefile+
        " annotation/df4.annotated.tfam > " +thisID+"/"+thisID+".anno.tfam annotation/clinical_head.txt "+phenofile)
testinfo.write(cmd4+"\n")
os.system(cmd4)

# running the test
cmd5 = (py + " fbatpie2/runSingle.py --analysis " + thisID
        + " --tfam " +thisID+ "/" +thisID+".anno.tfam "
        + " --tped " +tpedfile+ " --mode " + mode 
        + " --offset " + offset + " --afreq " +freq+ " > " +thisID+"/"+ "fbat.out")
testinfo.write(cmd5+"\n")
os.system(cmd5)

# formatting!
cmd6 = (py + " scripts/result_format.py " + vcffile + "  " + thisID+"/"+ "fbat.out  > " + thisID+"/"+ "fbat.formatted")
testinfo.write(cmd6+"\n")
os.system(cmd6)

# done with testing and analysis
testinfo.close()

cmd7 = ("mv " + thisID + ".info " + thisID)
os.system(cmd7)

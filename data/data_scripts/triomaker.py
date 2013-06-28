
import argparse

usage="Input a ped file, and it outputs a ped with just trios"
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--ins', required=True, help='The ped file.')
parser.add_argument('--out', required=True, help='The output file containing only trios.')

args, unknown = parser.parse_known_args() #prevent parsing errrors from additional golem-related arguments

infile=args.ins
outfile=args.out


fin = open(infile,'r')
fout = open(outfile,'w')

trio = 0
fam = []
newfam = 0
oldfam = -1

text = fin.read().strip().split("\n")
ped = map(lambda x: x.strip().split(" "), text)
famidx = map(lambda x: x[0], ped)
famset = set(famidx)


def findfam(a, famidx, famset):
    idx = []
    i = famidx.index(a)
    while i < len(famidx):
        if famidx[i] == a:
            idx.append(i)
        if famidx[i] != a:
            break
        i+=1
    return(idx)


def definetrios(famidx, famset):
    keep = []
    for i in famset:
        a = findfam(i, famidx, famset)
        print("family: " + str(i) + "  has: " + str(len(a)) + " members: " + str(a))
        if len(a) == 3:
            keep.append(i)
    return(keep)

keepers = definetrios(famidx, famset)

for f in ped:
    if f[0] in keepers:
        fout.write('\t'.join(f) + "\n")

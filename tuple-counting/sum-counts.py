#!/usr/bin/python

import argparse
from collections import Counter
from plrutils import *

usage = '''
Sum up "count" tables,
i.e. with columns
x y count
AT T 134
AT G 123
'''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--infile', '-i', nargs='*')
parser.add_argument('--outfile', '-o', nargs='?', default="-")
args = parser.parse_args()

outfile = fileopt(args.outfile, "w" )

tuplecount = Counter()

uncounted = 0

for inf in args.infile:
    infile = fileopt(inf, "r" )
    for line in infile:
        vals = line.strip().split()
        try:
            count = int(vals[2])
        except:
            # print "Uncounted: " + line
            uncounted += 1
            continue
        tuplecount[ tuple(vals[:2]) ] += count

print str(uncounted) + " uncounted lines.\n"

# write out
outfile.write("\t".join(["reference","derived","count"])+"\n")
for x in tuplecount:
    outfile.write("\t".join(list(x))+"\t"+str(tuplecount[x])+"\n")


outfile.close()
raise SystemExit

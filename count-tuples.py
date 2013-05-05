#!/usr/bin/python

from collections import Counter

import gzip
import sys

tuplecount = Counter()

if len(sys.argv)>1:
    infile = gzip.open(sys.argv[1],"r")
else:
    infile = sys.stdin

header = infile.readline().strip().split("\t")

# read in data
for x in infile:
    y =  x.strip().split("\t")
    tuplecount[ tuple( y[0:3] ) ] += 1

infile.close()

# write out data
outfile = sys.stdout
outfile.write( "\t".join( header[0:3] + ["count"] ) + "\n" )

for x in tuplecount:
    outfile.write( "\t".join(x) + "\t" + str(tuplecount[x]) + "\n" )

outfile.close()

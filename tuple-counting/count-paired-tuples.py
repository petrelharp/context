#!/usr/bin/python

import argparse
from collections import Counter
from plrutils import *
import re

usage = '''
Count paired tuples from a paired fasta or an axt file.

Assume infile is in a paired format as follows: 
  >04-A-M_0029015
  TTCCGATCTACTAGGACACATGAGGGCTGGAAAGCCACGTTTGGTGAGGCTTCAGGACACGGCTGTGTATTACTGTGCGAGACTTCAAGGGAGGCAGCAGCTAATGCCATTTGACTACTGGGGCCAGGGA
  >04-A-M_0029015_ref IGHV3-64*04 IGHJ4*02_F IGHD6-13*01
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGANNNNNNNNNNNGCAGCAGCTNNNNNNNTTTGACTACTGGGGCCAGGGA

Then count paired tuples of the form:
    WWWWWWWWW
    llMMMMrrr
where the  length of the W's is 'winlen', the length of the l's is 'lwin', and the length of the 'r's is 'rwin'.

If --reverse, also count in the other direction, i.e.
    llMMMMrrr
    WWWWWWWWW

The longer window (Ws) matches the second (reference) sequence, and the shorter window (M's) matches the first sequence.

If 'strict' then ignores anything but ACGT (case-sensitive).
'''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--infile', '-i', nargs='?', default="-")
parser.add_argument('--informat', '-f', nargs='?', default="fasta")
parser.add_argument('--outfile', '-o', nargs='?', default="-")
parser.add_argument('--reverse', '-v', action="store_true")
parser.add_argument('--winlen', '-w', nargs=1)
parser.add_argument('--lwin', '-l', nargs=1)
parser.add_argument('--rwin', '-r', nargs=1)
parser.add_argument('--strict', '-s', action="store_true")

# args.infile='../bcells/A-IGHV3-64_04-IGHJ4_02_F-42-aligned_pairs.fasta'
# args.outfile='../bcells/A-IGHV3-64_04-IGHJ4_02_F-42-aligned_pairs.fasta.counts'
# args.informat='fasta'
# args.winlen=3
# args.lwin=1
# args.rwin=1

args = parser.parse_args()
if not args.winlen or not args.lwin or not args.rwin :
    print "Must specify -w, -l, and -r."
    raise SystemExit
winlen = int(args.winlen[0])
lwin = int(args.lwin[0])
rwin = int(args.rwin[0])
midwin = winlen-lwin-rwin
if (midwin<=0):
    raise ValueError
strict = args.strict

outfile = fileopt(args.outfile, "w" )
tuplecount = Counter()

if args.reverse:
    revfile = fileopt("rev."+args.outfile, "w" )
    revcount = Counter()

if args.informat == "fasta":
    infile = PairedFastaFile( args.infile )
elif args.informat == "axt":
    infile = AxtFile( args.infile )
else:
    raise ValueError

while True:
    try:
        infile.next()
        for x in xrange(len(infile.seq1)-winlen+1):
            word1 = infile.seq2[x:(x+winlen)]
            word2 = infile.seq1[(x+lwin):(x+lwin+midwin)]
            if strict and ( re.search("[^ACGT]",word1) or re.search("[^ACGT]",word2) ):
                next
            else:
                tuplecount[ ( word1, word2 ) ] += 1
            if args.reverse:
                word1 = infile.seq1[x:(x+winlen)]
                word2 = infile.seq2[(x+lwin):(x+lwin+midwin)]
                if strict and ( re.search("[^ACGT]",word1) or re.search("[^ACGT]",word2) ):
                    next
                else:
                    revcount[ ( word1, word2 ) ] += 1
    except StopIteration:
        break

# write out
outfile.write("\t".join(["reference","derived","count"])+"\n")
for x in tuplecount:
    # omit counts with missing values
    if not any([ 'N' in y for y in x ]):
        outfile.write("\t".join(list(x))+"\t"+str(tuplecount[x])+"\n")
outfile.close()

if args.reverse:
    revfile.write("\t".join(["reference","derived","count"])+"\n")
    for x in revcount:
        # omit counts with missing values
        if not any([ 'N' in y for y in x ]):
            revfile.write("\t".join(list(x))+"\t"+str(revcount[x])+"\n")
    revfile.close()

infile.close()

raise SystemExit

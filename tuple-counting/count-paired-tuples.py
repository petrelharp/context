#!/usr/bin/python

import argparse
from collections import Counter
from plrutils import *

usage = '''
Count paired tuples from a paired fasta file.

Assume infile is in a paired format as follows: 
  >04-A-M_0029015
  TTCCGATCTACTAGGACACATGAGGGCTGGAAAGCCACGTTTGGTGAGGCTTCAGGACACGGCTGTGTATTACTGTGCGAGACTTCAAGGGAGGCAGCAGCTAATGCCATTTGACTACTGGGGCCAGGGA
  >04-A-M_0029015_ref IGHV3-64*04 IGHJ4*02_F IGHD6-13*01
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGANNNNNNNNNNNGCAGCAGCTNNNNNNNTTTGACTACTGGGGCCAGGGA

Then count paired tuples of the form:
    WWWWWWWWW
    llMMMMrrr
where the  length of the W's is 'winlen', the length of the l's is 'lwin', and the length of the 'r's is 'rwin'.

The longer window (Ws) matches the second (reference) sequence, and the shorter window (M's) matches the first sequence.

'''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--infile', '-i', nargs='?', default="-")
parser.add_argument('--outfile', '-o', nargs='?', default="-")
parser.add_argument('--winlen', '-w', nargs=1)
parser.add_argument('--lwin', '-l', nargs=1)
parser.add_argument('--rwin', '-r', nargs=1)

# args.infile='../bcells/A-IGHV3-64_04-IGHJ4_02_F-42-aligned_pairs.fasta'
# args.outfile='../bcells/A-IGHV3-64_04-IGHJ4_02_F-42-aligned_pairs.fasta.counts'
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

infile = fileopt( args.infile, "r" )
outfile = fileopt(args.outfile, "w" )

tuplecount = Counter()

while True:
    head1 = infile.readline().strip()
    if not head1:
        break
    seq1 = infile.readline().strip()
    head2 = infile.readline().strip()
    seq2 = infile.readline().strip()
    if not head1 == head2[0:len(head1)] :
        print "Uh-oh: " + head1 + " doesn't match " + head2 + "\n"
        raise ValueError
    for x in xrange(len(seq1)-winlen+1):
        tuplecount[ ( seq2[x:(x+winlen)], seq1[(x+lwin):(x+lwin+midwin)] ) ] += 1


# write out
outfile.write("\t".join(["reference","derived","count"])+"\n")
for x in tuplecount:
    # omit counts with missing values
    if not any([ 'N' in y for y in x ]):
        outfile.write("\t".join(list(x))+"\t"+str(tuplecount[x])+"\n")

infile.close()
outfile.close()

raise SystemExit

#!/usr/bin/python

import argparse
from plrutils import *

usage = '''
Write out the subset of a .axt file
given by an input file of start and end positions (two whitespace-separated columns).
Positions are in the "primary" organism, i.e. the one listed first.

e.g. positions [NO HEADER]
---------------
chrX 124656759 124710720
chrX 124839636 124921900
chrX 124972584 125032464
chrX 125085784 125179309
---------------

.axt
---------------
1 chrX 3025337 3025439 chr2 259145367 259145469 + 7008
AAGCAAGGAAGAGAGCACAAAGCTGCCCTTAAGTCATTGAAGACCTGGAGCACCTATGGGCCTCCCAAGAGGTTCTCACAAAACCTACAACCTCCATCCCAAG
AAGCCAGGAAGAGAGCACAAAACTGCCCTTAAGTCATTGCAGACCCGGAGCTCCGGTGAGCCTGCCACGAGGTTCTCGCGAGATCCGCCGCCTCCATCCCAAG

2 chrX 3025441 3025541 chr2 259146733 259146851 + 6109
GCCTGTGTCTGACAGAGCTGGTATGAGAAAAGGAGACAATTCTTCCCC-----------------AGGTTTTAGCCAGA-CTTTGAGTACACAATGGGATAGTGAGGAGCCCAGTGAAA
GCCTGTGTCTGGCAGAGCCGGTGTGAGAAGAAGAGAGAATCCTTCCCCGCCGCCGGGATAATCAAGAGTTTTGGCCGGACCTTTGAGTACACATCGGGATAGTGAGAAGCCCGGCGAAA
-------------------

Assumes both are in sorted order.
'''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--posfile', '-p', nargs='?', default="-")
parser.add_argument('--axtfile', '-a', nargs='?', default="-")
parser.add_argument('--outfile', '-o', nargs='?', default="-")

args = parser.parse_args()

# example
# tail -n +2 10_filter_regions.out | cut -f 1 | sed -e 's/\(.*\):\([0-9]*\)-\([0-9]*\).*/\1 \2 \3/' > 11_filter_regions.positions.out
args = parser.parse_args("-p test.out -a test.axt -o -".split())

posfile = fileopt( args.posfile, "r" )
axtfile = fileopt( args.axtfile, "r" )
outfile = fileopt( args.outfile, "w" )

posfile.close()
axtfile.close()
posfile = fileopt( args.posfile, "r" )
axtfile = fileopt( args.axtfile, "r" )

for posline in posfile:
    # current position
    posline = posline.strip().split()
    chrom = posline[0]
    pos = [ int(x) for x in posline[1:3] ]
    print "pos:"
    print chrom
    print pos
    # we've moved past that position
    while True :
        while True:
            axtinfo = axtfile.readline().strip().split()
            if not ( len(axtinfo)==0 or axtinfo[0][0] == "#" ) :
                break
        axtchrom = axtinfo[1]
        axtpos = [ int(x) for x in axtinfo[2:4] ]
        seq1 = axtfile.readline().strip()
        seq2 = axtfile.readline().strip()
        print "axt:"
        print axtchrom
        print axtpos
        # next alignment position
        if axtchrom > chrom or axtpos[0] > pos[1] :
            break
        overlap = [ max(pos[0],axtpos[0]), min(pos[1],axtpos[1]) ]
        if overlap[0] < overlap[1] :
            # overlap! output this.
            print "overlap:"
            print overlap
            newinfo = axtinfo[0:2] + overlap + axtinfo[4:5] + [ int(axtinfo[5])+overlap[0]-axtpos[0] + x for x in [0,overlap[1]-overlap[0]] ] + axtinfo[7:]
            outfile.write(" ".join(map(str,newinfo))+"\n")
            outfile.write( seq1[(overlap[0]-axtpos[0]):(overlap[1]-axtpos[1])] + "\n" )
            outfile.write( seq2[(overlap[0]-1):(overlap[1]-1)] + "\n" )


outfile.close()
axtfile.close()
posfile.close()

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
# tail -n +2 10_filter_regions.pos | cut -f 1 | sed -e 's/\(.*\):\([0-9]*\)-\([0-9]*\).*/\1 \2 \3/' > 11_filter_regions.positions.pos
# args = parser.parse_args("-p test.pos -a test.axt -o test.sub.axt".split())
pos = PosFile( args.posfile )
axt = AxtFile( args.axtfile )
outfile = fileopt( args.outfile, "w" )

pos.next()
axt.next()

while True:
    try:
        overlap = [ max(pos.pos[0],axt.pos[0]), min(pos.pos[1],axt.pos[1]) ]
        if pos.chrom == axt.chrom and overlap[0] < overlap[1] :
            # overlap! output this.
            newinfo = axt.info[0:2] + overlap + axt.info[4:5] + [ int(axt.info[5])+overlap[0]-axt.pos[0] + x for x in [0,overlap[1]-overlap[0]] ] + axt.info[7:]
            outfile.write(" ".join(map(str,newinfo))+"\n")
            outfile.write( axt.seq1[(overlap[0]-axt.pos[0]):(overlap[1]-axt.pos[0])] + "\n" )
            outfile.write( axt.seq2[(overlap[0]-axt.pos[0]):(overlap[1]-axt.pos[0])] + "\n" )
            outfile.write("\n")
        if axt.chrom > pos.chrom or axt.pos[1] > pos.pos[1] :
            pos.next()
        if axt.chrom < pos.chrom or axt.pos[1] < pos.pos[1] :
            axt.next()
    except StopIteration:
        break


outfile.close()
pos.file.close()
axt.file.close()

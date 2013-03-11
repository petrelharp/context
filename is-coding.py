#!/usr/bin/python2.7

import argparse
import sys, gzip, signal
from plrutils import *
import collections

parser = argparse.ArgumentParser(description='Write out a collection of files giving the locations of certain sequences.')
parser.add_argument('--infile', '-i', nargs='?', default="-")
parser.add_argument('--codingfile', '-c', nargs='?', help="Tab-separated file with starting & ending positions of coding sequence.")
parser.add_argument('--logfile', '-l', nargs='?', default="-")
parser.add_argument('--outfix', '-o', default="locations", action='store', help='output will be files prefixed by this')

args = parser.parse_args()
# example:
    # args = parser.parse_args("-c CDS-dmel-3R-r5.50.CDS.starts.ends.gz -i 4.rawTAC.gz -o 4.rawTAC -l 4.rawTAC.CDS.log".split())

infile = fileopt(args.infile,"r")
codingfile = fileopt(args.codingfile,"r")
orffile = fileopt(args.outfix + ".orfcoding.gz","w")
nonorffile = fileopt(args.outfix + ".nonorfcoding.gz","w")
noncodingfile = fileopt(args.outfix + ".noncoding.gz","w")

pos = int(infile.readline())
cdstxt = codingfile.readline().split()

while cds:
    cds = map( int, cds )
    while pos < cds[0]:
        noncodingfile.write(str(pos)+"\n")
        pos = int(infile.readline())
    while pos < cds[1]:
        if pos-cds[0]==0:
            orffile.write(str(pos)+"\n")
        else:
            nonorffile.write(str(pos)+"\n")
        pos = int(infile.readline())
    cdstxt = codingfile.readline().split()

# finish off
for pos in infile:
    noncodingfile.write(pos+"\n")

raise SystemExit

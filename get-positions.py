#!/usr/bin/python2.7

import argparse
import sys

parser = argparse.ArgumentParser(description='Write out a collection of files giving the locations of certain sequences.')
parser.add_argument('--in', '-i', type=file, nargs='?', default=sys.stdin)
parser.add_argument('--outfix', '-o', default="locations", action='store', help='output will be files prefixed by this')
parser.add_argument('--patterns', nargs="*", action='store')

args = parser.parse_args()

# keep track of history for this long
maxlen = max( [ len(x) for x in args.patterns ] )

state = [ args.in.read(1) for k in xrange(maxlen) ]

print state

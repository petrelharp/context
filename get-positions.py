#!/usr/bin/python2.7

import argparse
import sys, gzip, signal
from plrutils import *
import collections

parser = argparse.ArgumentParser(description='Write out a collection of files giving the locations of certain sequences.')
parser.add_argument('--infile', '-i', nargs='?', default="-")
parser.add_argument('--logfile', '-l', nargs='?', default="-")
parser.add_argument('--outfix', '-o', default="locations", action='store', help='output will be files prefixed by this')
parser.add_argument('--patterns', nargs="*", action='store')

args = parser.parse_args()
# example:
# args = parser.parse_args("-i 2LHet.raw.gz -o test --patterns GCT GTC ACTAGT".split())

infile = fileopt( args.infile, "r" )
logfile = fileopt( args.logfile, "w" )
outfiles = dict([ (x,fileopt( args.outfix + x + ".gz", "w" )) for x in args.patterns ])
patterns = dict([ (y,[ x for x in y ]) for y in args.patterns ])

# catch ctrl-c gracefully
_exitnow = []

def catch_int(signal,frame):
    _exitnow.append(True)
    if len(_exitnow)>1:
        logfile.write("Caught SIGINT twice, terminating immediately.\n")
        logfile.flush()
        raise SystemExit
    else:
        logfile.write("Caught SIGINT, exiting after this generation.  SIGINT again to terminate.\n")
        logfile.flush()
        pass

signal.signal( signal.SIGINT, catch_int )


# keep track of history for this long
maxlen = max( [ len(x) for x in args.patterns ] )

state = collections.deque([ infile.read(1) for k in xrange(maxlen) ])
pos = 1

while True:
    for patt in args.patterns:
        if all( [ patterns[patt][k] == state[k] for k in xrange(len(patterns[patt])) ] ):
            outfiles[patt].write( str(pos)+"\n" )
    nextchar = infile.read(1)
    if _exitnow or not nextchar:
        break
    state.popleft()
    state.append( nextchar )
    pos += 1

logfile.write("All done with " + str(pos) + " positions.\n")

for outf in outfiles:
    outfiles[outf].close()

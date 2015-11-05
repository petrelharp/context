#!/usr/bin/python2.7

import sys, gzip

def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print "Something not right here."
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj


class PosFile: 
    def __init__(self, filename):
        self.file = fileopt(filename,"r")
        self.chrom = None
        self.pos = [None,None]
    def __iter__(self):
        return self
    def next(self):
        posline = self.file.next().strip().split()
        if len(posline) < 3 :
            print "Bad position line:\n"
            print " ".join(posline)+"\n"
            raise ValueError
        self.chrom = posline[0]
        self.pos = [ int(x) for x in posline[1:3] ]
    def close(self):
        self.file.close()

class AxtFile: 
    def __init__(self, filename):
        self.file = fileopt(filename,"r")
        self.info = None
        self.chrom = None
        self.pos = [None,None]
        self.seq1 = None
        self.seq2 = None
    def __iter__(self):
        return self
    def next(self):
        while True:
            axtline = self.file.next()
            self.info = axtline.strip().split()
            if not ( len(self.info)==0 or self.info[0][0] == "#" ) :
                break
        if len(self.info) < 4:
            print "Bad info line:\n"
            print " ".join(self.info) + "\n"
            raise ValueError
        self.chrom = self.info[1]
        self.pos = [ int(x) for x in self.info[2:4] ]
        self.seq1 = self.file.next().strip()
        self.seq2 = self.file.next().strip()
    def close(self):
        self.file.close()

class PairedFastaFile: 
    def __init__(self, filename, checkfn=lambda s: s):
        self.file = fileopt(filename,"r")
        self.head1 = None
        self.seq1 = None
        self.head2 = None
        self.seq2 = None
        self.checkfn = checkfn
    def __iter__(self):
        return self
    def next(self):
        self.head1 = self.file.next().strip()
        self.seq1 = self.file.next().strip()
        self.head2 = self.file.next().strip()
        self.seq2 = self.file.next().strip()
        if not self.checkfn(self.head1) == self.checkfn(self.head2)[0:len(self.head1)] :
            print "Problem with paired fasta file: " + self.head1 + " doesn't match " + self.head2 + "\n"
            raise ValueError
    def close(self):
        self.file.close()

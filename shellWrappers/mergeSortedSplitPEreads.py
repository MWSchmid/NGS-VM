# -*- coding: utf-8 -*-
"""
usage:
python mergeAndSortSplitPEreads.py R1 R2 outR1 outR2 outUnpairable [--solid]

arguments:
R1 and R1: .fastq(.gz) files with forward/reverse reads (can be compressed, .gz only)
    IMPORTANT: the files must be sorted according to the read name.
outR1 and outR2: .fastq(.gz) files with the merged and sorted forward/reverse reads
outUnpairable: .fastq(.gz) file with the unpairable reads

--solid: if given, names are assumed to be @X_Y_Z_F3 and @X_Y_Z_F3-RNA/BC

requirements:
- four lines per read
- full read names in the first lines (@READNAME)
- f/r info after the name of the read
- space character will be used to separate the read name and the mate info
"""

import gzip
import sys

# check input
try:
    infileNameR1 = sys.argv[1]
    infileNameR2 = sys.argv[2]
    outfileNameR1 = sys.argv[3]
    outfileNameR2 = sys.argv[4]
    outfileNameUP = sys.argv[5]
except:
    print __doc__
    sys.exit(1)


class basicRead(object):
    """A read with the headers and sequences"""
    def __init__(self, seqHead, seq, qualHead, qual):
        self.seqHead = seqHead
        self.seq = seq
        self.qualHead = qualHead
        self.qual = qual
        self.name = self.getName()

    def __str__(self):
        return '\n'.join([self.seqHead, self.seq, self.qualHead, self.qual])

    def getName(self):
        return self.seqHead


class illuminaRead(basicRead):
    """An illumina read with the headers and sequences"""
    def getName(self):
        return self.seqHead.split(' ')[0]


class solidRead(basicRead):
    """A solid read with the headers and sequences"""
    def getName(self):
        return '_'.join(self.seqHead.split('_')[:-1])


def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)


def fastqIter(infileName, readObject):
    """file iterator (returns readObjects - e.g. solidRead or illuminaRead)"""
    with myopen(infileName) as infile:
        while True:
            seqHead = infile.readline().strip()
            if not seqHead:
                break
            seq = infile.readline().strip()
            qualHead = infile.readline().strip()
            qual = infile.readline().strip()
            yield readObject(seqHead, seq, qualHead, qual)

if __name__ == "__main__":
    readObject = solidRead if "--solid" in sys.argv else illuminaRead
    matched = 0
    counter = 1
    with myopen(outfileNameR1, 'w') as outR1, myopen(outfileNameR2, 'w') as outR2, myopen(outfileNameUP) as outUP:
        iterR1 = fastqIter(infileNameR1, readObject)
        iterR2 = fastqIter(infileNameR2, readObject)
        R1 = iterR1.next()
        R2 = iterR2.next()
        while True:
            while R2.name < R1.name:
                print >> outUP, R2
                try:
                    R2 = iterR2.next()
                    counter += 1
                except:
                    break
            if R1.name == R2.name:
                matched += 1
                print >> outR1, R1
                print >> outR2, R2
            else:
                print >> outUP, R1
            try:
                R1 = iterR1.next()
                counter += 1
            except:
                break
            print >> sys.stderr, "processed %d read ends..." % counter
    print >> sys.stderr, "matched %d pairs" % matched

# -*- coding: utf-8 -*-
"""
usage:
python mergeAndSortSplitPEreads.py R1 R2 outR1 outR2 outUnpairable [--solid] [--maxReads <int>]

arguments:
R1 and R1: .fastq(.gz) files with forward/reverse reads (can be compressed, .gz only)
outR1 and outR2: .fastq(.gz) files with the merged and sorted forward/reverse reads
outUnpairable: .fastq(.gz) file with the unpairable reads
--solid: if given, names are assumed to be @X_Y_Z_F3 and @X_Y_Z_F3-RNA/BC
--maxReads <int>: number of reads stored in memory (default is 1 million)

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

def readFileBlock(fastqIterator, numberOfReads):
    """read <numberOfReads> reads into a dictionary with the names (no mate info) as key"""
    notUnique = 0
    out = {}
    for i in xrange(0, numberOfReads, 1):
        try:
            read = fastqIterator.next()
        except:
            break
        if read.name in out:
            if read.seq != out[read.name].seq:
                notUnique += 1
                print >> sys.stderr, "And sequence is NOT THE SAME:\n", read
        out[read.name] = read
    return(out)

if __name__ == "__main__":
    maxReads = int(sys.argv[sys.argv.index("--maxReads")+1]) if "--maxReads" in sys.argv else int(1000000)
    readObject = solidRead if "--solid" in sys.argv else illuminaRead
    matched = 0
    counter = 0
    paired = set([])
    with myopen(outfileNameR1, 'w') as outR1, myopen(outfileNameR2, 'w') as outR2, myopen(outfileNameUP, 'w') as outUP:
        iterR1 = fastqIter(infileNameR1, readObject)
        while True:
            R1 = readFileBlock(iterR1, maxReads)
            counter += len(R1)
            print >> sys.stderr, "loaded %d R1 reads..." % counter
            if not R1:
                break
            iterR2 = fastqIter(infileNameR2, readObject)
            while True:
                try:
                    read = iterR2.next()
                except:
                    break
                if read.name in R1:
                    matched += 1
                    print >> outR1, R1.pop(read.name)
                    print >> outR2, read
                    paired.add(read.name)
                #else:
                    #print >> outUP, read # only works if all in memory
            for rn, read in R1.items():
                print >> outUP, read
        # write out the unpairable reverse reads
        iterR2 = fastqIter(infileNameR2, readObject)
        while True:
            try:
                read = iterR2.next()
            except:
                break
            if read.name not in paired:
                print >> outUP, read
        
    print >> sys.stderr, "matched %d pairs" % matched










#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python runsTestChXPseq.py bwTest bwControl

install dependencies:
sudo easy_install ngslib
sudo pip install skidmarks
btw - this is for deeptools:
sudo pip install pyBigWig

arguments:
\t-bwTest: bigWig of the test sample
\t-bwControl: bigWig of the control sample

optional arguments:
\t--baseSteps: take only every Xth base (default: 100)
\t--fPositive: minimal fraction of positive differences within a candidate region (default: 0.7)
\t--numReps: number of random sets (default: 5)
\t--randomized: do the same as always, but randomize the test set as well (default: off)

notes:
\t-the regions should not be overlapping anymore
\t-the bigWigs need to be normalized:
\t\thttps://github.com/fidelram/deepTools/wiki/Normalizations

idea taken from:
http://www.nature.com/nature/journal/v453/n7197/full/nature06947.html#online-methods
"""

import sys

try:
    bwTest = sys.argv[1]
    bwControl = sys.argv[2]
except:
    print >> sys.stderr, __doc__
    sys.exit(1)

fPositive = float(sys.argv[sys.argv.index("--fPositive")+1]) if "--fPositive" in sys.argv else float(0.7)
numReps = int(sys.argv[sys.argv.index("--numReps")+1]) if "--numReps" in sys.argv else 5
baseSteps = int(sys.argv[sys.argv.index("--baseSteps")+1]) if "--baseSteps" in sys.argv else 100
randomized = "--randomized" in sys.argv

GLOBAL_LEFT_RIGHT=int(1e5/baseSteps)

import wWigIO
from ngslib import BigWigFile
from skidmarks import wald_wolfowitz, serial_test
import numpy as np
import math
import random

# classes and functions
def intersect(a, b):
    return list(set(a) & set(b))


def sign(x):
    """A sign function - IEEE 754 standard."""
    if x > 0 or (x == 0 and math.atan2(x, -1.) > 0.):
        return 1
    else:
        return -1


def sign01(x):
    """Kind of a sign function, but returns 1 and 0 for the runs tests."""
    if x > 0 or (x == 0 and math.atan2(x, -1.) > 0.):
        return 1
    else:
        return 0

sign01vec = np.vectorize(sign01)

def negZeroPos(x):
    """Like a sign function but returning -1, 0 and 1. It's actually cmp(x, 0)."""
    return cmp(x, 0)

negZeroPosVec = np.vectorize(negZeroPos)

class genomicRegion(object):
    """A genomic region with chrom, start, end, values for the borders and
    the fraction of positive differences. Just for printing."""
    def __init__(self, chrom, start, end, LBval, RBval, posFrac):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.size = self.end - self.start
        self.LBval = LBval
        self.RBval = RBval
        self.posFrac = posFrac

    def __str__(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end), str(self.LBval),
                          str(self.RBval), str(self.posFrac),])


class chromosome(object):
    """A chromosome - does everything necessary."""
    def __init__(self, chrom, start, end, baseStep=1e3, contCov=np.zeros(0), testCov=np.zeros(0)):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.size = self.end - self.start
        self.contCov = contCov if contCov.size else np.zeros(int(self.size/baseStep))
        self.testCov = testCov if testCov.size else np.zeros(int(self.size/baseStep))
        self.baseStep = baseStep

    def __str__(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end)])

    def addRawContCov(self, bwConnection):
        """Add the control coverage given a bigWigConnection."""
        self.contCov = bwConnection.getRawCoverage(self.chrom, self.start, self.end)
        return None

    def addRawTestCov(self, bwConnection):
        """Add the test coverage given a bigWigConnection."""
        self.testCov = bwConnection.getRawCoverage(self.chrom, self.start, self.end)
        return None

    def addSmoothContCov(self, bwConnection, winLen=147, winType="flat"):
        """Add the control coverage given a bigWigConnection."""
        self.contCov = bwConnection.getSmoothCoverage(self.chrom, self.start, self.end,
                                                      winLen, winType)
        return None

    def addSmoothTestCov(self, bwConnection, winLen=147, winType="flat"):
        """Add the test coverage given a bigWigConnection."""
        self.testCov = bwConnection.getSmoothCoverage(self.chrom, self.start, self.end,
                                                      winLen, winType)
        return None

    def addRawContCovLargeMode(self, bwConnection):
        """Add the control coverage given a bigWigConnection."""
        self.contCov = bwConnection.getRawCoverageLargeMode(self.chrom, self.start, self.end,
                                                            self.baseStep)
        return None

    def addRawTestCovLargeMode(self, bwConnection):
        """Add the test coverage given a bigWigConnection."""
        self.testCov = bwConnection.getRawCoverageLargeMode(self.chrom, self.start, self.end,
                                                            self.baseStep)
        return None

    def addSmoothContCovLargeMode(self, bwConnection, winLen=1e5, winType="flat"):
        """Add the control coverage given a bigWigConnection."""
        self.contCov = bwConnection.getSmoothCoverageLargeMode(self.chrom, self.start, self.end,
                                                      winLen, winType, self.baseStep)
        return None

    def addSmoothTestCovLargeMode(self, bwConnection, winLen=1e5, winType="flat"):
        """Add the test coverage given a bigWigConnection."""
        self.testCov = bwConnection.getSmoothCoverageLargeMode(self.chrom, self.start, self.end,
                                                      winLen, winType, self.baseStep)
        return None

    def getBorders(self, floatDiff, threshold, randomized=False):
        """Identify edges in a vector of differences."""
        out = {"MP":[], "BT":[], "SC":[], "SKIP":False}
        diff = sign01vec(floatDiff)
        #diff = negZeroPosVec(floatDiff)
        if randomized:
            np.random.shuffle(diff)
        for midPoint in xrange(GLOBAL_LEFT_RIGHT, diff.size-GLOBAL_LEFT_RIGHT, 1):
            lef = np.sum(diff[(midPoint-GLOBAL_LEFT_RIGHT):midPoint])
            rig = np.sum(diff[midPoint:(midPoint+GLOBAL_LEFT_RIGHT)])
            if abs(lef-rig) < threshold:
                continue
            BT = "LB" if lef < rig else "RB"
            out["SC"].append(abs(lef-rig))
            out["MP"].append(midPoint)
            out["BT"].append(BT)
        numLB = sum([x=="LB" for x in out["BT"]])
        numRB = sum([x=="RB" for x in out["BT"]])
        print >> sys.stderr, "Identified %d left and %d right borders." % (numLB, numRB)
        if (numLB == 0) or (numRB == 0):
            out["SKIP"] = True
        return out

    def getBorderThreshold(self, floatDiff, numRep=10):
        """Shuffle everything <numRep> times and get the distribution for the border values."""
        diff = sign01vec(floatDiff)
        #diff = negZeroPosVec(floatDiff)
        temp = numRep*[0]
        for i in xrange(0, numRep, 1):
            np.random.shuffle(diff)
            ranVals = [abs(np.sum(diff[(midPoint-GLOBAL_LEFT_RIGHT):midPoint]) - np.sum(diff[midPoint:(midPoint+GLOBAL_LEFT_RIGHT)])) for midPoint in xrange(GLOBAL_LEFT_RIGHT, diff.size-GLOBAL_LEFT_RIGHT, 1)]
            temp[i] = np.percentile(np.absolute(ranVals), 99)
            print >> sys.stderr, "Perm", i, "- 99th percentile:", temp[i]
        out = sum(temp)/float(len(temp))
        print >> sys.stderr, "Borderthreshold:", out
        return out

    def getRegions(self, minPosFrac=0.7, numRep=10, randomized=False):
        """Compare the control and the test sample.
        minPosFrac: minimal positive fraction
        numRep: number of random sets for the cutoff
        return: a list of significant regions"""
        realDiff = self.testCov - self.contCov
        threshold = self.getBorderThreshold(realDiff, numRep)
        borders = self.getBorders(realDiff, threshold, randomized)
        if borders["SKIP"]:
            return []
        diff = sign01vec(realDiff)
        #diff = negZeroPosVec(realDiff)
        out = []
        prevBorderIdx = borders["BT"].index("LB")
        isOpen = True
        for i in xrange(prevBorderIdx+1, len(borders["MP"]), 1):
            posFrac = diff[prevBorderIdx:i].mean()
            prevBorderType = borders["BT"][prevBorderIdx]
            curBorderType = borders["BT"][i]
            if prevBorderType == curBorderType:
                if posFrac < minPosFrac:
                    prevBorderIdx = i
                continue
            if curBorderType == "RB" and (posFrac >= minPosFrac):
                curStart = borders["MP"][prevBorderIdx]*self.baseStep
                curEnd = borders["MP"][i]*self.baseStep
                sigReg = genomicRegion(self.chrom, curStart, curEnd,
                                       borders["SC"][prevBorderIdx], borders["SC"][i], posFrac)
                if out and (out[-1].start == curStart):
                    out.pop()
                out.append(sigReg)
            prevBorderIdx = i
        return out


class bigWigConnection(object):
    """A bigWig handler"""
    def __init__(self, path):
        self.path = path

    def getChromSizesWigIO(self):
        wWigIO.open(self.path)
        out = wWigIO.getChromSize(self.path)
        out = dict(zip(out[0], out[1]))
        wWigIO.close(self.path)
        return out

    def getChromSizesNGSLIB(self):
        infile = BigWigFile(self.path)
        out = infile.chromSizes()
        out = dict(zip(out[0], out[1]))
        infile.close()
        return out

    def smoothCoverageSave(self, covArray, winLen=147, winType="flat"):
        """Calculate smoothened coverage, see http://www.scipy.org/Cookbook/SignalSmooth.
        Note that winLen should be an odd integer. winType == flat meand moving average."""
        x = covArray
        if x.size < winLen:
            raise ValueError("Input vector needs to be bigger than window size.")
        if x.ndim != 1:
            raise ValueError("Smooth only accepts 1 dimension arrays.")
        if window not in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
            raise ValueError("Valid winTypes are flat, hanning, hamming, bartlett, and blackman.")
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        if window == "flat":
            w = np.ones(winLen, 'd')
        else:
            w = eval("np." + window + "(winLen)")
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]

    def smoothCoverageFlat(self, covArray, winLen=147):
        """Calculate smoothened coverage using moving average.
        Note that winLen should be an odd integer."""
        if (winLen % 2) == 0:
            winLen += 1
        x = covArray
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        w = np.ones(winLen, 'd')
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]

    def smoothCoverageSpec(self, covArray, winLen=147, winType="blackman"):
        """Calculate smoothened coverage with a specified window type.
        Note that winLen should be an odd integer."""
        if (winLen % 2) == 0:
            winLen += 1
        x = covArray
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        w = eval("np." + window + "(winLen)")
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]

    def getRawCoverage(self, chrom, start, end):
        """Retrieve an array with the genome coverage."""
        out = np.zeros(end-start)
        infile = BigWigFile(self.path)
        wigs = infile.fetch(chrom, start, end)
        for wig in wigs:
            out[wig[0]-start:wig[1]-start] = wig[2]
        infile.close()
        return out

    def getSmoothCoverage(self, chrom, start, end, winLen=147, winType="flat"):
        """Retrieve an array with the smoothened genome coverage.
        CHECK/TODO: one might have to extend the coverage array first"""
        raw = self.getRawCoverage(chrom, start, end)
        if winType == "flat":
            out = self.smoothCoverageFlat(raw, winLen)
        else:
            out = self.smoothCoverageSpec(raw, winLen, winType)
        return out

    def getRawCoverageLargeMode(self, chrom, start, end, largeMode=1e3):
        """Retrieve an array with the genome coverage."""
        out = np.zeros(end-start)
        infile = BigWigFile(self.path)
        wigs = infile.fetch(chrom, start, end)
        for wig in wigs:
            out[wig[0]-start:wig[1]-start] = wig[2]
        infile.close()
        return out[::largeMode]

    def getSmoothCoverageLargeMode(self, chrom, start, end, winLen=1e5, winType="flat", largeMode=1e3):
        """Retrieve an array with the smoothened genome coverage.
        CHECK/TODO: one might have to extend the coverage array first"""
        winLen = int(winLen/largeMode)
        raw = self.getRawCoverageLargeMode(chrom, start, end, largeMode)
        if winType == "flat":
            out = self.smoothCoverageFlat(raw, winLen)
        else:
            out = self.smoothCoverageSpec(raw, winLen, winType)
        return out


if (__name__ == "__main__"):
    contBW = bigWigConnection(bwControl)
    testBW = bigWigConnection(bwTest)
    contCS = contBW.getChromSizesNGSLIB()
    testCS = testBW.getChromSizesNGSLIB()
    if len(contCS) != len(testCS):
        print >> sys.stderr, "WARNING: The two files have different numbers of chromosomes."
    chromsToCheck = intersect(contCS.keys(), testCS.keys())
    for chrom in chromsToCheck:
        contSize = contCS[chrom]
        testSize = testCS[chrom]
        size = max([contSize, testSize])
        if size < 10*GLOBAL_LEFT_RIGHT*baseSteps:
            print >> sys.stderr, "skipping chromosome %s due to its size" % chrom
            continue
        print >> sys.stderr, "processing chromosome %s of length %d" % (chrom, size)
        curChrom = chromosome(chrom, 0, size, baseSteps)
        curChrom.addRawContCovLargeMode(contBW)
        curChrom.addRawTestCovLargeMode(testBW)
        sigRegs = curChrom.getRegions(fPositive, numReps, randomized)
        for reg in sigRegs:
            print reg

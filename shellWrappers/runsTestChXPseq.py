#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python runsTestChXPseq.py bwTest bwControl winLen winType fragSize minWinSize [--subWins 5] [--baseSteps 25] [--pCut 0.00001] [--diffCut 1]

install dependencies:
sudo easy_install ngslib
sudo pip install skidmarks
btw - this is for deeptools:
sudo pip install pyBigWig

arguments:
\t-bwTest: bigWig of the test sample
\t-bwControl: bigWig of the control sample
\t-winLen: size of the smoothing window (e.g. 147)
\t\tnote: must be an odd integer
\t-winType: type of the smoothing window (e.g. flat)
\t\tvalues: flat, hanning, hamming, blackman, bartlett
\t\tnote: flat is a moving average
\t-fragSize: the original fragment size (e.g. 1000000)
\t-minWinSize: the minimal window size (e.g. 1000)

optional arguments:
--subWins (5): number of segments a windows is split into
--baseSteps (25): take only every Xth base
--pCut (0.00001): cutoff for the P-value
--diffCut (1): cutoff for the LFC

notes:
\t-the bigWigs need to be normalized:
\t\thttps://github.com/fidelram/deepTools/wiki/Normalizations
"""

import sys

try:
    bwTest = sys.argv[1]
    bwControl = sys.argv[2]
    winLen = int(sys.argv[3])
    winType = sys.argv[4]
    fragSize = int(sys.argv[5])
    minWinSize = int(sys.argv[6])
#    posSetFile = sys.argv[7]
#    negSetFile = sys.argv[8]
except:
    print >> sys.stderr, __doc__
    sys.exit(1)

COMPsubWins = int(sys.argv[sys.argv.index("--subWins")+1]) if "--subWins" in sys.argv else int(5)
COMPbaseSteps = int(sys.argv[sys.argv.index("--baseSteps")+1]) if "--baseSteps" in sys.argv else int(25)
COMPpCutoff = float(sys.argv[sys.argv.index("--pCut")+1]) if "--pCut" in sys.argv else float(.00001)
COMPdiffCutoff = float(sys.argv[sys.argv.index("--diffCut")+1]) if "--diffCut" in sys.argv else int(1)

import wWigIO
from ngslib import BigWigFile
from skidmarks import wald_wolfowitz, serial_test
import numpy as np
import math

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
    """Kind of a sign function, but returns 1 and 0 for the tests."""
    if x > 0 or (x == 0 and math.atan2(x, -1.) > 0.):
        return 1
    else:
        return 0

sign01vec = np.vectorize(sign01)


class genomicRegion(object):
    """A genomic region with chrom, start, end, and a coverage array."""
    def __init__(self, chrom, start, end, contCov=np.zeros(0), testCov=np.zeros(0)):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.size = self.end - self.start
        self.contCov = contCov if contCov.size else np.zeros(self.end-self.start)
        self.testCov = testCov if testCov.size else np.zeros(self.end-self.start)
        self.p = 1
        self.aveDiff = 0

    def __str__(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end), str(self.p),
                          str(self.aveDiff)])

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

    def compare(self, subWins=5, baseSteps=25, pCutoff=.00001, diffCutoff=1, minSize=1e3):
        """Compare the control and the test sample.
        subWins: number of windows to split a region into
        baseSteps: to avoid some of the autocorrelation given by the reads, take only every
        baseSteps-th base
        pCutoff: the p-value cutoff to either follow up or not
        minSize: minimal size of a region to still be tested
        return: a list of significant regions"""
        sigRegs = []
        if self.size < minSize:
            return sigRegs
        realDiff = self.testCov[::baseSteps] - self.contCov[::baseSteps]
        self.aveDiff = realDiff.mean()
        diff = sign01vec(realDiff)
        try:
            testRes = wald_wolfowitz(diff)
            self.p = testRes['p']
        except ZeroDivisionError:  # this means that there are only zeroes (or ones)
            self.p = 0
            if abs(self.aveDiff) > diffCutoff:
                sigRegs = [self]
            return sigRegs
        if self.p > pCutoff:
            return sigRegs
        nextSize = int(self.size/subWins)
        for startPoint in xrange(self.start, self.end, nextSize):
            relStart = startPoint-self.start
            relEnd = relStart+nextSize
            reg = genomicRegion(self.chrom, startPoint, startPoint+nextSize,
                                self.contCov[relStart:relEnd], self.testCov[relStart:relEnd])
            sigRegs.extend(reg.compare(subWins, baseSteps, pCutoff, diffCutoff, minSize))
        if not sigRegs and abs(self.aveDiff) > diffCutoff:
            sigRegs = [self]
        return sigRegs


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

if __name__ == "__main__":
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
        size = max([contSize, testSize])    # in case of Raffaella and Giody I used only contSize
        print >> sys.stderr, chrom, size
        if size < winLen:
            print >> sys.stderr, "skipping it because of the size"
            continue
        for start in xrange(0, size, fragSize):
            if (start % 1e7) == 0:
                print >> sys.stderr, chrom, start
            for offSet in xrange(0, minWinSize, int(minWinSize/5)):
                print >> sys.stderr, "extracting region"
                reg = genomicRegion(chrom, start+offSet, start+offSet+fragSize)
                print >> sys.stderr, "smoothening coverage"
                reg.addSmoothContCov(contBW, winLen, winType)
                reg.addSmoothTestCov(testBW, winLen, winType)
                print >> sys.stderr, "searching segments"
                sigRegs = reg.compare(COMPsubWins, COMPbaseSteps, COMPpCutoff, COMPdiffCutoff, minWinSize)
                for sigReg in sigRegs:
                    print sigReg




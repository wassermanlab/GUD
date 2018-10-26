from array import array
import copy

# Initialize
_binOffsets = array(
    "i", [512+64+8+1, 64+8+1, 8+1, 1, 0])
_binOffsetsExtended = array(
    "i", [4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0])
_binFirstShift = 17
_binNextShift = 3
_binOffsetOldToExtended = 4681
_BINRANGE_MAXEND_512M = 512*1024*1024
_maxParentBin = 584

class BinRange(object):
    """
    Implements a {BinRange} object.
    """

    def binFromRangeStandard(self, start, end):
        """
        This function is a direct implementation in Python of
        Jim Kent's C code written for the UCSC genome browser
        in the source code module src/lib/binRange.c

        The following description is from the binRange.c module:

        Given a start and end in chromosome coordinates assign it
        a bin. There's a bin for each 128k segment, for each
        1M segment, for each 8M segment, for each 64M segment,
        and for each chromosome (which is assumed to be less than
        512M). For extended implementation the top level is 4Gb but
        note that, since start and end are int's, the practical
        limit is up to 2Gb-1, and thus, only four result bins on
        the second level.

        A range goes into the smallest bin it will fit in.
        """
        
        # Initialize
        startBin = start
        endBin = end - 1
        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in range(len(_binOffsets)):
            if startBin == endBin:
                return _binOffsets[i] + startBin

            startBin = startBin >> _binNextShift
            endBin = endBin >> _binNextShift

        raise ValueError(
            "start {0}, end {1} out of range (max is 512 Mb)".format(
            start, end))

    def binFromRangeExtended(self, start, end):

        # Initialize
        startBin = start
        endBin = end-1
        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in range(len(_binOffsetsExtended)):
            if startBin == endBin:
                return (_binOffsetOldToExtended + _binOffsetsExtended[i]
                    + startBin)

            startBin = startBin >> _binNextShift
            endBin = endBin >> _binNextShift

        raise ValueError(
            "start {0}, end {1} out of range (max is 2 Gb)".format(start, end))

    def allBinsInRangeStandard(self, start, end):
        """
        Return a list of all bins which fall into range.
        """

        # Initialize
        bins = []
        startBin = start
        endBin = end-1
        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in range(len(_binOffsets)):
            for b in range(startBin, endBin + 1):
                bins.append(_binOffsets[i] + b)

            if startBin == endBin:
                break

            startBin = startBin >> _binNextShift
            endBin = endBin >> _binNextShift

        return bins

    def allBinsInRangeExtended(self, start, end):
        """
        Return a list of all bins which fall into extended range.
        """

        # Initialize
        bins = []
        startBin = start
        endBin = end-1
        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in range(len(_binOffsetsExtended)):
            for b in range(startBin, endBin + 1):
                bins.append(
                    _binOffsetOldToExtended + _binOffsetsExtended[i] + b
                )

            if startBin == endBin:
                break

            startBin = startBin >> _binNextShift
            endBin = endBin >> _binNextShift

        return bins

    def binFromRange(self, start, end):
        if end <= _BINRANGE_MAXEND_512M:
            return self.binFromRangeStandard(start, end);
        else:
            return self.binFromRangeExtended(start, end);

    def allBinsInRange(self, start, end):
        if end <= _BINRANGE_MAXEND_512M:
            return self.allBinsInRangeStandard(start, end);
        else:
            return self.allBinsInRangeExtended(start, end);

    def childBins(self, pbin):
        """
        Return all immediate child bins for the given bin.
        """

        child_bins = []

        if pbin > _maxParentBin:
            return child_bins

        bin = pbin << _binNextShift

        for i in range(1, 9):
            child_bins.append(bin + i)

        return child_bins

    def parentBin(self, cbin):
        """
        Return all immediate parent bins for the given bin.
        """

        if cbin > 0:
            return cbin >> _binNextShift

        return None

    def descendentBins(self, pbin):
        """
        Return all descendent bins for the given bin.
        """

        descendent_bins = []

        if pbin > _maxParentBin:
            return descendent_bins

        for bin in self.childBins(pbin):
            descendent_bins.append(bin)

            for bin2 in self.descendentBins(bin):
                descendent_bins.append(bin2)

        return descendent_bins

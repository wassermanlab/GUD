#!/usr/bin/env python

""" Class which provides funtionality to facilitate the bin indexing system
as used within the databases underlying the UCSC genome browser.

The binFromRange function is essentially a direct Python implementation of
Jim Kent's C code written for the UCSC genome browser in src/lib/binRange.c

See also the wiki page here for details on the UCSC bin indexing system:
http://genomewiki.ucsc.edu/index.php/Bin_indexing_system

Additional, two functions, childBins and descendentBins, were created as it
seems these are necessary to retrieve all features which fall into a given
chromosome range which does not correspond to a leaf node in the "tree" of
bins.

"""

import argparse
import copy
from array import *


# Initial implementation
_binOffsets = array('i', [512+64+8+1, 64+8+1, 8+1, 1, 0])

# Extended implementation
_binOffsetsExtended = array(
    'i', [4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0])

_binFirstShift = 17
_binNextShift = 3

_binOffsetOldToExtended = 4681

_BINRANGE_MAXEND_512M = 512*1024*1024

_maxParentBin = 584


def range_intersection(start1, end1, start2, end2):
    s = max(start1, start2)
    e = min(end1, end2)

    return e - s


class BinRange(object):
    def binFromRangeStandard(self, start, end):
        """ This function is essentially a direct Python implementation of
        Jim Kent's C code written for the UCSC genome browser in the source
        code module src/lib/binRange.c

        The following description is given verbatim from the binRange.c
        module:

        Given start,end in chromosome coordinates assign it
        a bin. There's a bin for each 128k segment, for each
        1M segment, for each 8M segment, for each 64M segment,
        and for each chromosome (which is assumed to be less than
        512M.) For extended implementation the top level is 4Gb but
        Note, since start and end are int's, the practical limit
        is up to 2Gb-1, and thus, only four result bins on the second
        level.
        A range goes into the smallest bin it will fit in.
        """

        startBin = start
        endBin = end-1

        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in xrange(len(_binOffsets)):
            if startBin == endBin:
                return _binOffsets[i] + startBin

            startBin = startBin >> _binNextShift
            endBin   = endBin >> _binNextShift

        raise ValueError(
            "start {0}, end {1} out of range (max is 512 Mb)".format(
            start, end))


    def binFromRangeExtended(self, start, end):
        startBin = start
        endBin = end-1

        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        for i in xrange(len(_binOffsetsExtended)):
            if startBin == endBin:
                return (_binOffsetOldToExtended + _binOffsetsExtended[i]
                    + startBin)

            startBin = startBin >> _binNextShift
            endBin   = endBin >> _binNextShift

        raise ValueError(
            "start {0}, end {1} out of range (max is 2 Gb)".format(start, end))


    def allBinsInRangeStandard(self, start, end):
        """ Return a list of all bins which fall into range
        """

        startBin = start
        endBin = end-1

        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        bins = []
        for i in xrange(len(_binOffsets)):
            for b in range(startBin, endBin + 1):
                bins.append(_binOffsets[i] + b)

            if startBin == endBin:
                break

            startBin = startBin >> _binNextShift
            endBin   = endBin >> _binNextShift

        return bins


    def allBinsInRangeExtended(self, start, end):
        """ Return a list of all bins which fall into extended range
        """

        startBin = start
        endBin = end-1

        startBin = startBin >> _binFirstShift
        endBin = endBin >> _binFirstShift

        bins = []
        for i in xrange(len(_binOffsetsExtended)):
            for b in range(startBin, endBin + 1):
                bins.append(
                    _binOffsetOldToExtended + _binOffsetsExtended[i] + b
                )

            if startBin == endBin:
                break

            startBin = startBin >> _binNextShift
            endBin   = endBin >> _binNextShift

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
        ''' Return all immediate child bins for the given bin
        '''

        child_bins = []

        if pbin > _maxParentBin:
            return child_bins

        bin = pbin << _binNextShift

        for i in range(1, 9):
            child_bins.append(bin + i)

        return child_bins


    def parentBin(self, cbin):
        ''' Return parent bin of the given bin
        '''

        if cbin > 0:
            return cbin >> _binNextShift

        return None


    def descendentBins(self, pbin):
        ''' Return all descendent bins for the given bin
        '''

        descendent_bins = []

        if pbin > _maxParentBin:
            return descendent_bins

        for bin in self.childBins(pbin):
            descendent_bins.append(bin)

            for bin2 in self.descendentBins(bin):
                descendent_bins.append(bin2)

        return descendent_bins


#
# Started translating the following classes from Jim Kent's code, but this
# has not been completed. It is not clear these are even needed.
#
class binElement(object):
    def __init__(self, start=0, end=0, val=None):
        self.start = start
        self.end = end
        self.val = val


class binKeeper(object):
    def __init__(self, minPos = 0, maxPos = 0):
        if minPos < 0 or maxPos < 0 or minPos > maxPos:
            raise ValueError(
                "Bad range {0}, {1} in initializing binKeeper".format(
                minPos, maxPos))


        self.minPos = minPos;
        self.maxPos = maxPos;
        self.binLists = []

        self.binCount = self.binFromRangeBinKeeperExtended(maxPos-1, maxPos) + 1


    def binKeeperAdd(self, start, end, val):
        if start < self.minPos or end > self.maxPos or start > end:
            raise ValueError("({} {}) out of range ({} {}) in binKeeperAdd".format(start, end, self.minPos, self.maxPos))

        bin = self.binFromRangeBinKeeperExtended(start, end);

        if bin >= self.binCount:
            raise ValueError("bin {} not less than binCount {}".format(bin, self.binCount))

        el = binElement(start, end, val)
        self.binLists[bin].insert(0, el)


    def binKeeperFind(self, start, end):
        """ Return a list of all bins that intersect range.
        """

        if start < self.minPos:
            start = self.minPos
        if end > self.maxPos:
            end = self.maxPos
        if start >= end:
            return None

        startBin = start >> _binFirstShift
        endBin = (end-1) >> _binFirstShift;

        #bin_list = binElement()
        bin_list = []
        for i in xrange(len(_binOffsetsExtended)):
            print "i = {}".format(i)
            #offset = _binOffsetsExtended[i];
            offset = _binOffsets[i];
            for j in xrange(startBin+offset, endBin+offset+1):
                print "j = {}".format(j)
                #for el=bins[j]; el != NULL; el = el->next)
                el_list = self.binLists[j]
                for el in el_list:
                    if rangeIntersection(el.start, el.end, start, end) > 0:
                        newEl = copy.copy(el);
                        bin_list.insert(0, newEl)

            startBin >>= _binNextShift;
            endBin >>= _binNextShift;

        return bin_list;


    def binFromRangeBinKeeperExtended(self, start, end):
        """ This is just like binRange.binFromRangeExtended() above, but it
        doesn't limit the answers to the range from _binOffsetOldToExtended
        and up. It simply uses the whole new bin scheme as if it was the only
        one.
        """

        startBin = start
        endBin = end-1

        startBin >>= _binFirstShift;
        endBin >>= _binFirstShift;

        for i in xrange(len(_binOffsetsExtended)):
            if startBin == endBin:
                return _binOffsetsExtended[i] + startBin
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;

        raise ValueError(
            "Start {0}, end {0} out of range (max is 2Gb)".format(start, end))



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Return bin based on input start end coordinates'
    )

    parser.add_argument(
        '-s', '--start', type=int, required=True, help='Start position on chromosome'
    )

    parser.add_argument(
        '-e', '--end', type=int, required=True, help='End position on chromosome'
    )

    parser.add_argument(
        '-c', '--child', action='store_true', help='Print child bins'
    )

    parser.add_argument(
        '-d', '--descendent', action='store_true', help='Print descendent bins'
    )

    parser.add_argument(
        '-bir', '--bin_in_range', action='store_true', help='Print all bins which overlap the range'
    )

    parser.add_argument(
        '-p', '--parent', action='store_true', help='Print parent bin containing this range'
    )

    args = parser.parse_args()

    start = args.start
    end = args.end

    br = BinRange()

    bin_num = br.binFromRange(start, end)

    print "\nThe bin # = {0}\n".format(bin_num)

    if args.child:
        print "\nThe child bins are:"
        child_bins = br.childBins(bin_num)
        for b in child_bins:
            print b

    if args.descendent:
        print "\nThe descendent bins are:"
        descendent_bins = br.descendentBins(bin_num)
        for b in descendent_bins:
            print b

    #bk = binKeeper(0, 300000000)
    #bin_list = bk.binKeeperFind(start, end)
    #print bin_list

    if args.bin_in_range:
        bin_list = br.allBinsInRange(start, end)

        print("\nAll bins in range:")
        print bin_list
        print("\n")

    if args.parent:
        parent_bin = br.parentBin(bin_num)
        print "\nThe parent bin is: {}".format(parent_bin)

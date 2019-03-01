#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

__author__ = 'Nick Dickens'
__copyright__ = 'Copyright 2015, Nicholas J. Dickens'

'''
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#Python script to read and check bed files files using pybedtools interface to bedtools

import argparse
from pybedtools import BedTool
from numpy import median, max, min, array


parser = argparse.ArgumentParser(description='Reads the max from a region in wiggle file')
parser.add_argument('--bed', required=True, help='the name of the bedgraph file')
parser.add_argument('--chr', required=True, help='the chromosome')
parser.add_argument('--start', required=False, help='the start base')
parser.add_argument('--end', required=False, help='the end base')
parser.add_argument('--window', required=False, help='instead of using the supplied region take use a window from the centre base')
args = parser.parse_args()


def simpleStats (bedData):
    values = []

    for item in bedData:
        values.append(float(item.score))

    values = array(values)
   # print (values)
    dataset = {'min':min(values), 'max': max(values),'med':median(values) }

    return dataset


result = None



startPos = 0
endPos = 0

if args.start and args.end:
    startPos = int(args.start)
    endPos = int(args.end)

if args.window:
    window = int(args.window)
    centrePos = startPos + int((endPos-startPos)/2)
    startPos = int((centrePos-int(window/2))+0.5)
    endPos = int((centrePos+int(window/2))+0.5)
    if startPos < 1:
        startPos = 1

allFileData = BedTool(args.bed)
compareRegion = None
locationString = None
if args.start and args.end:
    locationString = str(args.chr) + ' ' + str(startPos) + ' ' + str(endPos)
    compareRegion = BedTool(locationString, from_string=True)
else:
    startPos = 1
    endPos = 32000000 # Leishmania genome size
    locationString = str(args.chr) + ' ' + str(startPos) + ' ' + str(endPos)
    compareRegion = BedTool(locationString, from_string=True)

#return everything from allFileData that overlaps the compare region
intersectRegions = allFileData.intersect(compareRegion)



result = simpleStats( intersectRegions )



print ("%s\t%s\t%s\t%s" % (locationString, result['min'], result['max'], result['med']))



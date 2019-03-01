#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

__author__ = 'Nick Dickens and Samantha Campbell'
__copyright__ = 'Copyright 2015, Nicholas J. Dickens and Samantha J. Campbell'

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

import pysam
import sys
import os.path
import argparse


from numpy import random
from numpy import mean
from numpy import std

#commandline input
parser = argparse.ArgumentParser(description='Generates multiple bam files with simulated origin regions')
parser.add_argument('--infile', required=True, help='the BAM file to use as input')
parser.add_argument('--prefix', required=False, help='the prefix of the output file (outfiles will be prefix_region.bam), defaults to simulated')
parser.add_argument('--window', required=False, help='the window size in bases (tile) for binning the read data, for speed defaults to 100')
parser.add_argument('--bumpsize', required=False, help='the size of an origin bump (rather than just the SSR), defaults to 30kb')



fileE = args.infile

originBumpSize = 30000
probabilityWindow = 100
outPrefix = 'simulated'

if arg.window:
    probabilityWindow = int(arg.window)

if arg.bumpsize:
    originBumpSize = int(arg.bumpsize)

if arg.prefix:
    outPrefix = str(arg.prefix)


# connect to the bam files

try:
    samfile = pysam.Samfile(fileE, "rb")
except IOError as e:
    print ("Error: cannot open the sam file: %d %s" % (e.errno, e.strerror ))
    sys.exit(1)
except:
    print ("Error opening the bam file: %s" % sys.exc_info()[0])
    sys.exit(1)

# check the file is indexed
# TO DO: Call the samtools indexing pysam function if it isn't indexed
if not os.path.exists(fileE + ".bai"):
    sys.exit('ERROR: There is no index for ' + fileE + '!')


def writeAlignments(alignmentList, bamFileOut):
    for alignment in alignmentList:
        bamFileOut.write(alignment)
    return


# do original classification

#TO DO: this was hard coded for the paper, add bedfile reader code and a comparison region
origin = {'chr' : 'LmjF.36', 'start': 1110127, 'end':1116528}
originLength = (origin.get('end') - origin.get('start')) + 1

#get all reads from the origin coordinates
#originAlignmentsIterator = samfile.fetch(origin.get('chr') , origin.get('start')-1, origin.get('end')-1)

def countAlignmentStarts(alignments, regionStart, regionEnd):
    counter = 0
    for alignment in alignments:
        if alignment.reference_start>=regionStart & alignment.reference_start<=regionEnd:
            counter+=1
    return counter



def probabilityReader (chr, start, end, samfile, window):
    theseAlignments = samfile.fetch(chr, start, end)
    totalStartsForRegion = countAlignmentStarts(theseAlignments, start, end)
    if totalStartsForRegion == 0:
        return 0

    probabilityList = []
    for base in range(start,end,window):
        baseAlignments = samfile.fetch(chr, base, base+window)
        thisCount = countAlignmentStarts(baseAlignments,base,base+window)
        probabilityList.append(thisCount/totalStartsForRegion)
    return probabilityList

def reweightProbabilityList (originProbabilities, genomeProbabilities):
    #check both lists are the same length
    if len(originProbabilities) != len(genomeProbabilities):
        print (">>>ERROR: " + str(len(originProbabilities)) + "!=" + str(len(genomeProbabilities)) + "")
        return 0
    newProbabilities = []
    #reweight to background genome
    for i,j in zip(originProbabilities, genomeProbabilities):
        #newProbabilities.append(i * j)
        newProbabilities.append(i + j) # adding the probabilities is a better way to do it

    newTotal = sum(newProbabilities)
    #rescale so the total for the region is always 1
    for x in range(1, len(newProbabilities),1):
        newProbabilities[x] = newProbabilities[x]/newTotal

    return newProbabilities

def createFakeAlignment (alignmentList,probabilityList,fakeRegionStart,window,outfile):
    totalAlignments = len(alignmentList)
    spareAlignment = alignmentList[0]
    print ("Name:" +spareAlignment.query_name)
    for pIndex in range(0,len(probabilityList),1):
        newCoordinate = fakeRegionStart + (pIndex*window)
        howManyToFix = int((probabilityList[pIndex]*totalAlignments)+0.5)
        for i in range(0,howManyToFix,1):
            newAlignment = spareAlignment
            newAlignment.reference_start=newCoordinate
            #used the outfile because I had a problem with the iterator
            #TO DO: fix this
            outfile.write(newAlignment)
    return



originCentre = int(0.5+(origin.get('start')+(origin.get('end')-origin.get('start'))/2))
originStart = originCentre - int(0.5+(originBumpSize/2))
originAlignmentsIterator = samfile.fetch(origin.get('chr') , originStart, originStart+originBumpSize)
originAlignmentsList = []

for alignment in originAlignmentsIterator:
    originAlignmentsList.append(alignment)

#print ("???? OR ???? " + str(len(originAlignmentsList)) + " ?????????")



#list of target (non-origin strand-switch) regions
nonOrigins = [
    {'chr' : 'LmjF.36', 'start': 155190, 'end':156279},
    {'chr' : 'LmjF.36', 'start': 778681, 'end' :779754},
    {'chr' : 'LmjF.36', 'start': 1413340, 'end':1414163},
    {'chr' : 'LmjF.36', 'start': 1607311, 'end':1608348},
    {'chr' : 'LmjF.36', 'start': 2078727, 'end':2082697},
    {'chr' : 'LmjF.36', 'start': 2468381, 'end':2470493}
]


#origin probabilites
originProbabilities = probabilityReader (origin.get('chr'),originStart, originStart+originBumpSize, samfile, probabilityWindow)


# iterate through target regions
for location in nonOrigins:
    print(location)
    #open an output file
    outfileName = 'simulated_' + location.get('chr') + '_' + str(location.get('start')) + '-' + str(location.get('end')) + '.bam'
    outfile = pysam.AlignmentFile(outfileName, "wb", template=samfile)

    thisCentre = int(0.5+(location.get('start')+(location.get('end')-location.get('start'))/2))
    newStart = thisCentre - int(0.5+(originBumpSize/2))
    if newStart < 0:
        newStart = 0
    newEnd = newStart + originBumpSize
    #create a list of new alignments

    realRegionProbabilities = probabilityReader (location.get('chr'), newStart, newEnd, samfile, probabilityWindow)
    probabilityList  = reweightProbabilityList (originProbabilities, realRegionProbabilities)
    createFakeAlignment (originAlignmentsList,probabilityList,newStart, probabilityWindow, outfile)

    for chr in samfile.references:
        if chr != origin.get('chr'):
            currentAlignments = samfile.fetch(chr)
            writeAlignments(currentAlignments, outfile)
            #pass

        else:
            currentAlignments = samfile.fetch(chr, 0, newStart)
            writeAlignments(currentAlignments , outfile )
            #writeAlignments(newAlignmentsList , outfile )


            chromEnd = samfile.lengths[samfile.references.index(chr)]
            if newEnd > chromEnd:
                newEnd = chromEnd
            currentAlignments = samfile.fetch(origin.get('chr'), newEnd, chromEnd)
            writeAlignments(currentAlignments , outfile )
    #close, sort and index the outfile
    outfile.close()
    pysam.sort(outfileName,outfileName+"srt")
    pysam.index(outfileName+"srt.bam")

#this is for single file output
#close it
#outfile.close()
#sort it
#pysam.sort(outfileName,outfileName+"srt")
#index it
#pysam.index(outfileName+"srt.bam")



#!/opt/bin/python

import sys
import re
import os
import fnmatch

#####################################
## Find all count files matching pattern entered by user (there should be one count file
## per sample) and create a matrix with one row for each gene and one column for each sample
#####################################

#####################################
## Usage: /opt/bin/python/ makeCountMatrix.py rootDir patternToSearch outputFileName
## Example: /opt/bin/python makeCountMatrix.py /ifs/res/liang/RNASeq/Proj2983_MassagueJ "*htseq.count*" Proj2983_MassagueJ_htseq.count_allSamples.txt
#####################################

## create a list of full paths to files that match pattern entered
def findFiles(rootDir,pattern):
    """
    create and return a list of full paths to files that match pattern entered
    """
    filepaths = []
    for path, dirs, files in os.walk(os.path.abspath(rootDir)):
        for dir in dirs:
            if fnmatch.filter(files, pattern):
                for file in fnmatch.filter(files, pattern):
                    filepaths.append(os.path.join(path,file))
            else:
                print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, dir)

    return filepaths

def getSampleID(filepath):
    """
    Extract and return sample ID from file path
    """
    ## this works for first example, but need to make sure file paths and names are
    ## standard before leaving this as is
    samp = filepath[filepath.index("-")+1:filepath.rindex("/")]
    return samp

def printMatrix(matrix,allSamps,outFile):
    """
    Print matrix of counts, with one row for each gene and one column
    for each sample
    """

    header = "Gene Symbol\t" + "\t".join(allSamps)
    with open(outFile,'w') as out:
        print>>out,header
        gns = sorted(matrix.keys())
        for gn in gns:
            row = gn
            for samp in allSamps:
                if matrix[gn].has_key(samp):
                    row += "\t" + matrix[gn][samp]
                else:
                    row += "\t" + "0"
            if not len(row.split("\t")) == len(allSamps)+1:
                print>>sys.stderr, "ERROR: data missing for gene ",gn
            else:
                print>>out,row
    return

def makeCountMatrix(args):
    """
    Find files to parse, create one matrix of all counts and print 
    matrix to file
    """

    rootDir,filePattern,outFile = args

    ## store all counts in a nested dictionary (indexed by gene then by sampleID)
    ## in order to account for genes that might be included in one count file 
    ## but not another
    matrix = {}

    ## store a list of all samples
    allSamps = []

    ## find all files matching pattern entered
    files = findFiles(rootDir,filePattern)

    if files:
        print>>sys.stderr, "\nCombining the following files:\n"
        for file in files:
            print>>sys.stderr, file
            if os.stat(file)[6]==0: ## file is empty
                print>>sys.stderr, "WARNING: This file is empty!"
            else:
                samp = getSampleID(file)
                allSamps.append(samp)
     
                with open(file,'r') as fl:
                    for line in fl:
                        gn,count = line.split()
                        if not matrix.has_key(gn):
                            matrix[gn] = {}
                        if matrix[gn].has_key(samp):
                            print>>sys.stderr, "ERROR: gene",gn,"already has sample",samp
                        else:
                            matrix[gn][samp] = count

        printMatrix(matrix,allSamps,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeCountMatrix(sys.argv[1:])

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
## Usage: /opt/bin/python rnaseq_count_matrix.py rootDir suffixToSearch outputFileName
## Example: /opt/bin/python rnaseq_count_matrix.py /ifs/res/liang/RNASeq/Proj2983_MassagueJ .htseq_count Proj2983_MassagueJ_htseq.count_allSamples.txt
#####################################

def usage():
    print "/opt/bin/python rnaseq_count_matrix.py rootDir patternToSearch outputFileName [idConversionFile]"
    return


## create a list of full paths to files that match pattern entered
def findFiles(rootDir,pattern):
    """
    create and return a list of full paths to files that match pattern entered
    """
    filepaths = []
    pattern = '*'+pattern
    for path, dirs, files in os.walk(os.path.abspath(rootDir)):
        if fnmatch.filter(files, pattern):
            for file in fnmatch.filter(files, pattern):
                filepaths.append(os.path.join(path,file))
        else:
            if "Proj" in path.split("/")[-1] and "-" in path.split("/")[-1]:
                print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, path)

    return filepaths

def getSampleID(filepath,pattern):
    """
    Extract and return sample ID from file path
    """
    if filepath[-1] == "/":
       filepath = filepath[:-1]
    file_name = os.path.basename(filepath)

    #pattern = pattern.replace("*","")
    samp = file_name.replace(pattern,"")

    return samp

def printMatrix(matrix,allSamps,id_idx,outFile):
    """
    Print matrix of counts, with one row for each gene and one column
    for each sample
    """

    if id_idx:
        header = "GeneID\tGeneSymbol\t" + "\t".join(allSamps)
    else:
        header = "GeneSymbol\t" + "\t".join(allSamps)
    with open(outFile,'w') as out:
        print>>out,header
        ids = sorted(matrix.keys())
        for id in ids:
            row = id
            for samp in allSamps:
                if matrix[id].has_key(samp):
                    row += "\t" + matrix[id][samp]
                else:
                    row += "\t" + "0"
            if not (len(row.split("\t")) == len(allSamps)+1 or len(row.split("\t")) == len(allSamps)+2):
                print>>sys.stderr, "ERROR: data missing for record with ID ",id
            else:
                print>>out,row
    return

def indexIdConversions(conversionFile):

    id_idx = {}
    with open(conversionFile,'r') as conv:
        for line in conv:
            id,gn = line.strip().split("\t")
            if not id in id_idx:
                id_idx[id] = gn
            else:
                print "WARNING: multiple gene symbols found for ID %s" %id
    return id_idx

def getGeneSymbol(id, id_idx):
    final_id = []
    part = ''
    if ":" in id:
        part = id[id.index(":"):]
        id = id.split(":")[0]
    ids = id.split("+")
    for id in ids:
        if id in id_idx:
            final_id.append(id_idx[id])
        else:
            print>>sys.stderr, "ERROR: NO GENE SYMBOL FOUND FOR ID %s" %id
    gen_sym = "+".join(final_id) + part
    return gen_sym

def makeCountMatrix(args):
    """
    Find files to parse, create one matrix of all counts and print 
    matrix to file
    """

    id_idx = {}

    if len(args) == 4:
        rootDir,filePattern,outFile,conversionFile = args
        id_idx = indexIdConversions(conversionFile)
    elif len(args) == 3:
        rootDir,filePattern,outFile = args
    else:
        usage()
        sys.exit(1)

    ## store all counts in a nested dictionary (indexed by gene then by sampleID)
    ## in order to account for genes that might be included in one count file 
    ## but not another
    matrix = {}

    ## store a list of all samples
    allSamps = []

    ## find all files matching pattern entered
    files = findFiles(rootDir,filePattern)
    fwg = None

    if files:
        print>>sys.stderr, "\nCombining the following files:\n"
        for file in files:
            if id_idx:
                file_w_gns = file + ".geneSymbols"
                fwg = open(file_w_gns,'w')
            print>>sys.stderr, file
            if os.stat(file)[6]==0: ## file is empty
                print>>sys.stderr, "WARNING: This file is empty!"
            else:
                samp = getSampleID(file,filePattern)
                allSamps.append(samp)
     
                with open(file,'r') as fl:
                    for line in fl:
                        line = line.strip()
                        id,count = line.split()
                        gen_sym = ''
                        if id_idx:
                            gen_sym = getGeneSymbol(id, id_idx)
                            id = "\t".join([id,gen_sym])
                            print>>fwg, "\t".join([id,count])
                        if not matrix.has_key(id):
                            matrix[id] = {}
                        if matrix[id].has_key(samp):
                            print>>sys.stderr, "ERROR: gene",id,"already has sample",samp
                        else:
                            matrix[id][samp] = count
            if id_idx:
               fwg.close()
        printMatrix(matrix,allSamps,id_idx,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeCountMatrix(sys.argv[1:])

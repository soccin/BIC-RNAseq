#!/opt/bin/python

import sys
from collections import OrderedDict

def usage():
    print>>sys.stderr,"%prog gtf_file raw_counts_file output_file [gene_class_to_keep]"
    sys.exit(-1)

if not len(sys.argv) in range(4,6):
    usage()

typeToKeep = None
key = {}

gtfFile = sys.argv[1]
countsFile = sys.argv[2] ## raw ensembl counts, with ensembl ids ONLY
outFile = sys.argv[3]

if len(sys.argv) == 5:
    typeToKeep = sys.argv[4] ## gene category to keep; must be typed exactly as it appears in gtf file

## index info from GTF file by ensembl ID
with open(gtfFile,'r') as gtf:
    for line in gtf:
        if not line[0] == "#":
            infoList = line.strip().split("\t")[8].split()
            infoDict = OrderedDict(zip(infoList[::2],infoList[1::2]))
            try:
                gnID = infoDict['gene_id'].replace("\"","").replace(";","").split(".")[0]
                gnNM = infoDict['gene_name'].replace("\"","").replace(";","")
                gnTP = infoDict['gene_type'].replace("\"","").replace(";","")
                if not gnID in key:
                    key[gnID] = [gnNM,gnTP]
                elif not [gnNM,gnTP] == key[gnID]:
                    key[gnID].append(gnTP)
                    #print>>sys.stderr,"WARNING: %s has multiple gene names and/or types:",[gnNM,tpID],"\t",key[gnID]
            except:
                print>>sys.stderr,"WARNING: gene_id, gene_name or gene_type is missing from record:\n\t",line

## modify counts file to include gene symbols and optionally
## filter for gene class
with open(countsFile,'r') as counts:
    header = counts.readline().strip()
    pts = header.split()
    with open(outFile,'w') as out:
        print>>out,"\t".join([pts[0]]+["GeneSymbol"]+pts[1:])
        for line in counts:
            pts = line.strip().split()
            id = pts[0].split(".")[0]
            counts = pts[1:]
            try:
                gnNm = key[id][0]
                if not typeToKeep or typeToKeep in key[id]:
                    print>>out,"\t".join([pts[0]]+[gnNm]+pts[1:])
            except:
                print>>sys.stderr,"WARNING: No info found for gene ID",id



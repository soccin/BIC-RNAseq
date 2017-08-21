import sys
import re
import os
import fnmatch
from collections import OrderedDict

def usage():
    print "/opt/bin/python merge_rseqc_stats.py stat rootDir patternToSearch outputFileName"
    return

class Stat:
    def __init__(self, name, metric, header_patterns, parse_method, sample_orientation):
        self.name = name
        self.header_patterns = header_patterns
        self.parse_method = parse_method
        self.metric = metric
        self.sample_orientation = sample_orientation # ["row"|"column"]

## create a list of full paths to files that match pattern entered
def find_files(rootDir,pattern):
    """
    create and return a list of full paths to files that match pattern entered
    """
    filepaths = []
    for path, dirs, files in os.walk(os.path.abspath(rootDir)):
        if fnmatch.filter(files, pattern):
            for file in fnmatch.filter(files, pattern):
                filepaths.append(os.path.join(path,file))
        else:
            if "Proj" in path.split("/")[-1] and "-" in path.split("/")[-1]:
                print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, path)

    return filepaths

def sample_id(filepath):
    """
    Extract sample ID from file path
    IMPORTANT: this script assumes that all files are under
    the 'intFiles' subdirectory in the output of the variants
    pipeline. Edit this method if that is not the case
    """
    return(filepath[filepath.find('intFiles')+9:].split("/")[0])

def file_empty(file):
    if(os.stat(file)[6] == 0):
        print "WARNING: File %s is empty!" %file
    return(os.stat(file)[6] == 0) 

def parse_insertion_profiles(line):
    metric,insnt,ninsnt = line.split("\t")
    value = float(insnt)/(float(insnt)+float(ninsnt))
    return(metric,value)

def parse_deletion_profiles(line):
    return(line.split("\t"))

def parse_gc_content(line):
    return(line.split("\t"))

def parse_clipping_profiles(line):
    metric,cnt,ncnt = line.split()
    value = float(cnt)/(float(cnt)+float(ncnt))
    return(metric,value)

def parse_bam_stats(line):
    cols = line.split()
    metric = " ".join([x.strip() for x in cols[0:-1]])
    value = cols[-1].strip()
    return(metric,value)   

def parse_duplication_rate(line):
    return(line.split("\t")) 

def parse_read_distribution(line):
    if line[0] == "=":
        return [None,None]
    metric,x,value,y = line.split()
    return(metric,value)

def print_matrix(matrix, all_y, metric_name, outfile):
    header = "\t".join([metric_name]+[str(x) for x in all_y])
    with open(outfile,'w') as out:
        print >> out, header
        for x in sorted(matrix.keys()):
            ln = [str(x)]
            for y in all_y:
                if not y in matrix[x]:
                    matrix[x][y] = 0
                ln.append(str(matrix[x][y])) ## round this?
            print >> out, "\t".join(ln)
    return
    
def merge_stats(args):
    stats = {"insertion_profile" : Stat("insertion_profile","Position",["Position"],parse_insertion_profiles,"column"),
             "deletion_profile" : Stat("deletion_profile","Position",["read_position"],parse_deletion_profiles,"column"),
             "gc_content" : Stat("gc_content","GC%",["GC%"],parse_gc_content,"column"),
             "clipping_profile" : Stat("clipping_profile","Position",["Position"],parse_clipping_profiles,"column"),
             "bam_stats" : Stat("bam_stats","Stat",["Total"],parse_bam_stats,"row"),
             "duplication_rate" : Stat("duplication_rate","Occurrence",["Occurrence"],parse_duplication_rate,"column"),
             "read_distribution" : Stat("read_distribution","Group",["Group"],parse_read_distribution,"row")
             }

    if len(args) == 4:
        stattype, rootdir, pattern, outfile = args
    else:
        usage()
        sys.exit(1)

    if not stattype in stats:
        print>>sys.stderr,"ERROR: unrecognized stat type: %s" %stattype
        sys.exit(1)

    stat = stats[stattype]
    pattern = '*' + pattern
    files = find_files(rootdir, pattern)

    if files:
        print "\nCombining the following files:\n"
        for hp in stat.header_patterns:
            samples = []
            all_y = []
            matrix = OrderedDict()
            of = outfile
            indivReads = False            

            if stat.name in ["insertion_profile","clipping_profile"]:
                matrix["Read-1"] = OrderedDict()
                matrix["Read-2"] = OrderedDict()
                indivReads = True

            for file in files:
                samp = sample_id(file)
                if file_empty(file):
                    continue
                if samp in samples:
                    print >> sys.stderr, "WARNING: Duplicate sample found! Skipping " + file
                    continue
                print file

                with open(file,'r') as fl:        
                    read = None
                    while not hp in fl.readline():
                        continue
                    if indivReads:
                        read = "Read-1"
                    for line in fl:
                        line = line.strip()
                        if line in ["Read-1:","Read-2:"]:
                            read = line[:-1]
                            continue
                        if not (line and len(line.split()) > 1):
                            continue
                        metric,value = stat.parse_method(line)
                        if metric and value:
                            print metric, value
                            try:
                                metric = int(float(metric))
                            except:
                                print "tried to change ",metric," to int"
                            if stat.sample_orientation == "row":
                                x = samp
                                y = metric
                            else:
                                x = metric
                                y = samp

                            if indivReads:
                                if not read:
                                    print "What the hell?"
                                    break
                                if not x in matrix[read]:
                                    matrix[read][x] = OrderedDict()
                                if not y in matrix[read][x]:
                                    matrix[read][x][y] = 0
                                matrix[read][x][y] = value

                            else:
                                if not x in matrix:
                                    matrix[x] = OrderedDict()
                                if not y in matrix[x]:
                                    matrix[x][y] = 0
                                matrix[x][y] = value
                            
                            if not y in all_y:
                                all_y.append(y)
            if indivReads:
                for r in matrix.keys():
                    if matrix[r]:
                        of = outfile + "_" + r + ".txt"
                        print_matrix(matrix[r], all_y, stat.metric, of)
            else:
                print_matrix(matrix, all_y, stat.metric, of)
        print "\n"
    else:
        print >> sys.stderr, "No files matching pattern " + pattern + " found."

if __name__ == '__main__':
    merge_stats(sys.argv[1:]) 

import sys
import os
import glob

## usage:
#    %prog project_id mapping_file metrics_img_dir
#    e.g., 
#    python qc_metrics.py 09928 /ifs/work/byrne/rnaseq/Proj_09928/Proj_09928_sample_mapping.txt /ifs/work/byrne/rnaseq/Proj_09928/results/metrics/images

projID = sys.argv[1]
sMap = sys.argv[2]
imgDir = sys.argv[3]


### get list of expected samples
samples = []
with open(sMap) as fl:
    for line in fl:
        samples.append(line.strip().split("\t")[1])


### check for per-sample pdf files in image directory
fileExts = {'Mismatch profile' : '.mismatch_profile.pdf',
            'Clipping profile R1' : '.clipping_profile.R1.pdf',
            'Clipping profile R2' : '.clipping_profile.R2.pdf',
            'Read quality boxplots' : '.qual.boxplot.pdf',
            'Read quality heatmap' : '.qual.heatmap.pdf',
            'Insertion profile R1' : '.insertion_profile.R1.pdf',
            'Insertion profile R2' : '.insertion_profile.R2.pdf',
            'Junction Saturation' : '.junctionSaturation_plot.pdf',
            'Splice events' : '.splice_events.pdf',
            'Splice junctions' : '.splice_junction.pdf',
            'Duplication rate' : '.DupRate_plot.pdf',
            'GC content': '.GC_plot.pdf',
            'NVC plot': '.NVC_plot.pdf',
            'Inner distance' : '.inner_distance_plot.pdf',
            'Deletion profile' : '.deletion_profile.pdf'
            }

for smp in samples:
    #print >> sys.stdout, "Checking RSEQC files for sample " + smp + " .... ",
    allFound = True
    
    for nm, fe in fileExts.items():
        expFile = os.path.abspath(imgDir) + "/rseqc_" + smp + fe
        if not os.path.exists(expFile) or os.path.getsize(expFile) == 0:
            allFound = False
            print >> sys.stderr, " "
            print >> sys.stderr, "ERROR: Missing " + nm + " for sample " + smp
            print >> sys.stderr, "   [" + expFile + "] does not exist or is empty."
            print >> sys.stderr, " "

    if allFound:
        print >> sys.stdout, "All files for sample " + smp + " are good." 
    #else:
    #    print >> sys.stdout, "WARNING: file(s) missing. See errors for details."


### check for project-wide files in image directory 
projFiles = {'5prime/3prime bias' : '_picard_5prime3prime_bias.pdf',
             'Alignment distribution' : '_picard_alignment_distribution.pdf',
             'Alignment distribution as percentages' : '_picard_alignment_distribution_percentage.pdf',
             'Alignment summary' : '_picard_alignment_summary.pdf',
             'Alignment summary as percentages' : '_picard_alignment_summary_percentage.pdf',
             'Coverage' : '_picard_coverage.pdf',
             'Clipping profiles' : '_rseqc_clipping_profiles.pdf', 
             'Deletion profiles' : '_rseqc_deletion_profiles.pdf',
             'GC content' : '_rseqc_gc_content.pdf',
             'Insertion profiles' : '_rseqc_insertion_profiles.pdf',
             'Read distribution' : '_rseqc_read_distribution.pdf',
             'Read distribution as percentages' : '_rseqc_read_distribution_percentage.pdf'
             }

print >> sys.stdout, ""
allFound = True
for pf, fe in projFiles.items():
    expFile = os.path.abspath(imgDir) + "/" + projID + fe
    if not os.path.exists(expFile) or os.path.getsize(expFile) == 0:
        allFound = False
        print >> sys.stderr, "ERROR: Missing " + pf 
        print >> sys.stderr, "    [" + expFile + "] does not exist or is empty."
if allFound:
    print >> sys.stdout, "All project PDFs good for " + projID
else:
    print >> sys.stdout, "One or more project PDF(s) is missing. See errors for details."

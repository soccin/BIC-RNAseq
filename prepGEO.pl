#!/usr/bin/perl

#use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use Schedule;
use Cluster;

## Set up everything needed for GEO submissions including meta data
## and a directory containing soft links to final delivery directory
## for easy upload of raw and processed data
#
## NOTE: Right now this script only supports standard differential
## gene analysis including alignment (with either STAR or Tophat2),
## counting with HTSeq and differential expression analysis with DESeq. 
## BEWARE if submitting data processed in any other way (e.g. dexseq,
## cufflinks, etc.)
#
## ALSO NOTE: The output of this file is meant to fit (for copy and paste)
## template downloaded here: 
##      http://www.ncbi.nlm.nih.gov/geo/info/examples/seq_template_v2.1.xls
## on 2017-03-09. 
## GEO claims that this format is updated every now and then. Although I have
## not seen it change in several years, be aware that it could be updated 
## without notice and the output from this script may not match the template 
## exactly.

my ($map, $pre, $output, $request, $config, $aligner, $species, $htseq, $deseq, $dexseq, $pass1, $strand);

$pre = 'TEMP';
$output = 'results';
$aligner = 'star';

my $curDir = `pwd`;
chomp $curDir;
my $uID = `/usr/bin/id -u -n`;
chomp $uID;

GetOptions ('map=s' => \$map,
            'pre=s' => \$pre,
            'output=s' => \$output,
            'request=s' => \$request,
            'config=s' => \$config,
            'aligner=s' => \$aligner,
            'species=s' => \$species,
            'dexseq' => \$dexseq,
            'htseq' => \$htseq,
            'deseq' => \$deseq,
            'pass1' => \$pass1, 
            'strand=s' => \$strand) or exit(1);

if(!$map || !$species || !$request || !$config || !$strand){
    print <<HELP;

    USAGE: prepGEO.pl -map MAP -pre PRE -species SPECIES -request REQUEST -config CONFIG -strand STRAND 
        * MAP file listing sample mapping information for processing (REQUIRED)
        * SPECIES: only hg19, mouse (mm10; default) and human-mouse hybrid (hybrid) currently supported (REQUIRED)
        * CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
        * REQUEST: project request file containing PI, Investigator and project IDs (REQUIRED)
        * STRAND: library strand; valid options are none, forward, reverse (REQUIRED)
        * ALIGNER: program used to align reads; valid options are star or tophat (default: star)
        * PRE: output prefix (default: TEMP)
        * OUTPUT: output results directory (default: results)
        * DESEQ: deseq was run; include this in list of data processing steps
        * HTSEQ: htseq was run; include this in list of data processing steps and use normalized counts file as processed data file
        * PASS1: star was run using the 1-pass method
HELP
exit;
}

if($output !~ /^\//){
    $output = "$curDir/$output";
}

chomp $aligner;
chomp $species;

if($species !~ /hg19|mm9|mm10|zv9|zv10|dm3|WBcel235|hybrid/i){
    die "Species must be hg19, mm9, mm10, zv9, zv10, dm3 or WBcel235";
}
if($aligner !~ /tophat|star/i){
    die "Aligner must be tophat or star";
}

my $htseq_stranded = '';
if($strand =~ /none/i){
    $htseq_stranded = 'no';
} elsif($strand =~ /forward/i){
    $htseq_stranded = 'yes';
} elsif($strand =~ /reverse/i){
    $htseq_stranded = 'reverse';
}

my $rawDir = "$output/geo/bicgeo/raw";
my $processedDir = "$output/geo/bicgeo/processed";
my $metaDir = "$output/geo/bicgeo/meta";

if($aligner =~ /tophat/i){
    $rawDir = "$output/geo/bicgeo/$aligner/raw";
    $processedDir = "$output/geo/bicgeo/$aligner/processed";
    $metaDir = "$output/geo/bicgeo/$aligner/meta";
}
my $geoFile = "$metaDir/meta.txt";

print "Creating geo directory...";
`/bin/mkdir -m 755 -p $rawDir`;
`/bin/mkdir -m 755 -p $processedDir`;
`/bin/mkdir -m 755 -p $metaDir`;
print "Done.\n";

### construct final delivery path using request file; we will link 
### to files in this path from geo directory

print "Parsing request file...";
my $pi = '';
my $inv = '';
my $proj = '';
my $run = '';

open(REQUEST, "$request") || die "Can't open request file $request $!";
while(<REQUEST>){
    chomp;
    my @req = split(/: /, $_);
    if($req[0] =~ /^pi$/i){
        $pi = $req[1];
        $pi =~ s/\@mskcc\.org//;
    }
    elsif($req[0] =~ /^investigator$/i){
        $inv = $req[1];
        $inv =~ s/\@mskcc\.org//;
    }
    elsif($req[0] =~ /projectid/i){
        $proj = $req[1];
    }
    elsif($req[0] =~ /runnumber/i){
        $run = sprintf("%03d",$req[1]);
        $run = "r_$run";    
    }
}

if(!$pi || !$inv || !$proj){
    die "Could not construct final delivery path. Need PI, Investigator and ProjectID from request file";
}
if(!$run){
    $run = "r_001";
}

my $delivery_root = "/ifs/res/seq";
my $delivery_path = "$delivery_root/$pi/$inv/$proj/$run";
print "Done.\n";

### get versions of all programs mentioned in list of "data processing steps"
print "Parsing config file...";
my $cutadapt_version = '[UNKNOWN]';
my $star_version = '[UNKNOWN]';
my $htseq_version = '[UNKNOWN]';
my $deseq_version = '[UNKNOWN]';
my $dexseq_version = '[UNKNOWN]';
my $r_version = '[UNKNOWN]';
my $tophat_version = '[UNKNOWN]';
my $cufflinks_version = '[UNKNOWN]';
my $bowtie_version = '[UNKNOWN]';
my $express_version = '[UNKNOWN]';
my $kallisto_version = '[UNKNOWN]';
my $picard_version = '[UNKNOWN]';

my $SAMTOOLS = '';
my $R = '';

open(CONFIG, "$config") || die "Can't open config file $config $!";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    my @path = split(/\//, $conf[1]);

    if($conf[0] =~ /samtools/i){
        $SAMTOOLS = $conf[1];
    }
    if($conf[0] =~ /cutadapt/i){
        $cutadapt_version = $path[5];
    }
    elsif($conf[0] =~ /^star$/i){
        $star_version = $path[5];
    }
    elsif($conf[0] =~ /htseq/i){
        $htseq_version = $path[5];
    }
    elsif($conf[0] =~ /^r$/i){
        $R = $conf[1];
        $r_version = $path[5];
    }
    elsif($conf[0] =~ /tophat/i){
        $tophat_version = $path[5];
    }
    elsif($conf[0] =~ /cufflinks/i){
        $cufflinks_version = $path[5];
    }
    elsif($conf[0] =~ /bowtie/i){
        $bowtie_version = $path[5];
    } 
    elsif($conf[0] =~ /express/i){
        $express_version = $path[5];
    }
    elsif($conf[0] =~ /kallisto_version/i){
        $kallisto_version = $path[5];
    }
    elsif($conf[0] =~ /picard/i){
        $picard_version = $path[5];
    }
    elsif($conf[0] =~ /dexseq/i){
        $dexseq_version = $path[5];
    }
}
close CONFIG;
print "Done.\n";

open(my $GEO, ">", $geoFile) or die "Could not open file '$geoFile' $!";

## get unique sample names & pe/se from mapping file so we 
## 1) have a list of samples for GEO template and
## 2) can construct bam file names and make links in geo directory

print "Getting sample info from mapping file...";    
my %mapping_samples = ();
open(MA, "$map") or die "Can't open mapping file $map $!";
while(<MA>){
    chomp;
    my @data = split(/\s+/, $_);
    my $ended = '';
    if($data[4] eq "PE"){
        $ended = 'paired-end';
    } elsif($data[4] eq "SE"){
        $ended = 'single-end';
    } else {
        die "Do not recognize 'ended-ness' $data[4]";
    }
    $mapping_samples{$data[1]} = $ended;
}         
close MA;
print "Done.\n";

open(my $GEO, ">", $geoFile) or die "Could not open file '$geoFile' $!";

print $GEO "The meta data in this file is meant to be copied\
and pasted into the template Excel file downloaded from
GEO. IMPORTANT: It does NOT contain all info needed for 
a complete and valid submission. This file simply provides
as much information as could be gathered automatically 
during the pipeline run. BE SURE TO REVIEW IT CAREFULLY 
AND COMPLETE IT AS NEEDED BEFORE ATTEMPTING TO SUBMIT TO 
GEO.\n\n";

##    
## write list of sample names & make links to bams
##
print "Creating soft links to final bam locations...";
print $GEO "[SAMPLE INFORMATION]\n";
print $GEO "title\n";

foreach my $ms(keys %mapping_samples){
    print $GEO $ms . "\n";
    #my $del_bam = $delivery_path . "/gene/alignments/" . $proj . "_" . $ms . ".bam";
    my $del_bam = "../../../gene/alignments/" . $proj . "_" . $ms . ".bam";
    my $geo_bam = $rawDir . "/" . $ms . ".bam";
    `cd $rawDir; ln -s $del_bam $geo_bam`;
}
print "Done.\n";

print $GEO "\n\n";

##
## write data processing steps
##
print "Writing data processing steps...";
print $GEO "[DATA PROCESSING PIPELINE]\n";

print $GEO "data processing step\t" . "CutAdapt version " . $cutadapt_version . " was used to trim and filter reads. Adapters were trimmed with minimum overlap of 10 bases and reads of less than half original read length were filtered out.\n";
if($aligner eq "tophat"){
    print $GEO "data processing step\t" . "Tophat version " . $tophat_version . " with Bowtie version " . $bowtie_version . " was used to align reads. Parameters used: -p 6 -r 70 --mate-std-dev 90\n"; ### TO DO: add GTF info????
    print $GEO "data processing step\t" . "Aligned reads were sorted contigs in reference file using PICARD version " . $picard_version . ".\n";
}
if($aligner eq "star"){
    if($pass1){
        print $GEO "data processing step\t" . "STAR version " . $star_version . " was run using the 1-pass method with the following parameters: --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within\n";
    } else {
        print $GEO "data processing step\t" . "STAR version " . $star_version . " was run using the 2-pass method. First pass was run with the following parameters: --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within. Genome database was then generated using the minimum read length of half desired read length and second pass was run against the new genome file with the same parameters as pass 1.\n";
    }
}
if($htseq){
    print $GEO "data processing step\t" . "Expression count matrix was computed using HTSeq version " . $htseq_version . ". Parameters set were: -m intersection-strict -s " . $htseq_stranded . " -t exon\n";
}
if($deseq){
    # get DESeq version (it is not in config file)
    my @deseq_info = `$R/Rscript $Bin/packageInfo.R "pkgs=c('deseq')" "fields=c('version')"`;
    foreach my $line(@deseq_info){
        chomp $line;
        if($line =~ /deseq/i){
            my @info = split(/ {1,}/, $line);
            $deseq_version = $info[1];
        }
    }
    print $GEO "data processing step\t" . "Raw count matrix was processed using the R/Bioconductor package DESeq, version " . $deseq_version . ", which is used to both normalize the full dataset and analyze differential expression between sample groups.\n";
}

##
## write genome build
##
print $GEO "genome build\t$species\n";

##
## write desc of processed data file
##
print $GEO "processed data files format and content\t" . "tab-delimited file containing normalized counts\n";

print "Done.\n";
print $GEO "\n\n";

##
## write processed data file info & link to final processed file location
##
print "Writing processed data file info...";
print $GEO "[PROCESSED DATA FILES]\n";
print $GEO "file name\tfile type\tfile checksum\n";
if($deseq){
    my $current_counts_file = $output . "/gene/counts_gene/counts_scaled_DESeq.xls";
    my $del_counts_file = $delivery_path . "/gene/counts_gene/counts_scaled_DESeq.xls";
    my $geo_counts_file = $processedDir . "/" . "counts_scaled_DESeq.txt";
    my $counts_md5 = `md5sum $current_counts_file | awk '{ print \$1 }'`;
    print $GEO "counts_scaled_DESeq.txt\ttab-delimited\t$counts_md5\n";

    `ln -s $del_counts_file $geo_counts_file`;
}
print "Done.\n";

print "Getting insert sizes...";
## store insert sizes, assuming at least one sample is paired end"
my $isfile = "$output/metrics/$pre\_InsertSizeMetrics.txt";
open(ISFILE, "$isfile") || die "Can't open insert size metrics file $isfile $!";
my $hflag = 0;
my $mean_idx = -1;
my $stddv_idx = -1;
my $samp_idx = -1;
my %ismets = ();
while(<ISFILE>){
    chomp;
    my @metrics = split(/\t/, $_);       
    if($hflag == 0){
        $hflag = 1;
        my $hidx = 0;
        foreach my $met (@metrics){
            if($met eq "MEAN_INSERT_SIZE"){
                $mean_idx = $hidx;
            } elsif($met eq "STANDARD_DEVIATION"){
                $stddv_idx = $hidx;
            } elsif($met eq "SAMPLE"){
                $samp_idx = $hidx;
            }
            $hidx+=1;
        }
    } else {
        my @mets = ($metrics[$mean_idx],$metrics[$stddv_idx]);
        $ismets{$metrics[$samp_idx]} = [@mets];
    }
}
close ISFILE;
print "Done.\n";

##
## write raw file info
##
print "Writing raw file info...";
print $GEO "[RAW FILES]\n";
print $GEO "file name\tfile type\tfile checksum\tinstrument model\tread length\tsingle or paired-end\n";

foreach my $ms(keys %mapping_samples){
    print "\n  $ms";
    my $pe = $mapping_samples{$ms};

    my $geobam = $ms . ".bam"; ## bam file name to be uploaded to geo
    my $curbam = $output . "/gene/alignments/" . $pre . "_" . $ms . ".bam";  ## path to bam in current working directory
    #my $curbam = "../../../gene/alignments/" . $pre . "_" . $ms . ".bam"; ## path to bam in current working directory
    print "\n    Getting md5sum...";
    my $bammd5 = `md5sum $curbam | awk '{ print \$1 }'`;
    chomp $bammd5;
    print "Done.";
    print "\n    Getting read length..."; 
    my $readlen = `$SAMTOOLS/samtools view $curbam | awk '{print length(\$10)}' | head -1000 | sort -u`; 
    chomp $readlen;
    ### TO DO: handle case where more than one value is returned from ^^^^
    print "Done.";
    print $GEO "$geobam\tbam\t$bammd5\t[UNKNOWN]\t$readlen\t$pe\n";
}
print "\nDone.\n"; 

print $GEO "\n\n";

print "Writing paired end info...";
print $GEO "[PAIRED-END EXPERIMENTS]\n";
print $GEO "file name 1\tfile name 2\taverage insert size\tstandard deviation\n";

foreach my $ms(keys %mapping_samples){
    if($mapping_samples{$ms} eq "paired-end"){
        ## get insert size metrics
        my $geobam = $ms . ".bam";
        print $GEO "$geobam\t$geobam\t$ismets{$ms}[0]\t$ismets{$ms}[1]\n";
    }
}
print "Done.\n";
close GEO;

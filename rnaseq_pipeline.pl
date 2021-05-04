#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

### takes a list of directories containing fastq.gz files
### 1) makes directories for each sample
### 2) soft lins to all fastq.gz files in directory

### ex. input list $file
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14028 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14278 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14314 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14348 SE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14512 SE
### etc.

### config file: list of programs
### PROGRAM_NAME PATH_TO_PROGRAM
### ex.
### STAR                    /opt/common/star/STAR_2.3.0e
### PICARD                  /opt/common/picard/picard-tools-1.104
### NOTE: MUST GIVE FULL PATH TO CONFIG FILE

### NOTE2: FOR FUSION, IF A SAMPLE CONSISTS OF PE AND SE READS e.g. HAS MULTIPLE RUNS, AND SOME ARE PE AND OTHERS ARE SE
###                    THIS WILL CAUSE FUSIONS TO NOT WORK BECAUSE OF UNEVEN READ FILES CAUSE DURING CAT OF ALL READS


my ($map, $pre, $config, $request, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1, $lncrna, $lincrna_BROAD, $output, $strand, $r1adaptor, $r2adaptor, $scheduler, $rsem, $kallisto, $express, $standard_gene, $clustering_only, $differential_gene, $standard_transcript, $differential_transcript, $alignment_only);

$pre = 'TEMP';
$output = "results";
my $priority_project = "ngs";
my $priority_group = "Pipeline";
my $rseqc_batch_size = 10;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $svnRev = `svn info $Bin | grep Revision | cut -d " " -f 2`;
chomp $svnRev;

my $email = "$uID\@cbio.mskcc.org";
my $rsync = "/juno/res/bic/$uID";

GetOptions ('map=s' => \$map,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
            'request=s' => \$request,
	    'samplekey=s' => \$samplekey,
	    #'comparisons=s' => \$comparisons,
	    'h|help' => \$help,
	    'star' => \$star,
	    'pass1' => \$pass1,
	    'tophat|tophat2' => \$tophat,
	    'cufflinks' => \$cufflinks,
	    'dexseq' => \$dexseq,
	    'htseq' => \$htseq,
	    'deseq' => \$deseq,
	    'chimerascan' => \$chimerascan,
	    'star_fusion' => \$star_fusion,
	    'mapsplice' => \$mapsplice,
	    'defuse' => \$defuse,
	    'fusioncatcher' => \$fusioncatcher,
	    'allfusions' => \$allfusions,
	    'kallisto' => \$kallisto,
	    'rsem' => \$rsem,
	    'express' => \$express,
            'species=s' => \$species,
            'strand=s' => \$strand,
            'standard_gene' => \$standard_gene,
            'differential_gene|diff_gene' => \$differential_gene,
            'clustering_only' => \$clustering_only,
            'standard_transcript|standard_trans' => \$standard_transcript,
            'differential_transcript|diff_trans' => \$differential_transcript,
            'lncrna' => \$lncrna,
 	    'output|out|o=s' => \$output,
	    'rsync=s' => \$rsync,
	    'email' => \$email,
	    'r1adaptor=s' => \$r1adaptor,
	    'r2adaptor=s' => \$r2adaptor,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
            'lincrna_BROAD' => \$lincrna_BROAD,
            'alignment_only' =>\$alignment_only) or exit(1);


if(!$map || !$species || !$strand || !$config || !$request || $help){
    print <<HELP;

    USAGE: rnaseq_pipeline.pl -map MAP -species SPECIES -strand STRAND -config CONFIG -pre PRE -samplekey SAMPLEKEY -comparisons COMPARISONS -request REQUEST
	* MAP: file listing sample mapping information for processing (REQUIRED)
	* SPECIES: only hg19, mouse (mm10; default) and human-mouse hybrid (hybrid), zebrafish (z11), fly (dm6) currently supported (REQUIRED)
	* STRAND: library strand; valid options are none, forward, reverse (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
        * REQUEST: file containing all request information including PI, Investigator and ProjectID (REQUIRED)
	* SCHEDULER: currently support for SGE, LUNA and JUNO  (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* STANDARD_GENE: standard analysis - star alignment, htseq gene count
	* DIFFERENTIAL_GENE: differential gene analysis - star alignment, htseq gene count, deseq, counts normalization, clustering, and gsa
	* STANDARD_TRANSCRIPT: standard analysis - rsem transcript counts
	* DIFFERENTIAL_TRANSCRIPT: differential transcript analysis - rsem transcript counts, deseq, counts normalization, and clustering
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B (if -deseq, REQUIRED)
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B (if -deseq, REQUIRED)
	* EMAIL: email to send notication of finished final job of pipeline (default: $uID\@cbio.mskcc.org)
	* RSYNC:  path to rsync data for archive (default: /ifs/res/$uID)
	* R1ADAPTOR/R2ADAPTOR: if provided, will trim adaptor sequences; NOTE: if provided for only one end, will also assign it to the other end
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless -pass1 specified; tophat2 (-tophat); if no aligner specifed, will default to STAR
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons); fusion callers chimerascan (-chimerascan), rna star (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs; transcript analysis using express (-express), kallisto (-kallisto) and rsem (-rsem), and all (-transcript)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* OUTPUT: output results directory (default: results)
        * OPTIONS: lncRNA analysis (-lncrna) runs all analyses based on lncRNA GTF (hg19 only); 
HELP
exit;
}

if($ENV{'LSF_ENVDIR'} eq "/common/lsf/conf"){
    $scheduler = 'luna';
}
elsif($ENV{'LSF_ENVDIR'} eq "/common/juno/OS7/conf" or $ENV{'LSF_ENVDIR'} eq "/admin/lsfjuno/lsf/conf"){
    $scheduler = 'juno';
}
elsif($ENV{'SGE_ROOT'} ne ""){
    $scheduler = 'sge';
}
else{
    die "unrecognized scheduler, valid scheduler [sge, luna, juno]";
}


if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

if($output !~ /^\//){
    $output = "$curDir/$output";
}

if($species !~ /human|hg19|mouse|mm9|mm10|hybrid|zebrafish|zv9|zv10|z11|dm6|fly|WBcel235|sacCer3|rat|Rnor6|rn6/i){
    die "Species must be human (hg19), mouse (mm9/mm10), human-mouse hybrid (hybrid), fly (dm6), C.elegans (WBcel235), yeast (sacCer3) or zebrafish (zv9/zv10/z11)";
}

if($r1adaptor){
    if(!$r2adaptor){
	$r2adaptor = $r1adaptor;
    }
}
else{
    if($r2adaptor){
	$r1adaptor = $r2adaptor;
    }
}
    

### dumb reconstruction of command line
my $commandLine = "$Bin/rnaseq_pipeline.pl";
if($pre){
    $commandLine .= " -pre $pre";
}
if($map){
    $commandLine .= " -map $map";
}
if($config){
    $commandLine .= " -config $config";
}
if($request){
    $commandLine .= " -request $request";
}
if($strand){
    $commandLine .= " -strand $strand";
}
if($samplekey){
    $commandLine .= " -samplekey $samplekey";
}
#if($comparisons){
#    $commandLine .= " -comparisons $comparisons";
#}
if($scheduler){
    $commandLine .= " -scheduler $scheduler";
}
if($r1adaptor){
    $commandLine .= " -r1adaptor $r1adaptor";
}
if($r2adaptor){
    $commandLine .= " -r2adaptor $r2adaptor";
}
if($standard_gene){
    $commandLine .= " -standard_gene";
}
if($standard_transcript){
    $commandLine .= " -standard_transcript";
}
if($star){
    $commandLine .= " -star";
}
if($pass1){
    $commandLine .= " -pass1";
}
if($tophat){
    $commandLine .= " -tophat";
}
if($cufflinks){
    $commandLine .= " -cufflinks";
}
if($dexseq){
    $commandLine .= " -dexseq";
}
if($htseq){
    $commandLine .= " -htseq";
}
if($deseq){
    $commandLine .= " -deseq";
}
if($chimerascan){
    $commandLine .= " -chimerascan";
}
if($star_fusion){
    $commandLine .= " -star_fusion";
}
if($mapsplice){
    $commandLine .= " -mapsplice";
}
if($defuse){
    $commandLine .= " -defuse";
}
if($fusioncatcher){
    $commandLine .= " -fusioncatcher";
}
if($allfusions){
    $commandLine .= " -allfusions";
}
if($kallisto){
    $commandLine .= " -kallisto";
}
if($rsem){
    $commandLine .= " -rsem";
}
if($express){
    $commandLine .= " -express";
}
if($species){
    $commandLine .= " -species $species";
}
if($lncrna){
    $commandLine .= " -lncrna";
}
if($rsync){
    $commandLine .= " -rsync $rsync";
}
if($email){
    $commandLine .= " -email $email";
}


$commandLine.= " -priority_project $priority_project";
$commandLine .= " -priority_group $priority_group";


my $numArgs = $#ARGV + 1;
foreach my $argnum (0 .. $#ARGV) {
    $commandLine .= " $ARGV[$argnum]";
}


if($standard_gene){
    $star = 1;
    $htseq = 1;
}

if($differential_gene){
    $star = 1;
    $htseq = 1;
    $deseq = 1;
}

if($clustering_only){
    $star = 1;
    $htseq = 1;
}

if($standard_transcript){
    $rsem = 1;
}

if($differential_transcript){
    $rsem = 1;
    $deseq = 1;
}


if($star && !$htseq && !$dexseq){
    $htseq = 1;
}

if($allfusions){
    $chimerascan = 1;
    $star_fusion = 1;
    $mapsplice = 1;
    $defuse = 1;
    $fusioncatcher = 1;
}

if($chimerascan || $star_fusion || $mapsplice || $defuse || $fusioncatcher){
    $detectFusions = 1;
}

if($deseq || $dexseq || $htseq || $cufflinks){
    if(!$star && !$tophat && !$rsem){
	$star = 1;

	my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tNO ALIGNER SPECIFIED SO DEFAULTING TO USING STAR\n";
    }
}

`/bin/mkdir -m 775 -p $output`;

## if no sample key specified, check for any that exist and if they do, 
## read them to check for samples that should be excluded
my @sample_keys = $samplekey;
my @excluded_samples = ();
if($deseq || $clustering_only){
    if(!$samplekey){
        @sample_keys = glob("*sample_key*.txt");
        foreach my $fl (@sample_keys){
            print "$fl\n";
        }
    }
    my $keyCount = scalar @sample_keys;
    if($keyCount == 0 && $deseq){
        die "No key file found. Can not run DESeq.\n";
    }
    foreach my $sample_key (@sample_keys){
        print "Checking for exclusions in $sample_key\n";
        open(SK, "$sample_key") or die "Can't open sample key $sample_key $!";
        while(<SK>){
            chomp;
            my @samp_group = split(/\s+/, $_);
            if($samp_group[1] =~ /_exclude_/i){
print "Excluding $samp_group[0]\t$samp_group[1]\n";
                push(@excluded_samples, $samp_group[0]);
            }
        }
        close SK;
    }
}

#die;

my %mapping_samples = ();
my %ism_samples = ();
open(MA, "$map") or die "Can't open mapping file $map $!";
while(<MA>){
    chomp;
    
    my @data = split(/\s+/, $_);
    
    if(grep(/^$data[1]$/, @excluded_samples)){
        print "EXCLUDING sample $data[1]\n";
        next;
    }

    $mapping_samples{$data[1]} = 1;
    if(!-d $data[3]){
	die "$data[3] does not exist\n";
    }
    ## keep track of PE/SE for InsertSizeMetrics
    if($data[4] eq "PE"){
        $ism_samples{$data[1]} = 1;
    }

    `/bin/cp $map $output/`;

}
close MA;

my $BOWTIE2 = '';
my $CUFFLINKS = '';
my $HTSEQ = '';
my $DEXSEQ = '';
my $EXPRESS = '';
my $KALLISTO = '';
my $PICARD = '';
my $CHIMERASCAN = '';
my $MAPSPLICE = '';
my $STAR = '';
my $DEFUSE = '';
my $FUSIONCATCHER = '';
my $TOPHAT = '';
my $PYTHON = '';
my $JAVA = '';
my $PERL = '';
my $R = '';
my $RSEM = '';
my $CUTADAPT = '';
my $STAR_FUSION = '';
my $SINGULARITY = '';
my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";


my %samp_libs_run = ();
my $slr_count = 0;
my %samp_pair = ();

open(LOG, ">$cd\_rnaseq_pipeline.log") or die "can't write to output log";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSTARTING RNASEQ PIPELINE FOR $pre\n";
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCOMMAND LINE: $commandLine\n";

### Check that all programs are available
&verifyConfig($config);

my $GTF = '';
my $DEXSEQ_GTF = '';
my $CHIMERASCAN_INDEX = '';
my $geneNameConversion = '';
my $starDB = '';
my $chrSplits = '';
my $BOWTIE_INDEX = '';
my $BOWTIE2_INDEX = '';
my $TRANS_FASTA_DEDUP = '';
my $TRANS_INDEX = '';
my $TRANS_INDEX_DEDUP = '';
my $REF_SEQ = '';
my $RIBOSOMAL_INTERVALS='';
my $REF_FLAT = '';
my $KALLISTO_INDEX = '';
my $STAR_FUSION_GENOME_LIB = '';
my $CHIMERASCAN_FP_FILTER = '';
my $RSEM_DB = '';
my $QC_BED = '';
my $GMT_DIR = '';

my $STAR_MAX_MEM = 100;

if($species =~ /human|hg19/i){
    $species = 'hg19';
    $REF_SEQ = '/juno/depot/assemblies/H.sapiens/hg19/hg19.fasta';
    $DEXSEQ_GTF = "$Bin/data/human/gencode.v18.annotation_dexseq.gtf";
    $QC_BED = "$Bin/data/hg19/hg19_GENCODE_GENE_V19_comprehensive.bed";
    $CHIMERASCAN_INDEX = '/juno/depot/assemblies/H.sapiens/hg19/index/chimerascan/0.4.5a';
    $BOWTIE_INDEX = '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/1.0.0/hg19_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/hg19_bowtie2';
    $chrSplits = '/juno/depot/assemblies/H.sapiens/hg19/chromosomes';
    $RIBOSOMAL_INTERVALS = "$Bin/data/hg19/ribosomal_hg19.interval_file";
    $REF_FLAT = "$Bin/data/hg19/refFlat__hg19.txt.gz";
    $TRANS_INDEX = '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/gencode.v18.annotation';
    $TRANS_INDEX_DEDUP = '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup';
    $TRANS_FASTA_DEDUP = '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup.fa';
    $geneNameConversion = "$Bin/data/human/gencode18IDToGeneName.txt";
    $KALLISTO_INDEX = '/juno/depot/assemblies/H.sapiens/hg19/index/kallisto/v0.42.1/gencode/v18/gencode.v18.annotation.gtf.fasta.idx';
    $RSEM_DB = '/opt/common/CentOS_6/rsem/RSEM-1.2.25/data/hg19/star/hg19_v19_gencode';
    $STAR_FUSION_GENOME_LIB = "$STAR_FUSION/Hg19_CTAT_resource_lib";
    $CHIMERASCAN_FP_FILTER = "$Bin/data/hg19/hg19_bodymap_false_positive_chimeras.txt";
    $GMT_DIR = "$Bin/data/human/MSigDB/v7.1";
    if($lncrna){
        $GTF = "$Bin/data/human/lncipedia.gtf";
        $starDB = '/juno/depot/assemblies/H.sapiens/hg19/index/star/2.3.0e_r291/LNCipedia';
    } else {
        $GTF = "$Bin/data/human/gencode.v18.annotation.gtf";
	if($r1adaptor){
	    $starDB = '/juno/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang49';
	}
	else{
	    $starDB = '/juno/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang74';
	}
    }

    $STAR_MAX_MEM = 35;
}
elsif($species =~ /mouse|mm10/i){
    $species = 'mm10';
    $REF_SEQ = '/juno/depot/assemblies/M.musculus/mm10/mm10.fasta';
    $QC_BED = "$Bin/data/mm10/mm10_GENCODE_VM9_basic.bed";
    $GTF = "$Bin/data/mm10/gencode.vM8.annotation.gtf";
    #$DEXSEQ_GTF = "$Bin/data/Mus_musculus.GRCm38.80_canonical_chromosomes.dexseq.gtf";
    $geneNameConversion = "$Bin/data/mm10/gencodeM8IDToGeneName.txt";
    $BOWTIE_INDEX = '/juno/depot/assemblies/M.musculus/mm10/index/bowtie/1.1.1/mm10_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/mm10_bowtie2';
    $chrSplits = '/juno/depot/assemblies/M.musculus/mm10/chromosomes';
    $TRANS_INDEX = '/juno/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/transcriptome/gencode/vM8/gencode.vM8.annotation';
    $TRANS_INDEX_DEDUP = '/juno/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/transcriptome/gencode/vM8/deduplicated/gencode.vM8.annotation.dedup';
    $TRANS_FASTA_DEDUP = '/juno/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/transcriptome/gencode/vM8/deduplicated/gencode.vM8.annotation.dedup.fa';
    $RIBOSOMAL_INTERVALS = "$Bin/data/mm10/ribosomal_mm10.interval_file";
    $REF_FLAT = "$Bin/data/mm10/refFlat__mm10.txt.gz";
    $KALLISTO_INDEX = '/juno/depot/assemblies/M.musculus/mm10/index/kallisto/v0.42.1/gencode/vM8/gencode.vM8.annotation.gtf.fasta.idx';
    $RSEM_DB = '/opt/common/CentOS_6/rsem/RSEM-1.2.25/data/GRCm38/star/GRCm38_vM8_gencode';
    $GMT_DIR = "$Bin/data/mouse/MSigDB/v7.1";
    if($r1adaptor){
	$starDB = '/juno/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/gencode/vM8/overhang49';
    }
    else{
	$starDB = '/juno/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/gencode/vM8/overhang74';
    }

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /mm9/i){
    $species = 'mm9';
    $REF_SEQ = '/juno/depot/assemblies/M.musculus/mm9/mm9.fasta';
    $GTF = "$Bin/data/mm9/Mus_musculus.NCBIM37.67_ENSEMBL.gtf";
    $DEXSEQ_GTF = "$Bin/data/mm9/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf";
    $geneNameConversion = "$Bin/data/mm9/mm9Ensembl67IDToGeneName.txt";
    $BOWTIE_INDEX = '/juno/depot/assemblies/M.musculus/mm9/index/bowtie/1.0.0/mm9_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/mm9_bowtie2';
    $chrSplits = '/juno/depot/assemblies/M.musculus/mm9/chromosomes';
    $TRANS_INDEX = '/juno/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/transcriptome/ensembl/vTBD/ensembl';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = "$Bin/data/mm9/ribosomal_MM9_assemblies.interval_file";
    $REF_FLAT = "$Bin/data/mm9/refFlat__mm9.txt.gz";
    $KALLISTO_INDEX = '';
    $RSEM_DB = '/opt/common/CentOS_6/rsem/RSEM-1.2.25/data/mm9/star/mm9_vm1_gencode';
    $GMT_DIR = "$Bin/data/mouse/MSigDB/v7.1";
    if($r1adaptor){
	$starDB = '/juno/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang49';
    }
    else{
	$starDB = '/juno/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang74';
    }

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /human-mouse|mouse-human|hybrid/i){
    $species = 'hybrid';
    $REF_SEQ = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta';
    $GTF = "$Bin/data/human/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/human/gencode.v18.annotation_dexseq.gtf";
    $geneNameConversion = "$Bin/data/human/gencode18IDToGeneName.txt";
    $RIBOSOMAL_INTERVALS = "$Bin/data/hg19_mm10/ribosomal_hg19_mm10.interval_file";
    $REF_FLAT = "$Bin/data/hg19/refFlat__hg19.txt.gz";
    $RSEM_DB = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/rsem/1.2.25/gencode/v18/hg19_mm10';
    $BOWTIE_INDEX = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/bowtie/1.1.1/hg19_mm10_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/bowtie/2.2.4/hg19_mm10_bowtie2';
    $KALLISTO_INDEX = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/kallisto/v0.42.1/gencode/v18/gencode.v18.annotation.gtf.fasta.idx';
    $STAR_FUSION_GENOME_LIB = "$STAR_FUSION/Hg19_CTAT_resource_lib";
    $CHIMERASCAN_INDEX = "/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/chimerascan/0.4.5a";
    $CHIMERASCAN_FP_FILTER = "$Bin/data/hg19/hg19_bodymap_false_positive_chimeras.txt";
    $chrSplits = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/chromosomes';
    $QC_BED = "$Bin/data/hg19/hg19_GENCODE_GENE_V19_comprehensive.bed";
    $GMT_DIR = "$Bin/data/human/MSigDB/v7.1";
    if($r1adaptor){
	$starDB = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang49';
    }
    else{
	$starDB = '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang74';
    }

    $STAR_MAX_MEM = 60;
}
elsif($species =~ /rat|Rn6/i){
    $species = 'Rn6';
    $REF_SEQ = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/Rnor_6.0.fasta';
    $GTF = "$Bin/data/rn6/Rattus_norvegicus.Rnor_6.0.92.gtf";

    $geneNameConversion = "$Bin/data/rn6/rn6Ensembl92IDToGeneName.txt";


    $BOWTIE_INDEX = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/bowtie/1.1.1/Rnor_6.0_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/bowtie/2.2.4/Rnor_6.0_bowtie2';
    $chrSplits = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/chromosomes';
    $TRANS_INDEX = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/bowtie/2.2.4/transcriptome/ensembl/v92/Rattus_norvegicus.Rnor_6.0.92';
    $TRANS_INDEX_DEDUP = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/bowtie/2.2.4/transcriptome/ensembl/v92/deduplicated/Rattus_norvegicus.Rnor_6.0.92.gtf.dedup';
    $TRANS_FASTA_DEDUP = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/bowtie/2.2.4/transcriptome/ensembl/v92/deduplicated/Rattus_norvegicus.Rnor_6.0.92.gtf.dedup.fa';
    

    $RIBOSOMAL_INTERVALS = "$Bin/data/rn6/ribosomal_Rnor6.interval_file";


    $REF_FLAT = "$Bin/data/rn6/refFlat__rn6.txt";
    $KALLISTO_INDEX = '';
    $RSEM_DB = '';

    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/star/2.4.1d/ensembl/v92/overhang49';
    }
    else{
        $starDB = '/juno/depot/assemblies/R.norvegicus/Rnor_6.0/index/star/2.4.1d/ensembl/v92/overhang74';
    }

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /zebrafish|z11|GRCz11/i){
    $species = 'z11';
    $REF_SEQ = '/juno/depot/assemblies/D.rerio/GRCz11/GRCz11.fasta';
    $GTF = '/juno/depot/assemblies/D.rerio/GRCz11/Danio_rerio.GRCz11.97.gtf';
    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/D.rerio/GRCz11/index/star/2.4.1d/ensembl/v97/overhang49';
    } else {
        $starDB = '/juno/depot/assemblies/D.rerio/GRCz11/index/star/2.4.1d/ensembl/v97/overhang74';
    }
    $chrSplits = '/juno/depot/assemblies/D.rerio/GRCz11/chromosomes';
    $geneNameConversion = "$Bin/data/GRCz11/IDToGeneName.txt";
    $TRANS_INDEX = '';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = "$Bin/data/GRCz11/ribosomal_GRCz11.97.interval_file";
    $REF_FLAT = "$Bin/data/GRCz11/refFlat__GRCz11.txt.gz";
    $RSEM_DB = "";
    $KALLISTO_INDEX = "";
    $QC_BED = "$Bin/data/GRCz11/Danio_rerio.GRCz11.97.bed";

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /zv10/i){
    $species = 'zv10';
    $REF_SEQ = '/juno/depot/assemblies/D.rerio/GRCz10/GRCz10.fasta';
    $GTF = "$Bin/data/zv10/GRCz10.gtf";
    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/D.rerio/GRCz10/index/star/2.5.0a/ensembl/v83/overhang49';
    } else {
        $starDB = '/juno/depot/assemblies/D.rerio/GRCz10/index/star/2.5.0a/ensembl/v83/overhang74';
    }
    $chrSplits = '/juno/depot/assemblies/D.rerio/GRCz10/chromosomes';
    $geneNameConversion = "$Bin/data/zv10/zv10EnsemblIDtoGeneName.txt";
    $TRANS_INDEX = '';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = "$Bin/data/zv10/ribosomal_zv10.interval_file";
    $REF_FLAT = "$Bin/data/zv10/refFlat__zv10.txt.gz";
    $RSEM_DB = '/opt/common/CentOS_6/rsem/RSEM-1.2.25/data/zv10/star/zv10_v83_ensembl';
    $KALLISTO_INDEX = '';
    $QC_BED = "$Bin/data/zv10/zv10_ENSEMBL_20180417.bed";

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /zv9/i){
    $species = 'zv9';
    $REF_SEQ = '/juno/depot/assemblies/D.rerio/zv9/zv9.fasta';
    $GTF = "$Bin/data/zv9/zv9.gtf";
    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/D.rerio/zv9/index/star/2.4.1d/ensembl/v79/overhang49';
    } else {
        $starDB = '/juno/depot/assemblies/D.rerio/zv9/index/star/2.4.1d/ensembl/v79/overhang74';
    }
    $chrSplits = '/juno/depot/assemblies/D.rerio/zv9/chromosomes';
    $geneNameConversion = "Bin/data/zv9/zv9EnsemblIDtoGeneName.txt";
    $TRANS_INDEX = '';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = '';
    $REF_FLAT = '';
    $KALLISTO_INDEX = '';

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /fly|dmel6|dm6/i){
    $species = 'dm6';
    $REF_SEQ = '/juno/depot/assemblies/D.melanogaster/dm6/dm6.fasta';
    $GTF = "/juno/depot/assemblies/D.melanogaster/dm6/dmel-all-r6.33.gtf";
    $QC_BED = "/juno/depot/assemblies/D.melanogaster/dm6/dm6_refseq_all.bed";
    $DEXSEQ_GTF = "";
    $geneNameConversion = "/juno/depot/assemblies/D.melanogaster/dm6/IDToGeneName.txt";
    $BOWTIE_INDEX = '/juno/depot/assemblies/D.melanogaster/dm6/index/bowtie/1.1.1/dm6_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/D.melanogaster/dm6/index/bowtie/2.2.4/dm6_bowtie2';
    $RSEM_DB = '/juno/depot/assemblies/D.melanogaster/dm6/index/rsem/1.3.3/flybase/6.33/star/dm6';
    $chrSplits = '/juno/depot/assemblies/D.melanogaster/dm6/chromosomes';
    $TRANS_INDEX = '';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = '/juno/depot/assemblies/D.melanogaster/dm6/ribosomal_dm6.interval_file';
    $REF_FLAT = '/juno/depot/assemblies/D.melanogaster/dm6/dm6_refFlat.txt';
    $KALLISTO_INDEX = '';

    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/D.melanogaster/dm6/index/star/2.4.1d/flybase/6.33/overhang49';
    }
    else{
        $starDB = '/juno/depot/assemblies/D.melanogaster/dm6/index/star/2.4.1d/flybase/6.33/overhang74';
    }

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /WBcel235/i){
    $species = 'WBcel235';
    $REF_SEQ = '/juno/depot/assemblies/custom/C.elegans/WBcel235/WBcel235.fasta';
    $GTF = "$Bin/data/WBcel235/ensembl_WBcel235.gtf";
    $DEXSEQ_GTF = '';
    $CHIMERASCAN_INDEX = '';
    $BOWTIE_INDEX = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/bowtie/1.1.1/WBcel235_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/bowtie/2.2.4/WBcel235_bowtie2';
    $chrSplits = '/juno/depot/assemblies/custom/C.elegans/WBcel235/chromosomes';
    $RIBOSOMAL_INTERVALS = "$Bin/data/WBcel235/ribosomal_WBcel235.interval_file";
    $REF_FLAT = "$Bin/data/WBcel235/refFlat__ce11.txt.gz";
    $TRANS_INDEX = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/bowtie/2.2.4/transcriptome/ensembl/v20151123/genes';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $geneNameConversion = "$Bin/data/WBcel235/WBcel235EnsemblIDToGeneName.txt";
    $KALLISTO_INDEX = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/kallisto/v0.42.1/ensembl/v20151123/genes.gtf.fasta.idx';
    $RSEM_DB = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/rsem/1.2.25/star/2.4.1d/WBcel235_ensembl';
    $STAR_FUSION_GENOME_LIB = '';
    $CHIMERASCAN_FP_FILTER = '';

    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/star/2.4.1d/ensembl/v20151123/overhang49';
    }
    else{
        $starDB = '/juno/depot/assemblies/custom/C.elegans/WBcel235/index/star/2.4.1d/ensembl/v20151123/overhang74';
    }

    $STAR_MAX_MEM = 30;
}
elsif($species =~ /yeast|s288c|R64|sacCer3/i){
    $species = 'S288C_R64';
    $REF_SEQ = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/S288C_R64.fasta';
    $GTF = '/juno/depot/annotation/S.cerevisiae/ensembl/v94/Saccharomyces_cerevisiae.R64-1-1.94.gtf';
    $QC_BED = "$Bin/data/S288C_R64/Saccharomyces_cerevisiae.R64-1-1.94.bed";
    $DEXSEQ_GTF = '';
    $CHIMERASCAN_INDEX = '';
    $BOWTIE_INDEX = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/bowtie/1.1.1/S288C_R64_bowtie';
    $BOWTIE2_INDEX = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/bowtie/2.2.4/S288C_R64_bowtie';
    $chrSplits = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/chromosomes';
    
    $RIBOSOMAL_INTERVALS = "$Bin/data/S288C_R64/ribosomal_S288C_R64_assemblies.interval_file";
    $REF_FLAT = "$Bin/data/S288C_R64/Saccharomyces_cerevisiae.R64-1-1.94.dedup.refflat.gz";


    $TRANS_INDEX = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/bowtie/2.2.4/transcriptome/ensembl/v94/Saccharomyces_cerevisiae.R64-1-1.94';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';

    $geneNameConversion = "$Bin/data/S288C_R64/S288C_R64EnsemblIDToGeneName.txt";

    $KALLISTO_INDEX = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/kallisto/v0.42.1/ensembl/v94/Saccharomyces_cerevisiae.R64-1-1.94';
    $RSEM_DB = '';

    $STAR_FUSION_GENOME_LIB = '';
    $CHIMERASCAN_FP_FILTER = '';

    if($r1adaptor){
        $starDB = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/star/2.4.1d/ensembl/v94/overhang49';
    }
    else{
        $starDB = '/juno/depot/assemblies/S.cerevisiae/S288C_R64/index/star/2.4.1d/ensembl/v94/overhang74';
    }

    $STAR_MAX_MEM = 30;
}

open(MA, "$map") or die "Can't open mapping file $map $!";
while(<MA>){
    chomp;

    my @data = split(/\s+/, $_);

    if(!-d $data[3]){
	die "$data[3] does not exist\n";
    }
}
close MA;

my $htseq_stranded = '';
my $picard_strand_specificity = '';
my $star_outSAMstrandField = ''; ### compatibility with cufflinks;
my $cufflinks_lib_type = '';
if($strand =~ /none/i){
    $htseq_stranded = 'no';
    $picard_strand_specificity = 'NONE';
    $star_outSAMstrandField = '--outSAMstrandField intronMotif';
    $cufflinks_lib_type = '--library-type fr-unstranded';
}
elsif($strand =~ /forward/i){
    $htseq_stranded = 'yes';
    $picard_strand_specificity = 'FIRST_READ_TRANSCRIPTION_STRAND';
    $cufflinks_lib_type = '--library-type fr-firststrand';
}
elsif($strand =~ /reverse/i){
    $htseq_stranded = 'reverse';
    $picard_strand_specificity = 'SECOND_READ_TRANSCRIPTION_STRAND';
    $cufflinks_lib_type = '--library-type fr-secondstrand';
}


`/bin/mkdir -m 775 -p $output/intFiles`; 
`/bin/mkdir -m 775 -p $output/progress`;
`/bin/mkdir -m 775 -p $output/gene`;
`/bin/mkdir -m 775 -p $output/gene/alignments`;
`/bin/mkdir -m 775 -p $output/metrics`;
`/bin/mkdir -m 775 -p $output/metrics/crm`;
`/bin/mkdir -m 775 -p $output/metrics/images`;

my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

my @syncJobs = ();
open(IN, "$map") or die "Can't open $map $!";
while(<IN>){
    chomp;

    my @data = split(/\s+/, $_);
    if(grep(/^$data[1]$/, @excluded_samples)){
        print "EXCLUDING sample $data[1]\n";
        next;
    }

    if($data[1] =~ /^\d+/){
	$data[1] = "s_" . "$data[1]";
    }

    ### sometimes in the mapping file, 
    ### the lib, sample, and run id isn't a unique identifier
    if($samp_libs_run{$data[1]}{$data[0]}{$data[2]}){
	$slr_count++;
my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tWARNING: $data[1]\t$data[0]\t$data[2] ISN'T UNIQUE; ";
	$data[2] = "$data[2]\_$slr_count";
	print LOG "WRITING INSTEAD TO $data[1]\/$data[0]\/$data[2]\n";
    }

    `/bin/mkdir -m 775 -p $output/intFiles/$data[1]`;
    `/bin/mkdir -m 775 -p $output/intFiles/$data[1]/$data[0]/$data[2]`;
    $samp_libs_run{$data[1]}{$data[0]}{$data[2]} = 1;
    `ln -s $data[3]/* $output/intFiles/$data[1]/$data[0]/$data[2]/`;
    chdir "$output/intFiles/$data[1]/$data[0]/$data[2]";

    opendir(workDir, "./");
    my @unsorted = readdir workDir;
    closedir workDir;
    my @files = sort @unsorted;
    
    open(OUT, ">files_$data[1]\_$data[0]\_$data[2]");
    foreach my $file (@files){
	if($file =~ /fastq\.gz$/ && $file =~ m/^(.*)(R\d+)(.*)$/){
	    if($2 eq 'R1'){
		my $file_R2 = $file;
		$file_R2 =~ s/^(.*)R1(.*)$/$1R2$2/;
		
		if($data[4] =~ /pe/i){
		    $samp_pair{$data[1]} = "PE";
		    print OUT "$file\t$file_R2\n";
		}
		elsif($data[4] =~ /se/i){
		    $samp_pair{$data[1]} = "SE";
		    print OUT "$file\n";
		}
		else{ 
		    die "CAN'T DETERMINE WHETHER RUN IS PAIRED OR SINGLE ENDED for sample $data[1]\n"; 
		}
	    }
	}
    }
    close OUT;
    chdir $curDir;
}
close IN;

## QC sample key(s) and comparison file(s)
if($deseq || $clustering_only){
    my $cluster_only = '';
    if($clustering_only){
        $cluster_only = '--clustering_only TRUE';
    }
    my $qc = `$singularityParams $R/Rscript $Bin/RunDE_inputQC.R --bin $Bin --Rlibs $Bin/lib/R-4.0.0 --pre $pre $cluster_only 2>&1`;
    my $ec = $? >> 8;
    print("$qc");
    if($ec == 1){
        die("Error(s) found in DESeq input files.");
    }
}

my @crm = ();
my @crm_tophat = ();
my @asm = ();
my @asm_tophat = ();
my @ism = ();
my @ism_tophat = ();
my $ran_tophat = 0;
my @cag_jids = ();
my $ran_cag = 0;
my @thro_jids = ();
my $ran_tophatcrm = 0;
my @tophatcrm_jids = ();
my $ran_tophatasm = 0;
my $ran_tophatism = 0;
my $ran_tophatism_merge = 0;
my @tophatasm_jids = ();
my @tophatism_jids = ();
my $ran_rseqc = 0;
my $ran_rseqc_merge = 0;
my $ran_plots = 0;
my @rseqc_jids = ();
my $ran_picard_metrics = 0;
my $ran_starcrm = 0;
my @starcrm_jids = ();
my $ran_starasm = 0;
my @starasm_jids = ();
my $ran_starism = 0;
my $ran_starism_merge = 0;
my @starism_jids = ();
my $ran_star = 0;
my @staraddrg_jids = ();
my $ran_tophathtseq = 0;
my @tophathtseq_jids = ();
my $ran_tophatdexseq = 0;
my @tophatdexseq_jids = ();
my $ran_starhtseq = 0;
my @starhtseq_jids = ();
my $ran_stardexseq = 0;
my @stardexseq_jids = ();
my $ran_rceg = 0;
my @rce_jids = ();
my $ran_kallisto = 0;
my @kallisto_jids = ();
my @qcplot_jids = ();
my @qcpdf_jids = ();

foreach my $sample (keys %samp_libs_run){
    my @R1 = ();
    my @R2 = ();
    my $readsFlag = 0;
    my $grl = 0;
    my $minReadLength = -1;
    my $sampReadLength = 0;
    foreach my $lib (keys %{$samp_libs_run{$sample}}){
	foreach my $run (keys %{$samp_libs_run{$sample}{$lib}}){	    
	    open(READS, "$output/intFiles/$sample/$lib/$run/files_$sample\_$lib\_$run") || die "Can't open $output/intFiles/$sample/$lib/$run/files_$sample\_$lib\_$run $!";
	    while(<READS>){
		chomp;
		    
		my @readPair = split(/\s+/, $_);
		
		### NOTE: JUST CHECKING TO SEE IF IT EXISTS
		###       HOWEVER DOES NOT GURANTEE THAT IT'S NON-EMPTY
		if(-e "$output/intFiles/$sample/$lib/$run/$readPair[0]"){
		    push @R1, "$output/intFiles/$sample/$lib/$run/$readPair[0]";
		    if($samp_pair{$sample} eq "PE"){
			if(!-e "$output/intFiles/$sample/$lib/$run/$readPair[1]"){
			    $readsFlag = 1;
			}
			push @R2, "$output/intFiles/$sample/$lib/$run/$readPair[1]";		    
		    }

		    if($grl == 0){
			my $readO = `/bin/gzip -cd $output/intFiles/$sample/$lib/$run/$readPair[0] | head -2`;
			### discard reads less than half of original read length by cutadapt
			#my $readO = `/usr/bin/head -2 $output/intFiles/$readPair[0]`;
			chomp $readO;
			my @dataO = split(/\n/, $readO);
			my $readLength = length($dataO[1]);
                        if($sampReadLength > 0 & $readLength != $sampReadLength){
                            print "WARNING: fastqs have variable read lengths\n";
                        }
                        $sampReadLength = $readLength;
			$minReadLength = int(0.5*$readLength);
			$grl = 1;
		    }
		}
		else{
		    $readsFlag = 1;
		}
	    }
	    close READS;
	}
    }

    if($readsFlag == 1){
	my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSKIPPING SAMPLE $sample ANALYSIS BECAUSE CAN'T LOCATE ALL LISTED READS";
	next;
    }

    my $r1_gz_files = join(",", @R1);
    my $r2_gz_files = join(",", @R2);

    my @gz_jids = ();
    my $ran_gz = 0;
    if($r1adaptor){
	my $r1_gz_files_TRIM = join(" ", @R1);
	my $r2_gz_files_TRIM = join(" ", @R2);
	my $ran_zcat = 0;
	my @zcat_jids = ();
	if(!-e "$output/progress/$pre\_$uID\_ZCAT_$sample\_R1.done"){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$sample\_R1", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT_$sample\_R1.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r1_gz_files_TRIM ">$output/intFiles/$sample/$sample\_R1.fastq"`;
	    `/bin/touch $output/progress/$pre\_$uID\_ZCAT_$sample\_R1.done`;
	    push @zcat_jids, "$pre\_$uID\_ZCAT_$sample\_R1";
	    $ran_zcat = 1;
	}

	if($samp_pair{$sample} eq "PE"){
	    if(!-e "$output/progress/$pre\_$uID\_ZCAT_$sample\_R2.done"){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$sample\_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT_$sample\_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r2_gz_files_TRIM ">$output/intFiles/$sample/$sample\_R2.fastq"`;
		###`$Bin/qSYNC $pre\_$uID\_ZCAT_$sample\_R2`;
		`/bin/touch $output/progress/$pre\_$uID\_ZCAT_$sample\_R2.done`;
		push @zcat_jids, "$pre\_$uID\_ZCAT_$sample\_R2";
		$ran_zcat = 1;
	    }
	}
	
	my $processR1 = "-a $r1adaptor -O 10";
	my $processR2 = "-a $r2adaptor -O 10";
	my $zcatj = join(",", @zcat_jids);
	my $ran_ca = 0;
	my @ca_jids = ();
	if($samp_pair{$sample} eq "PE"){
	    if(!-e "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.done" || $ran_zcat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$sample\_R1", job_hold => "$zcatj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR1 --paired-output $output/intFiles/$sample/$sample\_R2_TEMP.fastq -o $output/intFiles/$sample/$sample\_R1_TEMP.fastq $output/intFiles/$sample/$sample\_R1.fastq $output/intFiles/$sample/$sample\_R2.fastq >$output/intFiles/$sample/$sample\_R1\_CUTADAPT\_STATS.txt"`;
		`/bin/touch $output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.done`;
		push @ca_jids, "$pre\_$uID\_CUTADAPT_$sample\_R1";
		push @cag_jids, "$pre\_$uID\_CUTADAPT_$sample\_R1";
		$ran_ca = 1;
		$ran_cag = 1;
	    }

	    if(!-e "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R2.done" || $ran_zcat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$sample\_R2", job_hold => "$zcatj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR2 --paired-output $output/intFiles/$sample/$sample\_R1_CT.fastq -o $output/intFiles/$sample/$sample\_R2_CT.fastq $output/intFiles/$sample/$sample\_R2_TEMP.fastq $output/intFiles/$sample/$sample\_R1_TEMP.fastq >$output/intFiles/$sample/$sample\_R2\_CUTADAPT\_STATS.txt"`;
		`/bin/touch $output/progress/$pre\_$uID\_CUTADAPT_$sample\_R2.done`;
		push @ca_jids, "$pre\_$uID\_CUTADAPT_$sample\_R2";
		push @cag_jids, "$pre\_$uID\_CUTADAPT_$sample\_R2";
		$ran_ca = 1;
		$ran_cag = 1;
	    }
	   
	    my $caj = join(",", @ca_jids);
	    if(!-e "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.done" || $ran_ca){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GZIP_$sample\_R1", job_hold => "$caj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/gzip $output/intFiles/$sample/$sample\_R1_CT.fastq`;
		`/bin/touch $output/progress/$pre\_$uID\_GZIP_$sample\_R1.done`;
		push @gz_jids, "$pre\_$uID\_GZIP_$sample\_R1";
		$ran_gz = 1;
	    }

	    if(!-e "$output/progress/$pre\_$uID\_GZIP_$sample\_R2.done" || $ran_ca){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GZIP_$sample\_R2", job_hold => "$caj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GZIP_$sample\_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/gzip $output/intFiles/$sample/$sample\_R2_CT.fastq`;
		`/bin/touch $output/progress/$pre\_$uID\_GZIP_$sample\_R2.done`;
		push @gz_jids, "$pre\_$uID\_GZIP_$sample\_R2";
		$ran_gz = 1;
	    }
	}
	else{
	    if(!-e "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.done" || $ran_zcat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$sample\_R1", job_hold => "$zcatj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR1 -o $output/intFiles/$sample/$sample\_R1_CT.fastq $output/intFiles/$sample/$sample\_R1.fastq >$output/intFiles/$sample/$sample\_R1\_CUTADAPT\_STATS.txt"`;
		`/bin/touch $output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.done`;
		push @ca_jids, "$pre\_$uID\_CUTADAPT_$sample\_R1";
		push @cag_jids, "$pre\_$uID\_CUTADAPT_$sample\_R1";
		$ran_ca = 1;
		$ran_cag = 1;
	    }

    	    my $caj = join(",", @ca_jids);
	    if(!-e "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.done" || $ran_ca){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GZIP_$sample\_R1", job_hold => "$caj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/gzip $output/intFiles/$sample/$sample\_R1_CT.fastq`;
		`/bin/touch $output/progress/$pre\_$uID\_GZIP_$sample\_R1.done`;
		push @gz_jids, "$pre\_$uID\_GZIP_$sample\_R1";
		$ran_gz = 1;
	    }
	}

	$r1_gz_files = "$output/intFiles/$sample/$sample\_R1_CT.fastq.gz";
	$r2_gz_files = "$output/intFiles/$sample/$sample\_R2_CT.fastq.gz";
    }

    my $gzj = join(",", @gz_jids);
    if($tophat){
	`/bin/mkdir -m 775 -p $output/intFiles/tophat2`;
	`/bin/mkdir -m 775 -p $output/intFiles/tophat2/$sample`;
	`/bin/mkdir -m 775 -p $output/gene/alignments/tophat2`;
	`/bin/mkdir -m 775 -p $output/metrics/tophat2`;

	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}
	
	my $tophatj = '';
	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_$sample.done" || $ran_gz){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_$sample", job_hold => "$gzj", cpu => "6", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $TOPHAT/tophat2 -p 6 -r 70 --mate-std-dev 90 --GTF $GTF --transcriptome-index=$TRANS_INDEX -o $output/intFiles/tophat2/$sample --rg-id $sample\_1 --rg-sample $sample --rg-library $sample\_1 --rg-platform Illumina --rg-platform-unit $sample\_1 $BOWTIE2_INDEX $inReads`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_$sample.done`;
	    $tophatj = "$pre\_$uID\_TOPHAT_$sample";
	    $ran_tophat = 1;
	}

	my $ran_reorder = 0;
	my $reorderj = '';
        #if(!-e "$output/progress/$pre\_$uID\_TOPHAT_$sample.done" || $ran_tophat){   ## THIS LOOKS LIKE A BUG - NOT SURE IF I INTRODUCED IT, BUT I'M CHANGING IT - CJ
        if(!-e "$output/progress/$pre\_$uID\_REORDER_$sample.done" || $ran_tophat){
            sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_REORDER_$sample", job_hold => "$tophatj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_REORDER_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar ReorderSam I=$output/intFiles/tophat2/$sample/accepted_hits.bam O=$output/gene/alignments/tophat2/$pre\_$sample\.bam REFERENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true`;
	    `/bin/touch $output/progress/$pre\_$uID\_REORDER_$sample.done`; ## bug fix - CJ
            `ln -s $output/gene/alignments/tophat2/$pre\_$sample\.bai $output/gene/alignments/tophat2/$pre\_$sample\.bam.bai`;
	    $reorderj = "$pre\_$uID\_REORDER_$sample";
	    $ran_reorder = 1;
            push @thro_jids, $reorderj;
            push @syncJobs, $reorderj;
	}

        next if($alignment_only);
	
	if($cufflinks){
	    `/bin/mkdir -m 775 -p $output/cufflinks`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2/$sample`;
	    if(!-e "$output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.done" || $ran_tophat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUFFLINKS_TOPHAT_$sample", job_hold => "$tophatj", cpu => "5", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $CUFFLINKS/cufflinks -q -p 12 --no-update-check $cufflinks_lib_type -N -G $GTF -o $output/cufflinks/tophat2/$sample $output/intFiles/tophat2/$sample/accepted_hits.bam`;
		`/bin/touch $output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.done`;
		push @syncJobs, "$pre\_$uID\_CUFFLINKS_TOPHAT_$sample";
	    }
	}
	
	if($htseq || $dexseq){
	    my $ran_tophatqns = 0;
	    my $tophatqnsj = '';
	    if(!-e "$output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.done" || $ran_tophat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QNS_TOPHAT_$sample", job_hold => "$tophatj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles INPUT=$output/intFiles/tophat2/$sample/accepted_hits.bam OUTPUT=$output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam SORT_ORDER=queryname TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT`;
		`/bin/touch $output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.done`;
		$tophatqnsj = "$pre\_$uID\_QNS_TOPHAT_$sample";
		$ran_tophatqns = 1;
	    }
	    
	    if($htseq){
		`/bin/mkdir -m 775 -p $output/gene/counts_gene`;
		`/bin/mkdir -m 775 -p $output/gene/counts_gene/tophat2`;
		
		if(!-e "$output/progress/$pre\_$uID\_HT_TOPHAT_$sample.done" || $ran_tophatqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HT_TOPHAT_$sample", job_hold => "$tophatqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_HT_TOPHAT_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $PYTHON/python $HTSEQ/htseq-count -m intersection-strict -s $htseq_stranded -t exon $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $GTF > $output/gene/counts_gene/tophat2/$sample.htseq_count"`;
		    `/bin/touch $output/progress/$pre\_$uID\_HT_TOPHAT_$sample.done`;
		    push @tophathtseq_jids, "$pre\_$uID\_HT_TOPHAT_$sample";
		    $ran_tophathtseq = 1;
		}
	    }
	    
	    if($dexseq){
		`/bin/mkdir -m 775 -p $output/exon`;
		`/bin/mkdir -m 775 -p $output/exon/counts_exon`;
		`/bin/mkdir -m 775 -p $output/exon/counts_exon/tophat2`;
		
     		if(!-e "$output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.done" || $ran_tophatqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEXSEQ_TOPHAT_$sample", job_hold => "$tophatqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $output/exon/counts_exon/tophat2/$sample.dexseq_count`;
		    `/bin/touch $output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.done`;
		    push @tophatdexseq_jids, "$pre\_$uID\_DEXSEQ_TOPHAT_$sample";
		    $ran_tophatdexseq = 1;
		}
	    }
	}

	`/bin/mkdir -m 775 -p $output/metrics/tophat2/crm`;
	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.done" || $ran_reorder){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_CRM_$sample", job_hold => "$reorderj", cpu => "1", mem => "3", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/gene/alignments/tophat2/$pre\_$sample\.bam O=$output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/tophat2/crm/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=$picard_strand_specificity METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.done`;
	    push @tophatcrm_jids, "$pre\_$uID\_TOPHAT_CRM_$sample";
	    $ran_tophatcrm = 1;
	}	
	push @crm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.done" || $ran_reorder){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_ASM_$sample", job_hold => "$reorderj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/gene/alignments/tophat2/$pre\_$sample\.bam OUTPUT=$output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.done`;
	    push @tophatasm_jids, "$pre\_$uID\_TOPHAT_ASM_$sample";
	    $ran_tophatasm = 1;
	}
	push @asm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt";

        if(exists $ism_samples{$sample}){
            if(!-e "$output/progress/$pre\_$uID\_TOPHAT_ISM_$sample.done" || $ran_reorder){
                sleep(3);
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_ISM_$sample", job_hold => "$reorderj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_ISM_$sample.log");
                my $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE INPUT=$output/gene/alignments/tophat2/$pre\_$sample\.bam OUTPUT=$output/intFiles/tophat2/$pre\_$sample\_InsertSizeMetrics.txt HISTOGRAM_FILE=$output/intFiles/tophat2/$sample/$pre\_$sample\_InsertSizeHistogram.txt`;
                `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_ISM_$sample.done`;
                push @tophatism_jids, "$pre\_$uID\_TOPHAT_ISM_$sample";
                $ran_tophatism = 1;
            }
            push @ism_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_InsertSizeMetrics.txt";
        }

        if(!-e "$output/progress/$pre\_$uID\_RSEQC_TOPHAT_$sample.done" || $ran_reorder){
            `ln -s $output/gene/alignments/tophat2/$pre\_$sample\.bai $output/gene/alignments/tophat2/$pre\_$sample\.bam.bai`;
            sleep(3);
            # only allow a limited number of rseqc jobs to run at the same time, so added additional job holding
            my $job_count = scalar @rseqc_jids;
            my $batch_job_hold = "";
            if($job_count >= $rseqc_batch_size)
            {
                $batch_job_hold = $rseqc_jids[$job_count - $rseqc_batch_size];
            }
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_TOPHAT_$sample", job_hold => "$reorderj,$batch_job_hold", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_TOPHAT_$sample.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl -I $Bin/lib $Bin/qc/rseqc.pl -pre $pre -config $config -bam $output/gene/alignments/tophat2/$pre\_$sample\.bam -bed $QC_BED -sample $sample -intdir $output/intFiles/$sample -outdir $output/metrics/images -progdir $output/progress -scheduler $scheduler -readlen $sampReadLength -layout $samp_pair{$sample} -sync`;
            `/bin/touch $output/progress/$pre\_$uID\_RSEQC_TOPHAT_$sample.done`;
            push @rseqc_jids, "$pre\_$uID\_RSEQC_TOPHAT_$sample";
            push @qcpdf_jids, "$pre\_$uID\_RSEQC_TOPHAT_$sample";
            $ran_rseqc = 1;
        }
    }

    my $starOut = '';    
    if($star){
        $ran_star = 1;
	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}

	my $ran_star1p = 0;
	my $star1pj = '';
	if(!-e "$output/progress/$pre\_$uID\_STAR_1PASS_$sample.done" || $ran_gz){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_1PASS_$sample", job_hold => "$gzj", cpu => "12", mem => "$STAR_MAX_MEM", cluster_out => "$output/progress/$pre\_$uID\_STAR_1PASS_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_1PASS_ $star_outSAMstrandField --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_1PASS_$sample.done`;
	    $star1pj = "$pre\_$uID\_STAR_1PASS_$sample";
	    $ran_star1p = 1;
	}

	my $ran_sp = 0;
	my $spj = '';
	if($pass1){	    
	    if(!-e "$output/progress/$pre\_$uID\_SP_$sample.done" || $ran_star1p){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SP_$sample", job_hold => "$star1pj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_SP_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_1PASS_Aligned.out.sam`;
		`/bin/touch $output/progress/$pre\_$uID\_SP_$sample.done`;
		$spj = "$pre\_$uID\_SP_$sample";
		$ran_sp = 1;
	    }
	    $starOut = "$output/intFiles/$sample/$sample\_STAR_1PASS_Aligned.out.sam_filtered.sam";
	}
	else{
	    `/bin/mkdir -m 775 -p $output/intFiles/$sample/star2passGG`;	
	    my $mrl = $minReadLength - 1;
	    my $ran_sgg2 = 0;
	    my $sgg2j = '';
	    if(!-e "$output/progress/$pre\_$uID\_STAR_GG2_$sample.done" || $ran_star1p){
		sleep(3);
		chdir "$output/intFiles/$sample/star2passGG";
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_GG2_$sample", job_hold => "$star1pj", cpu => "12", mem => "240", cluster_out => "$output/progress/$pre\_$uID\_STAR_GG2_$sample.log");
                #my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_GG2_$sample", job_hold => "$star1pj", cpu => "12", mem => "$STAR_MAX_MEM", cluster_out => "$output/progress/$pre\_$uID\_STAR_GG2_$sample.log");
my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $STAR/STAR --runMode genomeGenerate --genomeDir $output/intFiles/$sample/star2passGG --genomeFastaFiles $REF_SEQ --sjdbFileChrStartEnd $output/intFiles/$sample/$sample\_STAR_1PASS_SJ.out.tab --sjdbOverhang $mrl --runThreadN 12`;
		`/bin/touch $output/progress/$pre\_$uID\_STAR_GG2_$sample.done`;
		$sgg2j = "$pre\_$uID\_STAR_GG2_$sample";
		$ran_sgg2 = 1;
		chdir $curDir;
	    }
	    
	    my $ran_star2p = 0;
	    my $star2pj = '';
	    if(!-e "$output/progress/$pre\_$uID\_STAR_2PASS_$sample.done" || $ran_sgg2){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_2PASS_$sample", job_hold => "$sgg2j", cpu => "12", mem => "$STAR_MAX_MEM", cluster_out => "$output/progress/$pre\_$uID\_STAR_2PASS_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $STAR/STAR --genomeDir $output/intFiles/$sample/star2passGG --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_2PASS_ $star_outSAMstrandField --outFilterIntronMotifs RemoveNoncanonical --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
		`/bin/touch $output/progress/$pre\_$uID\_STAR_2PASS_$sample.done`;
		$star2pj = "$pre\_$uID\_STAR_2PASS_$sample";
		$ran_star2p = 1;
	    }	    
	    
	    if(!-e "$output/progress/$pre\_$uID\_SP_$sample.done" || $ran_star2p){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SP_$sample", job_hold => "$star2pj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_SP_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam`;
		`/bin/touch $output/progress/$pre\_$uID\_SP_$sample.done`;
		$spj = "$pre\_$uID\_SP_$sample";
		$ran_sp = 1;
	    }
	    $starOut = "$output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.sam";
	}

	my $ran_staraddrg = 0;
	my $staraddrgj = '';
	if(!-e "$output/progress/$pre\_$uID\_STAR_AORRG_$sample.done" || $ran_sp){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_AORRG_$sample", job_hold => "$spj", cpu => "3", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_STAR_AORRG_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$starOut O=$output/gene/alignments/$pre\_$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true RGID=$sample\_1 RGLB=_1 RGPL=Illumina RGPU=$sample\_1 RGSM=$sample`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_AORRG_$sample.done`;
            `ln -s $output/gene/alignments/$pre\_$sample\.bai $output/gene/alignments/$pre\_$sample\.bam.bai`;
	    $staraddrgj = "$pre\_$uID\_STAR_AORRG_$sample";
	    $ran_staraddrg = 1;
            push @staraddrg_jids, $staraddrgj;
            push @syncJobs, $staraddrgj;
	}

        next if($alignment_only);

	if($cufflinks){
	    `/bin/mkdir -m 775 -p $output/cufflinks`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/$sample`;
	    if(!-e "$output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.done" || $ran_staraddrg){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUFFLINKS_STAR_$sample", job_hold => "$staraddrgj", cpu => "5", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $CUFFLINKS/cufflinks -q -p 12 --no-update-check $cufflinks_lib_type -N -G $GTF -o $output/cufflinks/$sample $output/gene/alignments/$pre\_$sample\.bam`;
		`/bin/touch $output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.done`;
		push @syncJobs, "$pre\_$uID\_CUFFLINKS_STAR_$sample";
	    }
	}

	if($htseq || $dexseq){
	    my $ran_starqns = 0;
	    my $starqnsj = '';
	    if(!-e "$output/progress/$pre\_$uID\_QNS_STAR_$sample.done" || $ran_sp){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QNS_STAR_$sample", job_hold => "$spj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_QNS_STAR_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles I=$starOut O=$starOut\_queryname_sorted.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID USE_THREADING=true`;
		`/bin/touch $output/progress/$pre\_$uID\_QNS_STAR_$sample.done`;
		$starqnsj = "$pre\_$uID\_QNS_STAR_$sample";
		$ran_starqns = 1;
	    }	    
	
	    if($htseq){
		`/bin/mkdir -m 775 -p $output/gene/counts_gene`;	   
		if(!-e "$output/progress/$pre\_$uID\_HT_STAR_$sample.done" || $ran_starqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HT_STAR_$sample", job_hold => "$starqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_HT_STAR_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $PYTHON/python $HTSEQ/htseq-count -m intersection-strict -s $htseq_stranded -t exon $starOut\_queryname_sorted.sam $GTF > $output/gene/counts_gene/$sample.htseq_count"`;
		    `/bin/touch $output/progress/$pre\_$uID\_HT_STAR_$sample.done`;
		    push @starhtseq_jids, "$pre\_$uID\_HT_STAR_$sample";
		    $ran_starhtseq = 1;
		}
	    }
	
	    if($dexseq){
		`/bin/mkdir -m 775 -p $output/exon/counts_exon`;
     		if(!-e "$output/progress/$pre\_$uID\_DEX_STAR_$sample.done" || $ran_starqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEX_STAR_$sample", job_hold => "$starqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DEX_STAR_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $starOut\_queryname_sorted.sam $output/exon/counts_exon/$sample.dexseq_count`;
		    `/bin/touch $output/progress/$pre\_$uID\_DEX_STAR_$sample.done`;
		    push @stardexseq_jids, "$pre\_$uID\_DEX_STAR_$sample";
		    $ran_stardexseq = 1;
		}
	    }
	}

	my @starcrm = ();
	if(!-e "$output/progress/$pre\_$uID\_STAR_CRM_$sample.done" || $ran_staraddrg){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_CRM_$sample", job_hold => "$staraddrgj", cpu => "1", mem => "3", cluster_out => "$output/progress/$pre\_$uID\_STAR_CRM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/gene/alignments/$pre\_$sample\.bam O=$output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/crm/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=$picard_strand_specificity METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_CRM_$sample.done`;
	    push @starcrm_jids, "$pre\_$uID\_STAR_CRM_$sample";
	    $ran_starcrm = 1;
	}
	push @crm, "-metrics $output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	my @starasm = ();
	if(!-e "$output/progress/$pre\_$uID\_STAR_ASM_$sample.done" || $ran_staraddrg){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_ASM_$sample", job_hold => "$staraddrgj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_STAR_ASM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/gene/alignments/$pre\_$sample\.bam OUTPUT=$output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_ASM_$sample.done`;
	    push @starasm_jids, "$pre\_$uID\_STAR_ASM_$sample";
	    $ran_starasm = 1;
	}
	push @asm, "-metrics $output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt";

        if(exists $ism_samples{$sample}){
           if(!-e "$output/progress/$pre\_$uID\_STAR_ISM_$sample.done" || $ran_staraddrg){
               sleep(3);
               my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_ISM_$sample", job_hold => "$staraddrgj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_STAR_ISM_$sample.log");
               my $standardParams = Schedule::queuing(%stdParams);
               `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE INPUT=$output/gene/alignments/$pre\_$sample\.bam OUTPUT=$output/intFiles/$pre\_$sample\_InsertSizeMetrics.txt HISTOGRAM_FILE=$output/intFiles/$pre\_$sample\_InsertSizeHistogram.txt`;
               `/bin/touch $output/progress/$pre\_$uID\_STAR_ISM_$sample.done`;
               push @starism_jids, "$pre\_$uID\_STAR_ISM_$sample";
               $ran_starism = 1;
           }
           push @ism, "-metrics $output/intFiles/$pre\_$sample\_InsertSizeMetrics.txt";
        }
        if(!-e "$output/progress/$pre\_$uID\_RSEQC_STAR_$sample.done" || $ran_staraddrg){
            sleep(3);
            # only allow a limited number of rseqc jobs to run at the same time, so added additional job holding
            my $job_count = scalar @rseqc_jids;
            my $batch_job_hold = "";
            if($job_count >= $rseqc_batch_size)
            {
                $batch_job_hold = $rseqc_jids[$job_count - $rseqc_batch_size];
            }
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_STAR_$sample", job_hold => "$staraddrgj,$batch_job_hold", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_STAR_$sample.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl -I $Bin/lib $Bin/qc/rseqc.pl -pre $pre -config $config -bam $output/gene/alignments/$pre\_$sample\.bam -bed $QC_BED -sample $sample -intdir $output/intFiles/$sample -outdir $output/metrics/images -progdir $output/progress -scheduler $scheduler -readlen $sampReadLength -layout $samp_pair{$sample} -sync`;
            `/bin/touch $output/progress/$pre\_$uID\_RSEQC_STAR_$sample.done`;
            push @rseqc_jids, "$pre\_$uID\_RSEQC_STAR_$sample";
            push @qcpdf_jids, "$pre\_$uID\_RSEQC_STAR_$sample";
            $ran_rseqc = 1;
        }

    }

    if($detectFusions){
	`/bin/mkdir -m 775 -p $output/fusion`;
	if($species =~ /hg19|human|hybrid/i){	
	    my @fusions = ();
	    my $r1_files = join(" ", @R1);
	    my $r2_files = join(" ", @R2);
	    my $ran_zcat2 = 0;
	    my @zcat2_jids = ();
	    my @fusion_jids = ();
	    my $ran_fusion = 0;
	    if(!-e "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R1.done"){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT2_$sample\_v2_R1", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r1_files ">$output/intFiles/$sample/$sample\_v2_R1.fastq"`;
		`/bin/touch $output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R1.done`;
		push @zcat2_jids, "$pre\_$uID\_ZCAT2_$sample\_v2_R1";
		$ran_zcat2 = 1;
		sleep(3);
	    }

	    if($samp_pair{$sample} eq "PE"){
		if(!-e "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.done"){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT2_$sample\_v2_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r2_files ">$output/intFiles/$sample/$sample\_v2_R2.fastq"`;
		    `/bin/touch $output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.done`;
		    push @zcat2_jids, "$pre\_$uID\_ZCAT2_$sample\_v2_R2";
		    $ran_zcat2 = 1;
		    sleep(3);
		}
	    }	

	    my $zcat2j = join(",", @zcat2_jids);
	    if($chimerascan && $species =~ /human|hg19/i){   ## do not run for hybrid genomes until we debug hanging issue
		### NOTE: CHIMERASCAN ONLY WORKS FOR PE READS
		if($samp_pair{$sample} eq "PE"){
		    my $ran_chimerascan = 0;
		    my $chimerascanj = '';
		    `/bin/mkdir -m 775 -p $output/fusion/chimerascan`;
		    `/bin/mkdir -m 775 -p $output/fusion/chimerascan/$sample`;
		
		    ### NOTE: CHIMERASCAN FAILS WHEN A READ PAIR IS OF DIFFERENT LENGTHS
		    ###       e.g. WHEN WE CLIP AND TRIM VARIABLE LENGTHS FROM
		    ###       SO HAVE TO USE UNPROCESSED READS		    
		    if(!-e "$output/progress/$pre\_$uID\_CHIMERASCAN_$sample.done" || $ran_zcat2){
			if(-d "$output/fusion/chimerascan/$sample/"){
			    `/bin/rm -rf $output/fusion/chimerascan/$sample/`;
			}
			
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CHIMERASCAN_$sample", job_hold => "$zcat2j", cpu => "6", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_CHIMERASCAN_$sample.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $CHIMERASCAN/chimerascan_run.py -p 6 --quals solexa --multihits=10 --filter-false-pos=$CHIMERASCAN_FP_FILTER $CHIMERASCAN_INDEX $output/intFiles/$sample/$sample\_v2_R1.fastq $output/intFiles/$sample/$sample\_v2_R2.fastq $output/fusion/chimerascan/$sample/`;
			`/bin/touch $output/progress/$pre\_$uID\_CHIMERASCAN_$sample.done`;
			$chimerascanj = "$pre\_$uID\_CHIMERASCAN_$sample";
			push @fusion_jids, "$pre\_$uID\_CHIMERASCAN_$sample";
			$ran_chimerascan = 1;
			$ran_fusion = 1;
		    }

		    if(!-e "$output/progress/$pre\_$uID\_CHIMERASCAN_$sample\_CLEANUP.done" || $ran_chimerascan){
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CHIMERASCAN_$sample\_CLEANUP", job_hold => "$chimerascanj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CHIMERASCAN_$sample\_CLEANUP.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/fusion/chimerascan/$sample/log $output/fusion/chimerascan/$sample/tmp`;
			`/bin/touch $output/progress/$pre\_$uID\_CHIMERASCAN_$sample\_CLEANUP.done`;
		    }
		    push @fusions, "--chimerascan $output/fusion/chimerascan/$sample/chimeras.bedpe";
		}
	    }

	    if($star_fusion){
		my $ran_star_fusion = 0;
		my $star_fusionj = '';
		my $inReads = "--left_fq $output/intFiles/$sample/$sample\_v2_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " --right_fq $output/intFiles/$sample/$sample\_v2_R2.fastq";
		}

		if(!-e "$output/progress/$pre\_$uID\_STAR_FUSION_$sample.done" || $ran_zcat2){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_FUSION_$sample", job_hold => "$zcat2j", cpu => "6", mem => "60", cluster_out => "$output/progress/$pre\_$uID\_STAR_FUSION_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 6 --outFileNamePrefix $output/fusion/star/$sample/$sample\_STAR_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --chimSegmentMin 20`;
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $STAR_FUSION/STAR-Fusion --genome_lib_dir $STAR_FUSION_GENOME_LIB $inReads --output_dir $output/fusion/star_fusion/$sample`;

		    `/bin/touch $output/progress/$pre\_$uID\_STAR_FUSION_$sample.done`;
		    $star_fusionj = "$pre\_$uID\_STAR_FUSION_$sample";
		    push @fusion_jids, "$pre\_$uID\_STAR_FUSION_$sample";
		    $ran_star_fusion = 1;
		    $ran_fusion = 1;
		}
		 
		if(!-e "$output/progress/$pre\_$uID\_STAR_FUSION_$sample\_CLEANUP.done" || $ran_star_fusion){
		    my %stdParams = (scheduler => "$scheduler", job_name => "pre\_$uID\_STAR_FUSION_$sample\_CLEANUP", job_hold => "$star_fusionj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_STAR_FUSION_$sample\_CLEANUP.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/fusion/star_fusion/$sample/_STARpass1/ $output/fusion/star_fusion/$sample/_STARgenome/ $output/fusion/star_fusion/$sample/_STARtmp/ $output/fusion/star_fusion/$sample/star-fusion.predict.intermediates_dir/ $output/fusion/star_fusion/$sample/star-fusion.filter.intermediates_dir/`;
		    `/bin/touch $output/progress/$pre\_$uID\_STAR_FUSION_$sample\_CLEANUP.done`;
		}
		push @fusions, "--star $output/fusion/star_fusion/$sample/star-fusion.fusion_candidates.final.abridged";
	    }

	    if($mapsplice && $species =~ /human|hg19/i){ ## do not run for hybrid genomes until we debug
		my $ran_mapsplice = 0;
		my $mapsplicej = '';
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice`;
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice/$sample`;

		my $inReads = "-1 $output/intFiles/$sample/$sample\_v2_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " -2 $output/intFiles/$sample/$sample\_v2_R2.fastq"
		}

		if(!-e "$output/progress/$pre\_$uID\_MAPSPLICE_$sample.done" || $ran_zcat2){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MAPSPLICE_$sample", job_hold => "$zcat2j", cpu => "6", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_MAPSPLICE_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $MAPSPLICE/mapsplice.py -p 6 --bam --fusion-non-canonical -c $chrSplits -x $BOWTIE_INDEX -o $output/fusion/mapsplice/$sample $inReads --gene-gtf $GTF`;
		    `/bin/touch $output/progress/$pre\_$uID\_MAPSPLICE_$sample.done`;
		    $mapsplicej = "$pre\_$uID\_MAPSPLICE_$sample";
		    push @fusion_jids, "$pre\_$uID\_MAPSPLICE_$sample";
		    $ran_mapsplice = 1;
		    $ran_fusion = 1;
		}

		if(!-e "$output/progress/$pre\_$uID\_MAPSPLICE_$sample\_CLEANUP.done" || $ran_mapsplice){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MAPSPLICE_$sample\_CLEANUP", job_hold => "$mapsplicej", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_MAPSPLICE_$sample\_CLEANUP.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/fusion/mapsplice/$sample/logs/`;
		    `/bin/touch $output/progress/$pre\_$uID\_MAPSPLICE_$sample\_CLEANUP.done`;
		}		
		push @fusions, "--mapsplice $output/fusion/mapsplice/$sample/fusions_well_annotated.txt";
	    }

	    if($defuse){
		### NOTE: DEFUSE ONLY WORKS FOR PE READS
		my $defuse_config = '';
		my $ran_defuse = 0;
		my $defusej = '';
		if($species =~ /human|hg19|hybrid/i){
		    $defuse_config = "$DEFUSE/scripts/config_homo_sapiens.txt";
		}
		elsif($species =~ /mouse|mm10/i){
		    $defuse_config = "$DEFUSE/scripts/config_mus_musculus.txt";
		}
		
		if($samp_pair{$sample} eq "PE"){
		    `/bin/mkdir -m 775 -p $output/fusion/defuse`;
		    `/bin/mkdir -m 775 -p $output/fusion/defuse/$sample`;
				
		    if(!-e "$output/progress/$pre\_$uID\_DEFUSE_$sample.done" || $ran_zcat2){
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEFUSE_$sample", job_hold => "$zcat2j", cpu => "12", mem => "600", cluster_out => "$output/progress/$pre\_$uID\_DEFUSE_$sample.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $DEFUSE/scripts/defuse.pl --config $defuse_config --output $output/fusion/defuse/$sample --parallel 12 --1fastq $output/intFiles/$sample/$sample\_v2_R1.fastq --2fastq $output/intFiles/$sample/$sample\_v2_R2.fastq`;
			`/bin/touch $output/progress/$pre\_$uID\_DEFUSE_$sample.done`;
			$defusej = "$pre\_$uID\_DEFUSE_$sample";
			push @fusion_jids, "$pre\_$uID\_DEFUSE_$sample";
			$ran_defuse = 1;
			$ran_fusion = 1;
		    }

		    if(!-e "$output/progress/$pre\_$uID\_DEFUSE_$sample\_CLEANUP.done" || $ran_defuse){
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEFUSE_$sample\_CLEANUP", job_hold => "$defusej", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_DEFUSE_$sample\_CLEANUP.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/fusion/defuse/$sample/jobs/ $output/fusion/defuse/$sample/log/`;
			`/bin/touch $output/progress/$pre\_$uID\_DEFUSE_$sample\_CLEANUP.done`;
		    }
		    push @fusions, "--defuse $output/fusion/defuse/$sample/results.filtered.tsv";
		}
	    }

	    if($fusioncatcher){
		`/bin/mkdir -m 775 -p $output/fusion/fusioncatcher`;
		`/bin/mkdir -m 775 -p $output/fusion/fusioncatcher/$sample`;

		my $inReads = "-i $output/intFiles/$sample/$sample\_v2_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= ",$output/intFiles/$sample/$sample\_v2_R2.fastq"
		}

		if(!-e "$output/progress/$pre\_$uID\_FC_$sample.done" || $ran_zcat2){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FC_$sample", job_hold => "$zcat2j", cpu => "6", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_FC_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $FUSIONCATCHER/bin/fusioncatcher -d $FUSIONCATCHER/data/ensembl_v77d $inReads -o $output/fusion/fusioncatcher/$sample -p 6 --skip-update-check --config=$FUSIONCATCHER/bin/configuration.cfg`;
		    `/bin/touch $output/progress/$pre\_$uID\_FC_$sample.done`;
		    push @fusion_jids, "$pre\_$uID\_FC_$sample";
		    $ran_fusion = 1;
		}
		push @fusions, "--fusioncatcher $output/fusion/fusioncatcher/$sample/final-list_candidate-fusion-genes.txt";
	    }

	    my $mergeFusions = join(" ", @fusions);
	    my $fusionj = join(",", @fusion_jids);
	    my $ran_merge_fusion = 0;
	    my $merge_fusionj = '';
	    if(!-e "$output/progress/$pre\_$uID\_MERGE_FUSION_$sample.done" || $ran_fusion){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_FUSION_$sample", job_hold => "$fusionj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_FUSION_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/MergeFusion $mergeFusions --out $output/fusion/$pre\_merged_fusions_$sample\.txt --normalize_gene $Bin/data/human/hugo_data_073013.tsv`;
		`/bin/touch $output/progress/$pre\_$uID\_MERGE_FUSION_$sample.done`;
		$merge_fusionj = "$pre\_$uID\_MERGE_FUSION_$sample";
		$ran_merge_fusion = 1;
	    }

            if(!-e "$output/progress/$pre\_$uID\_RANK_FUSION_$sample.done" || $ran_merge_fusion){
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RANK_FUSION_$sample", job_hold => "$merge_fusionj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_RANK_FUSION_$sample.log");
                my $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/FusionMetaCaller.R $output/fusion/$pre\_merged_fusions_$sample\.txt $output/fusion/$pre\_merged_fusions_$sample\_ranked\.txt`;
                `/bin/touch $output/progress/$pre\_$uID\_RANK_FUSION_$sample.done`;
		push @syncJobs, "$pre\_$uID\_RANK_FUSION_$sample";
            }

	}
	else{
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSKIPPING FUSION CALLING FOR SAMPLE $sample BECAUSE IT ISN'T SUPPORTED FOR $species; CURRENTLY ONLY SUPPORT FOR hg19";
	}
    }

    if($kallisto || $rsem || $express){
	my $r1_gz_files1 = join(" ", @R1);
	my $r2_gz_files1 = join(" ", @R2);
	my $r1_gz_files2 = join(",", @R1);
	my $r2_gz_files2 = join(",", @R2);
	
	if($express){
	    my $inReads = "-1 $r1_gz_files2";	
	    if($samp_pair{$sample} eq "PE"){
		$inReads .= " -2 $r2_gz_files2";
	    }
	    
	    my $ran_bowtie2 = 0;
	    my $bowtie2j = '';
	    if(!-e "$output/progress/$pre\_$uID\_BOWTIE2_$sample.done"){
		`/bin/mkdir -m 775 -p $output/intFiles/bowtie2`;
		`/bin/mkdir -m 775 -p $output/intFiles/bowtie2/$sample`;
		### express-recommended bowtie2 settings
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_BOWTIE2_$sample", cpu => "24", mem => "100", cluster_out => "$output/progress/$pre\_$uID\_BOWTIE2_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $BOWTIE2/bowtie2 --all --maxins 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed --threads 24 -x $TRANS_INDEX_DEDUP $inReads -S $output/intFiles/bowtie2/$sample/$sample\_bowtie2.sam`;
		`/bin/touch $output/progress/$pre\_$uID\_BOWTIE2_$sample.done`;
		$bowtie2j = "$pre\_$uID\_BOWTIE2_$sample";
		$ran_bowtie2 = 1;
	    }
	    
	    if(!-e "$output/progress/$pre\_$uID\_EXPRESS_$sample.done" || $ran_bowtie2){
		sleep(3);
		`/bin/mkdir -m 775 -p $output/transcript/express`;
		`/bin/mkdir -m 775 -p $output/transcript/express/counts_trans`;
		`/bin/mkdir -m 775 -p $output/transcript/express/counts_trans/$sample`;
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_EXPRESS_$sample", job_hold => "$bowtie2j", cpu => "5", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_EXPRESS_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $EXPRESS/express --output-dir $output/transcript/express/counts_trans/$sample --no-update-check $TRANS_FASTA_DEDUP $output/intFiles/bowtie2/$sample/$sample\_bowtie2.sam`;
		`/bin/touch $output/progress/$pre\_$uID\_EXPRESS_$sample.done`;
		push @syncJobs, "$pre\_$uID\_EXPRESS_$sample";
	    }
	}
    
	if($kallisto || $rsem){
	    my $ran_zcat3 = 0;
	    my @zcat3_jids = ();
	    my $rsem_mode = '';
	    my $kallisto_mode = '--single -l 250';
	    my $kinReads = "$output/intFiles/$sample/$sample\_v3_R1.fastq";
	    if(!-e "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.done"){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT3_$sample\_v3_R1", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r1_gz_files1 ">$output/intFiles/$sample/$sample\_v3_R1.fastq"`;
		`/bin/touch $output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.done`;
		push @zcat3_jids, "$pre\_$uID\_ZCAT3_$sample\_v3_R1";
		$ran_zcat3 = 1;
	    }
	    
	    if($samp_pair{$sample} eq "PE"){
		$rsem_mode = '--paired-end';
		$kallisto_mode = '';
		$kinReads .= " $output/intFiles/$sample/$sample\_v3_R2.fastq";
		if(!-e "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.done"){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT3_$sample\_v3_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/zcat $r2_gz_files1 ">$output/intFiles/$sample/$sample\_v3_R2.fastq"`;
		    `/bin/touch $output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.done`;
		    push @zcat3_jids, "$pre\_$uID\_ZCAT3_$sample\_v3_R2";
		    $ran_zcat3 = 1;
		}
	    }
	    
	    my $zcat3j = join(",", @zcat3_jids);

	    if($kallisto){
		if(!-e "$output/progress/$pre\_$uID\_KALLISTO_$sample.done" || $ran_zcat3){
		    sleep(3);
		    `/bin/mkdir -m 775 -p $output/transcript/kallisto`;
		    `/bin/mkdir -m 775 -p $output/transcript/kallisto/counts_trans`;
		    `/bin/mkdir -m 775 -p $output/transcript/kallisto/counts_trans/$sample`;
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_KALLISTO_$sample", job_hold => "$zcat3j", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_KALLISTO_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $KALLISTO/kallisto quant -i $KALLISTO_INDEX -o $output/transcript/kallisto/counts_trans/$sample -b 100 $kallisto_mode $kinReads`;
		    `/bin/touch $output/progress/$pre\_$uID\_KALLISTO_$sample.done`;
		    push @kallisto_jids, "$pre\_$uID\_KALLISTO_$sample";
		    $ran_kallisto = 1;
		}
	    }

	    if($rsem){
		my $ran_rsem_exp = 0;
		my $rsemj = '';
		if(!-e "$output/progress/$pre\_$uID\_RSEM_$sample.done" || $ran_zcat3){
		    sleep(3);
		    `/bin/mkdir -m 775 -p $output/transcript/rsem`;
		    `/bin/mkdir -m 775 -p $output/transcript/rsem/counts_trans`;
		    `/bin/mkdir -m 775 -p $output/transcript/rsem/counts_trans/$sample`;
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEM_$sample", job_hold => "$zcat3j", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_RSEM_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    #`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $RSEM/rsem-calculate-expression -p 8 $rsem_mode --star --star-path $STAR --estimate-rspd --append-names --output-genome-bam $kinReads $RSEM_DB $output/transcript/rsem/counts_trans/$sample/$sample\_RSEM`;
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $RSEM/rsem-calculate-expression -p 8 $rsem_mode --star --star-path $STAR --estimate-rspd --append-names --output-genome-bam $kinReads $RSEM_DB $output/transcript/rsem/counts_trans/$sample/$sample\_RSEM`;
		    `/bin/touch $output/progress/$pre\_$uID\_RSEM_$sample.done`;
		    $rsemj = "$pre\_$uID\_RSEM_$sample";
		    push @rce_jids, "$pre\_$uID\_RSEM_$sample";
		    $ran_rsem_exp = 1;
		    $ran_rceg = 1;
		}
		
		`/bin/mkdir -m 775 -p $output/metrics/rsem`;
		if(!-e "$output/progress/$pre\_$uID\_RSEM_PLOT_$sample.done" || $ran_rsem_exp){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEM_PLOT_$sample", job_hold => "$rsemj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEM_PLOT_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $RSEM/rsem-plot-model $output/transcript/rsem/counts_trans/$sample/$sample\_RSEM $output/metrics/rsem/$pre\_$sample\_RSEM.pdf`;
		    `/bin/touch $output/progress/$pre\_$uID\_RSEM_PLOT_$sample.done`;
		    push @syncJobs, "$pre\_$uID\_RSEM_PLOT_$sample";
		}		
	    }
	}
    }
}

if($alignment_only) {
    &rsyncJob();
    exit(0);
}


my $ran_shmatrix = 0;
my $shmatrixj = '';
my $ran_deseq = 0;
if($star){
    if($htseq){
	my $starhtseqj = join(",", @starhtseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.done" || $ran_starhtseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_HTSEQ_STAR", job_hold => "$starhtseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/rnaseq_count_matrix.py $output/gene/counts_gene .htseq_count $output/gene/counts_gene/$pre\_htseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.done`;
	    $shmatrixj = "$pre\_$uID\_MATRIX_HTSEQ_STAR";
	    push @syncJobs, "$pre\_$uID\_MATRIX_HTSEQ_STAR";
	    $ran_shmatrix = 1;
	}
    }
    if($dexseq){
	my $stardexseqj = join(",", @stardexseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_DEX_STAR.done" || $ran_stardexseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_DEX_STAR", job_hold => "$stardexseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_DEX_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/rnaseq_count_matrix.py $output/exon/counts_exon .dexseq_count $output/exon/counts_exon/$pre\_dexseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_DEX_STAR.done`;
	    push @syncJobs, "$pre\_$uID\_MATRIX_DEX_STAR";
	}
    }

    my $crmfiles = join(" ", @crm);
    my $scrmj = join(",", @starcrm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_STAR.done" || $ran_starcrm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_STAR", job_hold => "$scrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_STAR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $crmfiles ">$output/metrics/$pre\_CollectRnaSeqMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_STAR.done`;
	push @syncJobs, "$pre\_$uID\_MERGE_CRM_STAR";
        push @qcplot_jids, "$pre\_$uID\_MERGE_CRM_STAR";
        $ran_picard_metrics = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_HIST_STAR.done" || $ran_starcrm){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_HIST_STAR", job_hold => "$scrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_HIST_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/mergeCollectRnaSeqHistograms.py $output/intFiles _CollectRnaSeqMetrics.txt $output/metrics/$pre\_CollectRnaSeqHistograms.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_HIST_STAR.done`;
        push @syncJobs, "$pre\_$uID\_MERGE_CRM_HIST_STAR";
        push @qcplot_jids, "$pre\_$uID\_MERGE_CRM_HIST_STAR";
        $ran_picard_metrics = 1;
    }

    my $asmfiles = join(" ", @asm);
    my $sasmj = join(",", @starasm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM_STAR.done" || $ran_starasm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ASM_STAR", job_hold => "$sasmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ASM_STAR.log");
	my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $asmfiles ">$output/metrics/$pre\_AlignmentSummaryMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM_STAR.done`;
	push @syncJobs, "$pre\_$uID\_MERGE_ASM_STAR";
        push @qcplot_jids, "$pre\_$uID\_MERGE_ASM_STAR";
        $ran_picard_metrics = 1;
    }

    my $ismfiles = join(" ", @ism);
    my $sismj = join(",", @starism_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ISM_STAR.done" || $ran_starism){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ISM_STAR", job_hold => "$sismj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ISM_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $ismfiles ">$output/metrics/$pre\_InsertSizeMetrics.txt"`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_ISM_STAR.done`;
        $ran_starism_merge = 1;
        push @syncJobs, "$pre\_$uID\_MERGE_ISM_STAR";
        push @qcplot_jids, "$pre\_$uID\_MERGE_ISM_STAR";
    }

    if (!-e "$output/progress/$pre\_$uID\_GEO_STAR.done" || $ran_starism_merge){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GEO_STAR", job_hold => "$pre\_$uID\_MERGE_ISM_STAR", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GEO_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my $prepGeoCmd = "$Bin/prepGEO.pl -map $map -pre $pre -species $species -config $config -strand $strand -request $request";
        if($htseq){
            $prepGeoCmd .= " -htseq";
        }
        if($deseq){
            $prepGeoCmd .= " -deseq";
        }
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $prepGeoCmd`; 
        `/bin/touch $output/progress/$pre\_$uID\_GEO_STAR.done`;
        push @syncJobs, "$pre\_$uID\_GEO_STAR";
    }

    my $sargj = join(",",@staraddrg_jids);


    ##### THE GENE_BODY_COVERAGE MODULE TAKES >24 HOURS, SO NOT SURE IF IT'S WORTH EVER KEEPING IT,
    ##### BUT AT LEAST LEAVING IT OUT FOR NOW

    #if(!-e "$output/progress/$pre\_$uID\_RSEQC_GBC_STAR.done" || $ran_star){
    #    sleep(3);
    #    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_GBC_STAR", job_hold => "$sargj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_GBC_STAR.log");
    #    my $standardParams = Schedule::queuing(%stdParams);
    #    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/geneBody_coverage.py -r $QC_BED -i $output/gene/alignments/ -o $output/intFiles/$pre\_RSEQC_STAR_geneBody_coverage`;
    #    `/bin/touch $output/progress/$pre\_$uID\_RSEQC_GBC_STAR.done`;

        ## move only the output pdf file to delivered results
    #    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_GBC_STAR_MV", job_hold => "$pre\_$uID\_RSEQC_GBC_STAR", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_GBC_STAR_MV.log");
    #    my $standardParams = Schedule::queuing(%stdParams);
    #    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams mv $output/intFiles/$pre\_RSEQC_STAR_geneBody_coverage*.pdf $output/metrics/images`;
    #    push @qcpdf_jids, "$pre\_$uID\_RSEQC_GBC_STAR_MV";
    #}

    ## merge read distribution
    my $rseqcj = join(",",@rseqc_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_RD_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_RD_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_RD_STAR.log");
         my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py read_distribution $output/intFiles _read_distribution.txt $output/metrics/$pre\_rseqc_read_distribution.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_RD_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_RD_STAR";
        $ran_rseqc_merge = 1;
    }

    ## merge clipping profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_CP_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_CP_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_CP_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py clipping_profile $output/intFiles .clipping_profile.xls $output/metrics/$pre\_rseqc_clipping_profiles`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_CP_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_CP_STAR";
        $ran_rseqc_merge = 1;
    }

    ## merge GC content
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_GC_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_GC_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_GC_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py gc_content $output/intFiles GC.xls $output/metrics/$pre\_rseqc_gc_content.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_GC_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_GC_STAR";
        $ran_rseqc_merge = 1;
    }

    ## merge bam stats
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_BS_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_BS_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_BS_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py bam_stats $output/intFiles bam_stat.txt $output/metrics/$pre\_rseqc_bam_stats.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_BS_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_BS_STAR";
        $ran_rseqc_merge = 1;
    }

    ## merge insertion profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_IP_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_IP_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_IP_STAR.log");
         my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py insertion_profile $output/intFiles .insertion_profile.xls $output/metrics/$pre\_rseqc_insertion_profiles`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_IP_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_IP_STAR";
        $ran_rseqc_merge = 1;
    }

    ## merge deletion profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_DP_STAR.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_DP_STAR", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_DP_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py deletion_profile $output/intFiles .deletion_profile.txt $output/metrics/$pre\_rseqc_deletion_profiles.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_DP_STAR.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_DP_STAR";
        $ran_rseqc_merge = 1;
    }

    ## Plot merged RSEQC files
    if(!-e "$output/progress/$pre\_$uID\_QC_PLOT_STAR.done" || $ran_rseqc_merge || $ran_picard_metrics){
        sleep(3);
        my $qcplotj = join(",",@qcplot_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QC_PLOT_STAR", job_hold => "$qcplotj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QC_PLOT_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/qc/plot_metrics.R $output/metrics $output/metrics/images $pre $Bin/lib/R-4.0.0 $Bin`;
        `/bin/touch $output/progress/$pre\_$uID\_QC_PLOT_STAR.done`;
        push @qcpdf_jids, "$pre\_$uID\_QC_PLOT_STAR";
        $ran_plots = 1;
    }

    ## check for missing rseqc files
    if(!-e "$output/progress/$pre\_$uID\_RSEQC_CHECK_STAR.done" || $ran_rseqc_merge){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_CHECK_STAR", job_hold => "$pre\_$uID\_QC_PLOT_STAR", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_CHECK_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/check_rseqc.py $pre $map $output/metrics/images`;
        `/bin/touch $output/progress/$pre\_$uID\_RSEQC_CHECK_STAR.done`;
        push @qcpdf_jids, "$pre\_$uID\_RSEQC_CHECK_STAR";
    }

    ## QCPDF containing plots of merged data
    if(!-e "$output/progress/$pre\_$uID\_QCPDF_STAR.done" || $ran_rseqc_merge || $ran_plots){
        sleep(3);
        my $qcpdfj = join(",",@qcpdf_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QCPDF_STAR", job_hold => "$qcpdfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QCPDF_STAR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -cp .:lib/*:lib/java/* -jar $Bin/lib/java/PDFReport.jar -rf $request -v $svnRev -d $output/metrics -o $output/metrics`;
        `/bin/touch $output/progress/$pre\_$uID\_QCPDF_STAR.done`;
        push @syncJobs, "$pre\_$uID\_QCPDF_STAR";
    }
}

my $ran_thmatrix = 0;
my $thmatrixj = '';
if($tophat){
    if($htseq){
	my $tophathtseqj = join(",", @tophathtseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_TOPHAT.done" || $ran_tophathtseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_HTSEQ_TOPHAT", job_hold => "$tophathtseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/rnaseq_count_matrix.py $output/gene/counts_gene/tophat2 .htseq_count $output/gene/counts_gene/tophat2/$pre\_htseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_HTSEQ_TOPHAT.done`;
	    $thmatrixj = "$pre\_$uID\_MATRIX_HTSEQ_TOPHAT";
	    push @syncJobs, "$pre\_$uID\_MATRIX_HTSEQ_TOPHAT";
	    $ran_thmatrix = 1;
	}
    }

    if($dexseq){
	my $tophatdexseqj = join(",", @tophatdexseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.done" || $ran_tophatdexseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_DEX_TOPHAT", job_hold => "$tophatdexseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/rnaseq_count_matrix.py $output/exon/counts_exon/tophat2 .dexseq_count $output/exon/counts_exon/tophat2/$pre\_dexseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.done`;
	    push @syncJobs, "$pre\_$uID\_MATRIX_DEX_TOPHAT";
	}
    }

    my $crmfiles_tophat = join(" ", @crm_tophat);
    my $tcrmj = join(",", @tophatcrm_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.done" || $ran_tophatcrm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_TOPHAT", job_hold => "$tcrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $crmfiles_tophat ">$output/metrics/tophat2/$pre\_CollectRnaSeqMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.done`;
	push @syncJobs, "$pre\_$uID\_MERGE_CRM_TOPHAT";
        push @qcpdf_jids, "$pre\_$uID\_MERGE_CRM_TOPHAT";
        $ran_picard_metrics = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_HIST_TOPHAT.done" || $ran_starcrm){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_HIST_TOPHAT", job_hold => "$tcrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_HIST_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/mergeCollectRnaSeqHistograms.py $output/intFiles _CollectRnaSeqMetrics.txt $output/metrics`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_HIST_TOPHAT.done`;
        push @syncJobs, "$pre\_$uID\_MERGE_CRM_HIST_TOPHAT";
        push @qcplot_jids, "$pre\_$uID\_MERGE_CRM_HIST_TOPHAT";
        $ran_picard_metrics = 1;
    }

    my $asmfiles_tophat = join(" ", @asm_tophat);
    my $tasmj = join(",", @tophatasm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.done" || $ran_tophatasm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ASM_TOPHAT", job_hold => "$tasmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $asmfiles_tophat ">$output/metrics/tophat2/$pre\_AlignmentSummaryMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.done`;
	push @syncJobs, "$pre\_$uID\_MERGE_ASM_TOPHAT";
        push @qcplot_jids, "$pre\_$uID\_MERGE_ASM_TOPHAT";
        $ran_picard_metrics = 1;
    }

    my $ismfiles_tophat = join(" ", @ism_tophat);
    my $tismj = join(",", @tophatism_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ISM_TOPHAT.done" || $ran_tophatism){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ISM_TOPHAT", job_hold => "$tismj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ISM_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/mergePicardMetrics.pl $ismfiles_tophat ">$output/metrics/tophat2/$pre\_InsertSizeMetrics.txt"`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_ISM_TOPHAT.done`;
        push @syncJobs, "$pre\_$uID\_MERGE_ISM_TOPHAT";
        $ran_tophatism_merge = 1;
    }

    if (!-e "$output/progress/$pre\_$uID\_GEO_TOPHAT.done" || $ran_tophatism_merge){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GEO_TOPHAT", job_hold => "$pre\_$uID\_MERGE_ISM_TOPHAT", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GEO_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my $prepGeoCmd = "$Bin/prepGEO.pl -map $map -pre $pre -species $species -config $config -strand $strand -request $request -aligner tophat";
        if($htseq){
            $prepGeoCmd .= " -htseq";
        }
        if($deseq){
            $prepGeoCmd .= " -deseq";
        }        
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $prepGeoCmd`; 
        `/bin/touch $output/progress/$pre\_$uID\_GEO_TOPHAT.done`;
        push @syncJobs, "$pre\_$uID\_GEO_TOPHAT";
    }

    my $throj = join(",",@thro_jids);
    #if(!-e "$output/progress/$pre\_$uID\_RSEQC_GBC_TOPHAT.done" || $ran_tophat){
    #    sleep(3);
    #    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_GBC_TOPHAT", job_hold => "$throj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_GBC_TOPHAT.log");
    #    my $standardParams = Schedule::queuing(%stdParams);
    #    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/geneBody_coverage.py -r $QC_BED -i $output/gene/alignments/ -o $output/intFiles/$pre\_RSEQC_TOPHAT_geneBody_coverage`;
    #    `/bin/touch $output/progress/$pre\_$uID\_RSEQC_GBC_TOPHAT.done`;

    #    ## move only the output pdf file to delivered results
    #    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_GBC_TOPHAT_MV", job_hold => "$pre\_$uID\_RSEQC_GBC_TOPHAT", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSEQC_GBC_TOPHAT_MV.log");
    #    my $standardParams = Schedule::queuing(%stdParams);
    #    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams mv $output/intFiles/$pre\_RSEQC_TOPHAT_geneBody_coverage*.pdf $output/metrics/images`;
    #    push @qcpdf_jids, "$pre\_$uID\_RSEQC_GBC_TOPHAT_MV";
    #}

    ## merge read distribution
    my $rseqcj = join(",",@rseqc_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_RD_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_RD_TOPHAT", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_RD_TOPHAT.log");
         my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py read_distribution $output/intFiles _read_distribution.txt $output/metrics/$pre\_rseqc_read_distribution.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_RD_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_RD_TOPHAT";
    }

    ## merge clipping profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_CP_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_CP_TOPHAT", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_CP_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py clipping_profile $output/intFiles .clipping_profile.xls $output/metrics/$pre\_rseqc_clipping_profiles`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_CP_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_CP_TOPHAT";
    }

    ## merge GC content
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_GC_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_GC_TOPHAT", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_GC_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py gc_content $output/intFiles GC.xls $output/metrics/$pre\_rseqc_gc_content.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_GC_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_GC_TOPHAT";
    }

    ## merge bam stats
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_BS_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_BS_TOPHAT", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_BS_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py bam_stats $output/intFiles bam_stat.txt $output/metrics/$pre\_rseqc_bam_stats.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_BS_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_BS_TOPHAT";
    }

    ## merge insertion profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_IP_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_IP_TOPHAT", job_hold => "$rseqcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_IP_TOPHAT.log");
         my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py insertion_profile $output/intFiles .insertion_profile.xls $output/metrics/$pre\_rseqc_insertion_profiles`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_IP_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_IP_TOPHAT";
    }

    ## merge deletion profiles
    if(!-e "$output/progress/$pre\_$uID\_MERGE_RSEQC_DP_TOPHAT.done" || $ran_rseqc){
        sleep(3);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_RSEQC_DP_TOPHAT", job_hold => "$rseqcj", DPu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_RSEQC_DP_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/merge_rseqc_stats.py deletion_profile $output/intFiles .deletion_profile.txt $output/metrics/$pre\_rseqc_deletion_profiles.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_RSEQC_DP_TOPHAT.done`;
        push @qcplot_jids, "$pre\_$uID\_MERGE_RSEQC_DP_TOPHAT";
    }

    if(!-e "$output/progress/$pre\_$uID\_QC_PLOT_TOPHAT.done" || $ran_rseqc || $ran_picard_metrics){
        sleep(3);
        my $qcplotj = join(",",@qcplot_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QC_PLOT_TOPHAT", job_hold => "$qcplotj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QC_PLOT_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/qc/plot_metrics.R $output/metrics $output/metrics $pre $Bin/lib/R`;
        `/bin/touch $output/progress/$pre\_$uID\_QC_PLOT_TOPHAT.done`;
        push @qcpdf_jids, "$pre\_$uID\_QC_PLOT_TOPHAT";
        $ran_plots = 1;
    }

    ## QC PDF
    if(!-e "$output/progress/$pre\_$uID\_QCPDF_TOPHAT.done" || $ran_rseqc || $ran_plots){
        sleep(3);
        my $qcpdfj = join(",",@qcpdf_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QCPDF_TOPHAT", job_hold => "$qcpdfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QCPDF_TOPHAT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -jar $Bin/lib/java/PDFReport.jar -rf $request -v $svnRev -d $output/metrics -o $output/metrics`;
        `/bin/touch $output/progress/$pre\_$uID\_QCPDF_TOPHAT.done`;
        push @syncJobs, "$pre\_$uID\_QCPDF_TOPHAT";
    }
 
}


my $ran_kmatrix = 0;
my $kmatrixj = '';
if($kallisto){
    my $kallistoj = join(",", @kallisto_jids);
    if(!-e "$output/progress/$pre\_$uID\_MATRIX_KALLISTO.done" || $ran_kallisto){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_KALLISTO", job_hold => "$kallistoj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_KALLISTO.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/mergeKallistoAbundance.py $output/transcript/kallisto/counts_trans abundance.txt $output/transcript/kallisto/counts_trans/$pre\_kallisto_all_samples.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MATRIX_KALLISTO.done`;
	$kmatrixj = "$pre\_$uID\_MATRIX_KALLISTO";
	push @syncJobs, "$pre\_$uID\_MATRIX_KALLISTO";
	$ran_kmatrix = 1;
    }
}

my $ran_rmatrix = 0;
my $rmatrixj = '';
if($rsem){
    my $rsemj = join(",", @rce_jids);
    if(!-e "$output/progress/$pre\_$uID\_MATRIX_RSEM.done" || $ran_rceg){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_RSEM", job_hold => "$rsemj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_RSEM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/mergeRSEMcounts.py $output/transcript/rsem/counts_trans isoforms.results $output/transcript/rsem/counts_trans/$pre\_rsem_all_samples.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MATRIX_RSEM.done`;
	$rmatrixj = "$pre\_$uID\_MATRIX_RSEM";
	push @syncJobs, "$pre\_$uID\_MATRIX_RSEM";
	$ran_rmatrix = 1;
    }
}

if($clustering_only){
    if($star){
        `/bin/mkdir -m 775 -p $output/gene/clustering`;
        if(!-e "$output/progress/$pre\_$uID\_SAMPLE_CLUSTERING_STAR.done" || $ran_shmatrix){
            sleep(3);
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SAMPLE_CLUSTERING_STAR", job_hold => "$shmatrixj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_SAMPLE_CLUSTERING_STAR.log");
            my $standardParams = Schedule::queuing(%stdParams);

            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/RunDE.R --bin $Bin --counts.file $output/gene/counts_gene/$pre\_htseq_all_samples.txt --Rlibs $Bin/lib/R-4.0.0 --pre $pre --counts.dir $output/gene/counts_gene --clustering.dir $output/gene/clustering --diff.exp FALSE --GSA FALSE`; 

            `/bin/touch $output/progress/$pre\_$uID\_SAMPLE_CLUSTERING_STAR.done`;
            push @syncJobs, "$pre\_$uID\_SAMPLE_CLUSTERING_STAR";
        }
    }

    #### TODO: ADD CLUSTERING FOR TOPHAT, RSEM, KALLISTO AND DELETE OLD VERSION BELOW

}


if($deseq){   
    my $deseq_species = $species;
    if($species =~ /hybrid/i){
        $deseq_species = "human";
    }
 
    if($star){
	`/bin/mkdir -m 775 -p $output/gene/differentialExpression_gene`;
	`/bin/mkdir -m 775 -p $output/gene/clustering`;
	`/bin/mkdir -m 775 -p $output/gene/gsa`;
	if(!-e "$output/progress/$pre\_$uID\_DESeq_STAR.done" || $ran_shmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_STAR", job_hold => "$shmatrixj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_DESeq_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);

            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/RunDE.R --bin $Bin --counts.file $output/gene/counts_gene/$pre\_htseq_all_samples.txt --species $deseq_species --gmt.dir $GMT_DIR --Rlibs $Bin/lib/R-4.0.0 --pre $pre --diff.exp.dir $output/gene/differentialExpression_gene --all.gene.dir $output/gene/all_gene --counts.dir $output/gene/counts_gene --clustering.dir $output/gene/clustering --gsa.dir $output/gene/gsa --java $JAVA/java --javacp .:lib/*:lib/java/* --pdfjar $Bin/lib/java/PDFReport.jar --request $request --svnrev $svnRev`;

	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_STAR.done`;
	    push @syncJobs, "$pre\_$uID\_DESeq_STAR";
            $ran_deseq = 1;
	}
    }

    if($tophat){
	`/bin/mkdir -m 775 -p $output/gene/differentialExpression_gene/tophat2`;
	`/bin/mkdir -m 775 -p $output/gene/clustering/tophat2`;
	`/bin/mkdir -m 775 -p $output/gene/gsa/tophat2`;

	if(!-e "$output/progress/$pre\_$uID\_DESeq_TOPHAT.done" || $ran_thmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_TOPHAT", job_hold => "$thmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DESeq_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    #`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -pre $pre -diff_out $output/gene/differentialExpression_gene/tophat2 -count_out $output/gene/counts_gene/tophat2 -cluster_out $output/gene/clustering/tophat2 -gsa_out $output/gene/gsa/tophat2 -config $config -bin $Bin -species $species -counts $output/gene/counts_gene/tophat2/$pre\_htseq_all_samples.txt -samplekey $samplekey -comparisons $comparisons $reps -Rlibs $Bin/lib/R`;
	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_TOPHAT.done`;
	    push @syncJobs, "$pre\_$uID\_DESeq_TOPHAT";
	}
    }


    if($kallisto){
	`/bin/mkdir -m 775 -p $output/transcript/kallisto/differentialExpression_trans`;
	`/bin/mkdir -m 775 -p $output/transcript/kallisto/clustering`;
	
	if(!-e "$output/progress/$pre\_$uID\_DESeq_KALLISTO.done" || $ran_kmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_KALLISTO", job_hold => "$kmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DESeq_KALLISTO.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    #`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -pre $pre -diff_out $output/transcript/kallisto/differentialExpression_trans -count_out $output/transcript/kallisto/counts_trans -cluster_out $output/transcript/kallisto/clustering -config $config -bin $Bin -species $species -counts $output/transcript/kallisto/counts_trans/$pre\_kallisto_all_samples.txt -samplekey $samplekey -comparisons $comparisons $reps -Rlibs $Bin/lib/R`;
	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_KALLISTO.done`;
	    push @syncJobs, "$pre\_$uID\_DESeq_KALLISTO";
	}
    }

    if($rsem){
	`/bin/mkdir -m 775 -p $output/transcript/rsem/differentialExpression_trans`;
	`/bin/mkdir -m 775 -p $output/transcript/rsem/clustering`;
	
	if(!-e "$output/progress/$pre\_$uID\_DESeq_RSEM.done" || $ran_rmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_RSEM", job_hold => "$rmatrixj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_DESeq_RSEM.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    #`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -pre $pre -diff_out $output/transcript/rsem/differentialExpression_trans -count_out $output/transcript/rsem/counts_trans -cluster_out $output/transcript/rsem/clustering -config $config -bin $Bin -species $species -counts $output/transcript/rsem/counts_trans/$pre\_rsem_all_samples.txt -samplekey $samplekey -comparisons $comparisons $reps -Rlibs $Bin/lib/R`;
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $R/Rscript $Bin/RunDE.R --bin $Bin --counts.file $output/transcript/rsem/counts_trans/$pre\_rsem_all_samples.txt --species $deseq_species --Rlibs $Bin/lib/R-4.0.0 --pre $pre --diff.exp.dir $output/transcript/rsem/differentialExpression_trans --all.gene.dir $output/transcript/rsem/all_trans --counts.dir $output/transcript/rsem/counts_trans --clustering.dir $output/transcript/rsem/clustering --GSA FALSE --java $JAVA/java --javacp .:lib/*:lib/java/* --pdfjar $Bin/lib/java/PDFReport.jar --request $request --svnrev $svnRev`;

	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_RSEM.done`;
	    push @syncJobs, "$pre\_$uID\_DESeq_RSEM";
	}
    }
}
else{
=begin CLUSTER
    if($star && $htseq){
	`/bin/mkdir -m 775 -p $output/gene/clustering`;
	if(!-e "$output/progress/$pre\_$uID\_CLUSTERING_STAR.done" || $ran_shmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTERING_STAR", job_hold => "$shmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CLUSTERING_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -count_out $output/gene/counts_gene -cluster_out $output/gene/clustering -config $config -bin $Bin -counts $output/gene/counts_gene/$pre\_htseq_all_samples.txt -clusterOnly`;
	    `/bin/touch $output/progress/$pre\_$uID\_CLUSTERING_STAR.done`;
	}
    }

    if($tophat && $htseq){
	`/bin/mkdir -m 775 -p $output/gene/clustering/tophat2`;
	if(!-e "$output/progress/pre\_$uID\_CLUSTERING_TOPHAT.done" || $ran_thmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "pre\_$uID\_CLUSTERING_TOPHAT", job_hold => "$thmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/pre\_$uID\_CLUSTERING_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -count_out $output/gene/counts_gene/tophat2 -cluster_out $output/gene/clustering/tophat2 -config $config -bin $Bin -counts $output/gene/counts_gene/tophat2/$pre\_htseq_all_samples.txt -clusterOnly`;	    
	}
    }

    if($kallisto){
	`/bin/mkdir -m 775 -p $output/transcript/kallisto/clustering`;	
	if(!-e "$output/progress/$pre\_$uID\_CLUSTERING_KALLISTO.done" || $ran_kmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTERING_KALLISTO", job_hold => "$kmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CLUSTERING_KALLISTO.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -count_out $output/transcript/kallisto/counts_trans -cluster_out $output/transcript/kallisto/clustering -config $config -bin $Bin -counts $output/transcript/kallisto/counts_trans/$pre\_kallisto_all_samples.txt -clusterOnly`;
	    `/bin/touch $output/progress/$pre\_$uID\_CLUSTERING_KALLISTO.done`;
	}
    }

    if($rsem){
	`/bin/mkdir -m 775 -p $output/transcript/rsem/clustering`;	
	if(!-e "$output/progress/$pre\_$uID\_CLUSTERING_RSEM.done" || $ran_rmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTERING_RSEM", job_hold => "$rmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CLUSTERING_RSEM.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/run_DESeq_wrapper.pl -count_out $output/transcript/rsem/counts_trans -cluster_out $output/transcript/rsem/clustering -config $config -bin $Bin -counts $output/transcript/rsem/counts_trans/$pre\_rsem_all_samples.txt -clusterOnly`;
	    `/bin/touch $output/progress/$pre\_$uID\_CLUSTERING_RSEM.done`;
	}
    }
=end CLUSTER
=cut
}

if($r1adaptor){
    my $cagj = join(",", @cag_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_CAS.done" || $ran_cag){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", job_hold => "$cagj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/mergeCutAdaptStats.py . '*CUTADAPT_STATS.txt' $output/metrics/$pre\_CutAdaptStats.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CAS.done`;
	push @syncJobs, "$pre\_$uID\_MERGE_CAS";
    }
}
close LOG;


&rsyncJob();



sub rsyncJob{
    my $sj = join(",", @syncJobs);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC", job_hold => "$sj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_3.log");
    my $standardParams = Schedule::queuing(%stdParams);
    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1", mail => "$email");
    my $additionalParams = Schedule::additionalParams(%addParams);
    $ENV{'LSB_JOB_REPORT_MAIL'} = 'Y';
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' $curDir $rsync`;
    `/bin/touch $output/progress/$pre\_$uID\_RSYNC.done`;
}


sub verifyConfig{
    my $paths = shift;

    my $test_path = $ENV{'PATH'};

    $singularityenv_prepend_path .= ":/opt/common/CentOS_6/gcc/gcc-4.9.3/bin";

    print LOG "PATH BEFORE VERIFY CONFIG...   $test_path\n\n";

    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
	chomp;
	
	my @conf = split(/\s+/, $_);
        if($conf[0] =~ /singularity/i){
            if(!-e "$conf[1]/singularity"){
                die "CAN'T FIND singularity IN $conf[1] $!";
            }
            $SINGULARITY = $conf[1];
        }
	elsif($conf[0] =~ /picard/i){
	    if(!-e "$conf[1]/picard.jar"){
		die "CAN'T FIND picard.jar IN $conf[1] $!";
	    }
	    $PICARD = $conf[1];
	}
	elsif($conf[0] =~ /cufflinks/i){
	    if(!-e "$conf[1]/cufflinks"){
		die "CAN'T FIND cufflinks IN $conf[1] $!";
	    }
	    $CUFFLINKS = $conf[1];
	}
	elsif($conf[0] =~ /htseq/i){
	    if(!-e "$conf[1]/htseq-count"){
		die "CAN'T FIND htseq-count IN $conf[1] $!";
	    }
	    $HTSEQ = $conf[1];
	}
	elsif($conf[0] =~ /dexseq/i){
	    if(!-e "$conf[1]/dexseq_count.py"){
		die "CAN'T FIND dexseq_count.py IN $conf[1] $!";
	    }
	    $DEXSEQ = $conf[1];
	}
	elsif($conf[0] =~ /chimerascan/i){
	    if(!-e "$conf[1]/chimerascan_run.py"){
		die "CAN'T FIND chimerascan_run.py IN $conf[1] $!";
	    }
	    $CHIMERASCAN = $conf[1];
	}
	elsif($conf[0] =~ /samtools/i){
	    if(!-e "$conf[1]/samtools"){
		die "CAN'T FIND samtools IN $conf[1] $!";
	    }
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /^star$/i){
	    if(!-e "$conf[1]/STAR"){
		die "CAN'T FIND STAR IN $conf[1] $!";
	    }
	    $STAR = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /star_fusion/i){
	    if(!-e "$conf[1]/STAR-Fusion"){
		die "CAN'T FIND STAR-Fusion IN $conf[1] $!";
	    }
	    $STAR_FUSION = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /mapsplice/i){
	    if(!-e "$conf[1]/mapsplice.py"){
		die "CAN'T FIND mapsplice.py IN $conf[1] $!";
	    }
	    $MAPSPLICE = $conf[1];
	}
	elsif($conf[0] =~ /defuse/i){
	    if(!-e "$conf[1]/scripts/defuse.pl"){
		die "CAN'T FIND scripts/defuse.pl IN $conf[1] $!";
	    }
	    $DEFUSE = $conf[1];
	}
	elsif($conf[0] =~ /fusioncatcher/i){
	    if(!-e "$conf[1]/bin/fusioncatcher"){
		die "CAN'T FIND bin/fusioncatcher IN $conf[1] $!";
	    }
	    $FUSIONCATCHER = $conf[1];
	}
	elsif($conf[0] =~ /tophat/i){
	    if(!-e "$conf[1]/tophat2"){
		die "CAN'T FIND tophat2 IN $conf[1] $!";
	    }
	    $TOPHAT = $conf[1];
	}
	elsif($conf[0] =~ /express/i){
	    if(!-e "$conf[1]/express"){
		die "CAN'T FIND express IN $conf[1] $!";
	    }
	    $EXPRESS = $conf[1];
	}
	elsif($conf[0] =~ /kallisto/i){
	    if(!-e "$conf[1]/kallisto"){
		die "CAN'T FIND kallisto IN $conf[1] $!";
	    }
	    $KALLISTO = $conf[1];
	}
	elsif($conf[0] =~ /^bowtie2$/i){
	    ### need bowtie2 in path for tophat to run
	    if(!-e "$conf[1]/bowtie2"){
		die "CAN'T FIND bowtie2 IN $conf[1] $!";
	    }
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	    $BOWTIE2 = $conf[1];
	}
	elsif($conf[0] =~ /^bowtie$/i){
	    ### need bowtie in path for chimerascan to run
	    if(!-e "$conf[1]/bowtie"){
		die "CAN'T FIND bowtie IN $conf[1] $!";
	    }
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /java/i){
	    if(!-e "$conf[1]/java"){
		die "CAN'T FIND java IN $conf[1] $!";
	    }
	    $JAVA = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /perl/i){
	    if(!-e "$conf[1]/perl"){
		die "CAN'T FIND perl IN $conf[1] $!";
	    }
	    $PERL = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	}
 	elsif($conf[0] =~ /python/i){
	    if(!-e "$conf[1]/python"){
		die "CAN'T FIND python IN $conf[1] $!";
	    }
	    $PYTHON = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /^r$/i){
	#    if(!-e "$conf[1]/R"){
#		die "CAN'T FIND R IN $conf[1] $!";
#	    }
	    $R = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /rsem/i){
	    if(!-e "$conf[1]/rsem-calculate-expression"){
		die "CAN'T FIND rsem-calculate-expression IN $conf[1] $!";
	    }
	    $RSEM = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
            $singularityenv_prepend_path .= ":$conf[1]";
	}
	elsif($conf[0] =~ /cutadapt/i){
	    if(!-e "$conf[1]/cutadapt"){
		die "CAN'T FIND cutadapt IN $conf[1] $!";
	    }
	    $CUTADAPT = $conf[1];
	}
    }

    $test_path = $ENV{'PATH'};
    print LOG "PATH AFTER VERIFY CONFIG...   $test_path\n\n";

    my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$Bin/rnaseq_pipeline_singularity_prod.simg");
    $singularityParams = Schedule::singularityParams(%sinParams);
    $singularityBind = Schedule::singularityBind($scheduler);

    $ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
    $ENV{'SINGULARITY_BINDPATH'} = $singularityBind;
    $ENV{'SINGULARITYENV_LD_LIBRARY_PATH'} = "/opt/common/CentOS_6/gcc/gcc-4.9.3/lib64:$ENV{'LD_LIBRARY_PATH'}";


    print LOG "SINGULARITYENV_PREPEND_PATH...\n  $ENV{'SINGULARITYENV_PREPEND_PATH'}\n\n";
    print LOG "SINGULARITY_BINDPATH...\n         $ENV{'SINGULARITY_BINDPATH'}\n\n";
    print LOG "SINGULARITY_LD_LIBRARY_PATH...\n  $ENV{'SINGULARITYENV_LD_LIBRARY_PATH'}\n\n";

    close CONFIG;
}

sub getTime(){

    my @timeData = localtime();

    if($timeData[0] < 10){
	$timeData[0] = "0" . $timeData[0];
    }
    if($timeData[1] < 10){
	$timeData[1] = "0" . $timeData[1];
    }
    if($timeData[2] < 10){
	$timeData[2] = "0" . $timeData[2];
    }
    if($timeData[3] < 10){
	$timeData[3] = "0" . $timeData[3];
    }
    
    my $month = $timeData[4] + 1;
    $timeData[4] += 1;
    if($timeData[4] < 10){
	$timeData[4] = "0" . $timeData[4];
    }
    
    $timeData[5] += 1900;

    return(@timeData);
}

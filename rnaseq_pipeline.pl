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


my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1, $lncrna, $lincrna_BROAD, $output, $strand, $r1adaptor, $r2adaptor, $scheduler, $transcript, $no_replicates);

$pre = 'TEMP';
$output = "results";
my $priority_project = "ngs";
my $priority_group = "Pipeline";

GetOptions ('map=s' => \$map,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'samplekey=s' => \$samplekey,
	    'comparisons=s' => \$comparisons,
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
	    'transcript' => \$transcript,
            'species=s' => \$species,
            'strand=s' => \$strand,
            'lncrna' => \$lncrna,
 	    'output|out|o=s' => \$output,
	    'r1adaptor=s' => \$r1adaptor,
	    'r2adaptor=s' => \$r2adaptor,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
            'lincrna_BROAD' => \$lincrna_BROAD,
            'no_replicates' => \$no_replicates) or exit(1);


if(!$map || !$species || !$strand || !$config || !$scheduler || $help){
    print <<HELP;

    USAGE: rnaseq_pipeline.pl -map MAP -species SPECIES -strand STRAND -config CONFIG -pre PRE -samplekey SAMPLEKEY -comparisons COMPARISONS -scheduler SCHEDULER
	* MAP: file listing sample mapping information for processing (REQUIRED)
	* SPECIES: only hg19, mouse (mm10; default) and human-mouse hybrid (hybrid) currently supported (REQUIRED)
	* STRAND: library strand; valid options are none, forward, reverse (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B (if -deseq, REQUIRED)
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B (if -deseq, REQUIRED)
	* R1ADAPTOR/R2ADAPTOR: if provided, will trim adaptor sequences; NOTE: if provided for only one end, will also assign it to the other end
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless -pass1 specified; tophat2 (-tophat); if no aligner specifed, will default to STAR
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons); fusion callers chimerascan (-chimerascan), rna star (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs, transcript analysis using express and kallisto (-transcript)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* OUTPUT: output results directory (default: results)
        * OPTIONS: lncRNA analysis (-lncrna) runs all analyses based on lncRNA GTF (hg19 only); 
HELP
exit;
}

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($species !~ /human|hg19|mouse|mm9|mm10|hybrid|zebrafish|zv9|dm3|fly/i){
    die "Species must be human (hg19), mouse (mm9/mm10), human-mouse hybrid (hybrid), fly (dm3), or zebrafish (zv9)";
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
if($strand){
    $commandLine .= " -strand $strand";
}
if($samplekey){
    $commandLine .= " -samplekey $samplekey";
}
if($comparisons){
    $commandLine .= " -comparisons $comparisons";
}
if($scheduler){
    $commandLine .= " -scheduler $scheduler";
}
if($r1adaptor){
    $commandLine .= " -r1adaptor $r1adaptor";
}
if($r2adaptor){
    $commandLine .= " -r2adaptor $r2adaptor";
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
if($transcript){
    $commandLine .= " -transcript";
}
if($species){
    $commandLine .= " -species $species";
}
if($lncrna){
    $commandLine .= " -lncrna";
}
if($no_replicates){
    $no_replicates .= " -no_replicates";
}

$commandLine.= " -priority_project $priority_project";
$commandLine .= " -priority_group $priority_group";


my $numArgs = $#ARGV + 1;
foreach my $argnum (0 .. $#ARGV) {
    $commandLine .= " $ARGV[$argnum]";
}


my %mapping_samples = ();
open(MA, "$map") or die "Can't open mapping file $map $!";
while(<MA>){
    chomp;
    
    my @data = split(/\s+/, $_);
    
    $mapping_samples{$data[1]} = 1;
    if(!-d $data[3]){
	die "$data[3] does not exist\n";
    }
}
close MA;

if($comparisons || $samplekey || $deseq){
    my %sample_comparisons = ();
    open(SC, "$comparisons") or die "Can't open comparisons file $comparisons $!";
    while(<SC>){
	chomp;
	
	my @data = split(/\s+/, $_);
	$sample_comparisons{$data[0]} = 1;
	$sample_comparisons{$data[1]} = 1;
    }
    close SC;
    
    my %samplekey_samples = ();
    my %samplekey_conditions = ();
    open(SK, "$samplekey") or die "Can't open key file $samplekey $!";
    while(<SK>){
	chomp;
	
	my @data = split(/\s+/, $_);
	$samplekey_samples{$data[0]} = 1;
	$samplekey_conditions{$data[1]} = 1;
	print "data[0]: $data[0]\tdata[1]: $data[1]\n";
	if(!$mapping_samples{$data[0]} || !$sample_comparisons{$data[1]}){
	    die "either sample $data[0] cannot be found in $map and/or condition $data[1] cannot be found in $comparisons $!";
	}
    }
    close SK;

    foreach my $ms (keys %mapping_samples){
	if(!$samplekey_samples{$ms}){
	    die "sample $ms is in mapping file $map, but not in sample key file $samplekey $!";
	}
    }

    foreach my $sc (keys %sample_comparisons){
	if(!$samplekey_conditions{$sc}){
	    die "condition $sc is in sample comparisons file $comparisons, but not in sample key file $samplekey $!";
	}
    }
}

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

if($species =~ /human|hg19/i){
    $species = 'hg19';
    $REF_SEQ = '/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta';
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $CHIMERASCAN_INDEX = '/ifs/depot/assemblies/H.sapiens/hg19/index/chimerascan/0.4.5a';
    $BOWTIE_INDEX = '/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/1.0.0/hg19_bowtie';
    $BOWTIE2_INDEX = '/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/hg19_bowtie2';
    $chrSplits = '/ifs/depot/assemblies/H.sapiens/hg19/chromosomes';
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_hg19.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__hg19.txt.gz";
    $TRANS_INDEX = '/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/gencode.v18.annotation';
    $TRANS_INDEX_DEDUP = '/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup';
    $TRANS_FASTA_DEDUP = '/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup.fa';
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
    $KALLISTO_INDEX = '/ifs/depot/assemblies/H.sapiens/hg19/index/kallisto/v0.42.1/gencode/v18/gencode.v18.annotation.gtf.fasta.idx';

    if($lncrna){
        $GTF = "$Bin/data/lncipedia.gtf";
        $starDB = '/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.3.0e_r291/LNCipedia';
    } else {
        $GTF = "$Bin/data/gencode.v18.annotation.gtf";
	if($r1adaptor){
	    $starDB = '/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang49';
	}
	else{
	    $starDB = '/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang74';
	}
    }
}
elsif($species =~ /mouse|mm10/i){
    $species = 'mm10';
    $REF_SEQ = '/ifs/depot/assemblies/M.musculus/mm10/mm10.fasta';
    $GTF = "$Bin/data/Mus_musculus.GRCm38.80_canonical_chromosomes.gtf";
    $DEXSEQ_GTF = "$Bin/data/Mus_musculus.GRCm38.80_canonical_chromosomes.dexseq.gtf";
    $geneNameConversion = "$Bin/data/mm10Ensembl80IDToGeneName.txt";
    $BOWTIE_INDEX = '/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/1.1.1/mm10_bowtie';
    $BOWTIE2_INDEX = '/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/mm10_bowtie2';
    $chrSplits = '/ifs/depot/assemblies/M.musculus/mm10/chromosomes';
    $TRANS_INDEX = '/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/transcriptome/ensembl/v80/Mus_musculus.GRCm38.80_canonical_chromosomes';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_mm10.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__mm10.txt.gz";
    $KALLISTO_INDEX = '';

    if($r1adaptor){
	$starDB = '/ifs/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/ensembl/v80/overhang49';
    }
    else{
	$starDB = '/ifs/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/ensembl/v80/overhang74';
    }
}
elsif($species =~ /mm9/i){
    $species = 'mm9';
    $REF_SEQ = '/ifs/depot/assemblies/M.musculus/mm9/mm9.fasta';
    $GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.gtf";
    $DEXSEQ_GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf";
    $geneNameConversion = "$Bin/data/mm9Ensembl67IDToGeneName.txt";
    $BOWTIE_INDEX = '/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/1.0.0/mm9_bowtie';
    $BOWTIE2_INDEX = '/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/mm9_bowtie2';
    $chrSplits = '/ifs/depot/assemblies/M.musculus/mm9/chromosomes';
    $TRANS_INDEX = '/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/transcriptome/ensembl/vTBD/ensembl';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_MM9_assemblies.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__mm9.txt.gz";
    $KALLISTO_INDEX = '';

    if($r1adaptor){
	$starDB = '/ifs/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang49';
    }
    else{
	$starDB = '/ifs/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang74';
    }
}
elsif($species =~ /human-mouse|mouse-human|hybrid/i){
    $species = 'hybrid';
    $REF_SEQ = '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_hg19.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__hg19.txt.gz";

    if($r1adaptor){
	$starDB = '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang49';
    }
    else{
	$starDB = '/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang74';
    }
}
elsif($species =~ /zebrafish|zv9/i){
    $species = 'zv9';
    $REF_SEQ = '';
    $GTF = "$Bin/data/zv9.gtf";
    $starDB = '';
    $chrSplits = '';
    $geneNameConversion = "Bin/data/zv9EnsemblIDtoGeneName.txt";
}
elsif($species =~ /fly|dm3/i){
    $species = 'dm3';
    $REF_SEQ = '/ifs/depot/assemblies/D.melanogaster/dm3/dm3.fasta';
    $GTF = "$Bin/data/dm3.flybase_more150bp_CollapseGenes_20140925.gtf";
    $DEXSEQ_GTF = "";
    $geneNameConversion = "";
    $BOWTIE_INDEX = '/ifs/depot/assemblies/D.melanogaster/dm3/index/bowtie/1.1.1/dm3_bowtie';
    $BOWTIE2_INDEX = '/ifs/depot/assemblies/D.melanogaster/dm3/index/bowtie/2.2.4/dm3_bowtie2'; 
    $chrSplits = '/ifs/depot/assemblies/D.melanogaster/dm3/chromosomes';
    $TRANS_INDEX = '';
    $TRANS_INDEX_DEDUP = '';
    $TRANS_FASTA_DEDUP = '';
    $RIBOSOMAL_INTERVALS = '';
    $REF_FLAT = '';
    $KALLISTO_INDEX = '';

    if($r1adaptor){
        $starDB = '/ifs/depot/assemblies/D.melanogaster/dm3/index/star/2.4.1d/flybase/custom20140925/overhang49'; 
    }
    else{
        $starDB = '/ifs/depot/assemblies/D.melanogaster/dm3/index/star/2.4.1d/flybase/custom20140925/overhang74'; 
    }

}


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
my $CUTADAPT = '';

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my %samp_libs_run = ();
my $slr_count = 0;
my %samp_pair = ();

open(LOG, ">$cd\_rnaseq_pipeline.log") or die "can't write to output log";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSTARTING RNASEQ PIPELINE FOR $pre\n";
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCOMMAND LINE: $commandLine\n";

### Check that all programs are available
&verifyConfig($config);

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


if($comparisons || $samplekey){
    if(-e $comparisons && -e $samplekey){
	### GOING TO ASSUME THAT YOU WANT TO RUN DESEQ IF COMPARISONS & SAMPLEKEY FILE PROVIDED
	### IN ORDER TO RUN DESEQ, MUST ALSO RUN HTSEQ
	$deseq = 1;
	$htseq = 1;
    }
    else{
	die "COMPARISONS FILE $comparisons AND/OR SAMPLEKEY FILE $samplekey DON'T EXIST\nMUST PROVIDE BOTH COMPARISONS AND SAMPLEKEY FILES IN ORDER TO RUN DESEQ\n";
    }
}   

if($deseq || $dexseq || $htseq || $cufflinks){
    if(!$star && !$tophat){
	$star = 1;

	my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tNO ALIGNER SPECIFIED SO DEFAULTING TO USING STAR\n";
    }
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
if($strand =~ /none/i){
    $htseq_stranded = 'no';
    $picard_strand_specificity = 'NONE';
}
elsif($strand =~ /forward/i){
    $htseq_stranded = 'yes';
    $picard_strand_specificity = 'FIRST_READ_TRANSCRIPTION_STRAND';
}
elsif($strand =~ /reverse/i){
    $htseq_stranded = 'reverse';
    $picard_strand_specificity = 'SECOND_READ_TRANSCRIPTION_STRAND';
}


`/bin/mkdir -m 775 -p $output`; 
`/bin/mkdir -m 775 -p $output/intFiles`; 
`/bin/mkdir -m 775 -p $output/progress`;
`/bin/mkdir -m 775 -p $output/alignments`;
`/bin/mkdir -m 775 -p $output/metrics`;

my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

open(IN, "$map") or die "Can't open $map $!";
while(<IN>){
    chomp;

    my @data = split(/\s+/, $_);
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
	    }
	}
    }
    close OUT;
    chdir $curDir;
}
close IN;

my @crm = ();
my @crm_tophat = ();
my @asm = ();
my @asm_tophat = ();
my @cag_jids = ();
my $ran_cag = 0;
my $ran_tophatcrm = 0;
my @tophatcrm_jids = ();
my $ran_tophatasm = 0;
my @tophatasm_jids = ();
my $ran_starcrm = 0;
my @starcrm_jids = ();
my $ran_starasm = 0;
my @starasm_jids = ();
my $ran_tophathtseq = 0;
my @tophathtseq_jids = ();
my $ran_tophatdexseq = 0;
my @tophatdexseq_jids = ();
my $ran_starhtseq = 0;
my @starhtseq_jids = ();
my $ran_stardexseq = 0;
my @stardexseq_jids = ();

foreach my $sample (keys %samp_libs_run){
    my @R1 = ();
    my @R2 = ();
    my $readsFlag = 0;
    my $grl = 0;
    my $minReadLength = -1;
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
			`/bin/zcat $output/intFiles/$sample/$lib/$run/$readPair[0] >$output/intFiles/$readPair[0]`;
			### discard reads less than half of original read length by cutadapt
			my $readO = `/usr/bin/head -2 $output/intFiles/$readPair[0]`;
			chomp $readO;
			my @dataO = split(/\n/, $readO);
			my $readLength = length($dataO[1]);
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
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r1_gz_files_TRIM ">$output/intFiles/$sample/$sample\_R1.fastq"`;
	    `/bin/touch $output/progress/$pre\_$uID\_ZCAT_$sample\_R1.done`;
	    push @zcat_jids, "$pre\_$uID\_ZCAT_$sample\_R1";
	    $ran_zcat = 1;
	}

	if($samp_pair{$sample} eq "PE"){
	    if(!-e "$output/progress/$pre\_$uID\_ZCAT_$sample\_R2.done"){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$sample\_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT_$sample\_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r2_gz_files_TRIM ">$output/intFiles/$sample/$sample\_R2.fastq"`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR1 --paired-output $output/intFiles/$sample/$sample\_R2_TEMP.fastq -o $output/intFiles/$sample/$sample\_R1_TEMP.fastq $output/intFiles/$sample/$sample\_R1.fastq $output/intFiles/$sample/$sample\_R2.fastq >$output/intFiles/$sample/$sample\_R1\_CUTADAPT\_STATS.txt"`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR2 --paired-output $output/intFiles/$sample/$sample\_R1_CT.fastq -o $output/intFiles/$sample/$sample\_R2_CT.fastq $output/intFiles/$sample/$sample\_R2_TEMP.fastq $output/intFiles/$sample/$sample\_R1_TEMP.fastq >$output/intFiles/$sample/$sample\_R2\_CUTADAPT\_STATS.txt"`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/gzip $output/intFiles/$sample/$sample\_R1_CT.fastq`;
		`/bin/touch $output/progress/$pre\_$uID\_GZIP_$sample\_R1.done`;
		push @gz_jids, "$pre\_$uID\_GZIP_$sample\_R1";
		$ran_gz = 1;
	    }

	    if(!-e "$output/progress/$pre\_$uID\_GZIP_$sample\_R2.done" || $ran_ca){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GZIP_$sample\_R2", job_hold => "$caj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GZIP_$sample\_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/gzip $output/intFiles/$sample/$sample\_R2_CT.fastq`;
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
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUTADAPT_$sample\_R1 -hold_jid $zcatj -pe alloc 1 -l virtual_free=1G $Bin/qCMD "$PYTHON/python $CUTADAPT/cutadapt -f fastq -m $minReadLength $processR1 -o $output/intFiles/$sample/$sample\_R1_CT.fastq $output/intFiles/$sample/$sample\_R1.fastq >$output/intFiles/$sample/$sample\_R1\_CUTADAPT\_STATS.txt"`;
		`/bin/touch $output/progress/$pre\_$uID\_CUTADAPT_$sample\_R1.done`;
		push @ca_jids, "$pre\_$uID\_CUTADAPT_$sample\_R1";
		$ran_ca = 1;
	    }

    	    my $caj = join(",", @ca_jids);
	    if(!-e "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.done" || $ran_ca){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GZIP_$sample\_R1", job_hold => "$caj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_GZIP_$sample\_R1.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/gzip $output/intFiles/$sample/$sample\_R1_CT.fastq`;
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
	`/bin/mkdir -m 775 -p $output/alignments/tophat2`;
	`/bin/mkdir -m 775 -p $output/metrics/tophat2`;

	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}
	
	my $ran_tophat = 0;
	my $tophatj = '';
	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_$sample.done" || $ran_gz){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_$sample", job_hold => "$gzj", cpu => "6", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $TOPHAT/tophat2 -p 6 -r 70 --mate-std-dev 90 --GTF $GTF --transcriptome-index=$TRANS_INDEX -o $output/intFiles/tophat2/$sample --rg-id $sample\_1 --rg-sample $sample --rg-library $sample\_1 --rg-platform Illumina --rg-platform-unit $sample\_1 $BOWTIE2_INDEX $inReads`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_$sample.done`;
	    $tophatj = "$pre\_$uID\_TOPHAT_$sample";
	    $ran_tophat = 1;
	}

	my $ran_reorder = 0;
	my $reorderj = '';
	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_$sample.done" || $ran_tophat){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_REORDER_$sample", job_hold => "$tophatj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_REORDER_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar ReorderSam I=$output/intFiles/tophat2/$sample/accepted_hits.bam O=$output/alignments/tophat2/$pre\_$sample\.bam REFERENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_$sample.done`;
	    $reorderj = "$pre\_$uID\_REORDER_$sample";
	    $ran_reorder = 1;
	}
	
	if($cufflinks){
	    `/bin/mkdir -m 775 -p $output/cufflinks`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2/$sample`;
	    if(!-e "$output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.done" || $ran_tophat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUFFLINKS_TOPHAT_$sample", job_hold => "$tophatj", cpu => "5", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o $output/cufflinks/tophat2/$sample $output/intFiles/tophat2/$sample/accepted_hits.bam`;
		`/bin/touch $output/progress/$pre\_$uID\_CUFFLINKS_TOPHAT_$sample.done`;
	    }
	}
	
	if($htseq || $dexseq){
	    my $ran_tophatqns = 0;
	    my $tophatqnsj = '';
	    if(!-e "$output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.done" || $ran_tophat){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QNS_TOPHAT_$sample", job_hold => "$tophatj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles INPUT=$output/intFiles/tophat2/$sample/accepted_hits.bam OUTPUT=$output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam SORT_ORDER=queryname TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT`;
		`/bin/touch $output/progress/$pre\_$uID\_QNS_TOPHAT_$sample.done`;
		$tophatqnsj = "$pre\_$uID\_QNS_TOPHAT_$sample";
		$ran_tophatqns = 1;
	    }
	    
	    if($htseq){
		`/bin/mkdir -m 775 -p $output/counts_gene`;
		`/bin/mkdir -m 775 -p $output/counts_gene/tophat2`;
		
		if(!-e "$output/progress/$pre\_$uID\_HT_TOPHAT_$sample.done" || $ran_tophatqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HT_TOPHAT_$sample", job_hold => "$tophatqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_HT_TOPHAT_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $HTSEQ/htseq-count -m intersection-strict -s $htseq_stranded -t exon $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $GTF > $output/counts_gene/tophat2/$sample.htseq_count"`;
		    `/bin/touch $output/progress/$pre\_$uID\_HT_TOPHAT_$sample.done`;
		    push @tophathtseq_jids, "$pre\_$uID\_HT_TOPHAT_$sample";
		    $ran_tophathtseq = 1;
		}
	    }
	    
	    if($dexseq){
		`/bin/mkdir -m 775 -p $output/counts_exon`;
		`/bin/mkdir -m 775 -p $output/counts_exon/tophat2`;
		
     		if(!-e "$output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.done" || $ran_tophatqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEXSEQ_TOPHAT_$sample", job_hold => "$tophatqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $output/counts_exon/tophat2/$sample.dexseq_count`;
		    `/bin/touch $output/progress/$pre\_$uID\_DEXSEQ_TOPHAT_$sample.done`;
		    push @tophatdexseq_jids, "$pre\_$uID\_DEXSEQ_TOPHAT_$sample";
		    $ran_tophatdexseq = 1;
		}
	    }
	}

	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.done" || $ran_reorder){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_CRM_$sample", job_hold => "$reorderj", cpu => "1", mem => "3", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/alignments/tophat2/$pre\_$sample\.bam O=$output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/tophat2/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=NONE METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_CRM_$sample.done`;
	    push @tophatcrm_jids, "$pre\_$uID\_TOPHAT_CRM_$sample";
	    $ran_tophatcrm = 1;
	}	
	push @crm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	if(!-e "$output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.done" || $ran_reorder){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TOPHAT_ASM_$sample", job_hold => "$reorderj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/alignments/tophat2/$pre\_$sample\.bam OUTPUT=$output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_TOPHAT_ASM_$sample.done`;
	    push @tophatasm_jids, "$pre\_$uID\_TOPHAT_ASM_$sample";
	    $ran_tophatasm = 1;
	}
	push @asm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt";
    }

    my $starOut = '';    
    if($star){
	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}

	my $ran_star1p = 0;
	my $star1pj = '';
	if(!-e "$output/progress/$pre\_$uID\_STAR_1PASS_$sample.done" || $ran_gz){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_1PASS_$sample", job_hold => "$gzj", cpu => "12", mem => "40", cluster_out => "$output/progress/$pre\_$uID\_STAR_1PASS_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_1PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_1PASS_Aligned.out.sam`;
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
	    if(!-e "$curDir/$output/progress/$pre\_$uID\_STAR_GG2_$sample.done" || $ran_star1p){
		sleep(3);
		chdir "$curDir/$output/intFiles/$sample/star2passGG";
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_GG2_$sample", job_hold => "$star1pj", cpu => "12", mem => "40", cluster_out => "$curDir/$output/progress/$pre\_$uID\_STAR_GG2_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STAR/STAR --runMode genomeGenerate --genomeDir $curDir/$output/intFiles/$sample/star2passGG --genomeFastaFiles $REF_SEQ --sjdbFileChrStartEnd $curDir/$output/intFiles/$sample/$sample\_STAR_1PASS_SJ.out.tab --sjdbOverhang $mrl --runThreadN 12`;
		`/bin/touch $curDir/$output/progress/$pre\_$uID\_STAR_GG2_$sample.done`;
		$sgg2j = "$pre\_$uID\_STAR_GG2_$sample";
		$ran_sgg2 = 1;
		chdir $curDir;
	    }
	    
	    my $ran_star2p = 0;
	    my $star2pj = '';
	    if(!-e "$output/progress/$pre\_$uID\_STAR_2PASS_$sample.done" || $ran_sgg2){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_2PASS_$sample", job_hold => "$sgg2j", cpu => "12", mem => "40", cluster_out => "$output/progress/$pre\_$uID\_STAR_2PASS_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STAR/STAR --genomeDir $output/intFiles/$sample/star2passGG --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_2PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
		`/bin/touch $output/progress/$pre\_$uID\_STAR_2PASS_$sample.done`;
		$star2pj = "$pre\_$uID\_STAR_2PASS_$sample";
		$ran_star2p = 1;
	    }	    
	    
	    if(!-e "$output/progress/$pre\_$uID\_SP_$sample.done" || $ran_star2p){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SP_$sample", job_hold => "$star2pj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_SP_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam`;
		`/bin/touch $output/progress/$pre\_$uID\_SP_$sample.done`;
		$spj = "$pre\_$uID\_SP_$sample";
		$ran_sp = 1;
	    }
	    $starOut = "$output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.sam";
	}

	my $ran_starmerge = 0;
	my $starmergej = '';
	if(!-e "$output/progress/$pre\_$uID\_STAR_MERGE_$sample.done" || $ran_sp){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_MERGE_$sample", job_hold => "$spj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_STAR_MERGE_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles I=$starOut O=$starOut\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_MERGE_$sample.done`;
	    $starmergej = "$pre\_$uID\_STAR_MERGE_$sample";
	    $ran_starmerge = 1;
	}

	my $ran_staraddrg = 0;
	my $staraddrgj = '';
	if(!-e "$output/progress/$pre\_$uID\_STAR_AORRG_$sample.done" || $ran_starmerge){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_AORRG_$sample", job_hold => "$starmergej", cpu => "3", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_STAR_AORRG_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$starOut\.bam O=$output/alignments/$pre\_$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true RGID=$sample\_1 RGLB=_1 RGPL=Illumina RGPU=$sample\_1 RGSM=$sample`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_AORRG_$sample.done`;
	    $staraddrgj = "$pre\_$uID\_STAR_AORRG_$sample";
	    $ran_staraddrg = 1;	
	}

	if($cufflinks){
	    `/bin/mkdir -m 775 -p $output/cufflinks`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/$sample`;
	    if(!-e "$output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.done" || $ran_starmerge){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUFFLINKS_STAR_$sample", job_hold => "$starmergej", cpu => "5", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o $output/cufflinks/$sample $starOut\.bam`;
		`/bin/touch $output/progress/$pre\_$uID\_CUFFLINKS_STAR_$sample.done`;
	    }
	}

	if($htseq || $dexseq){
	    my $ran_starqns = 0;
	    my $starqnsj = '';
	    if(!-e "$output/progress/$pre\_$uID\_QNS_STAR_$sample.done" || $ran_sp){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QNS_STAR_$sample", job_hold => "$spj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_QNS_STAR_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles I=$starOut O=$starOut\_queryname_sorted.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID USE_THREADING=true`;
		`/bin/touch $output/progress/$pre\_$uID\_QNS_STAR_$sample.done`;
		$starqnsj = "$pre\_$uID\_QNS_STAR_$sample";
		$ran_starqns = 1;
	    }	    
	
	    if($htseq){
		`/bin/mkdir -m 775 -p $output/counts_gene`;	   
		if(!-e "$output/progress/$pre\_$uID\_HT_STAR_$sample.done" || $ran_starqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HT_STAR_$sample", job_hold => "$starqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_HT_STAR_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $HTSEQ/htseq-count -m intersection-strict -s $htseq_stranded -t exon $starOut\_queryname_sorted.sam $GTF > $output/counts_gene/$sample.htseq_count"`;
		    `/bin/touch $output/progress/$pre\_$uID\_HT_STAR_$sample.done`;
		    push @starhtseq_jids, "$pre\_$uID\_HT_STAR_$sample";
		    $ran_starhtseq = 1;
		}
	    }
	
	    if($dexseq){
		`/bin/mkdir -m 775 -p $output/counts_exon`;
     		if(!-e "$output/progress/$pre\_$uID\_DEX_STAR_$sample.done" || $ran_starqns){
		    sleep(3);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEX_STAR_$sample", job_hold => "$starqnsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DEX_STAR_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $starOut\_queryname_sorted.sam $output/counts_exon/$sample.dexseq_count`;
		    `/bin/touch $output/progress/$pre\_$uID\_DEX_STAR_$sample.done`;
		    push @stardexseq_jids, "$pre\_$uID\_DEX_STAR_$sample";
		    $ran_stardexseq = 1;
		}
	    }
	}

	my $ran_starcrm = 0;
	my @starcrm = ();
	if(!-e "$output/progress/$pre\_$uID\_STAR_CRM_$sample.done" || $ran_staraddrg){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_CRM_$sample", job_hold => "$staraddrgj", cpu => "1", mem => "3", cluster_out => "$output/progress/$pre\_$uID\_STAR_CRM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/alignments/$pre\_$sample\.bam O=$output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=$picard_strand_specificity METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_CRM_$sample.done`;
	    push @starcrm_jids, "$pre\_$uID\_STAR_CRM_$sample";
	    $ran_starcrm = 1;
	}
	push @crm, "-metrics $output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	my $ran_starasm = 0;
	my @starasm = ();
	if(!-e "$output/progress/$pre\_$uID\_STAR_ASM_$sample.done" || $ran_staraddrg){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_ASM_$sample", job_hold => "$staraddrgj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_STAR_ASM_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/alignments/$pre\_$sample\.bam OUTPUT=$output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_STAR_ASM_$sample.done`;
	    push @starasm_jids, "$pre\_$uID\_STAR_ASM_$sample";
	    $ran_starasm = 1;
	}
	push @asm, "-metrics $output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt";
    }

    if($detectFusions){
	`/bin/mkdir -m 775 -p $output/fusion`;
	if($species =~ /hg19|human/i){	
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r1_files ">$output/intFiles/$sample/$sample\_v2_R1.fastq"`;
		`/bin/touch $output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R1.done`;
		push @zcat2_jids, "$pre\_$uID\_ZCAT2_$sample\_v2_R1";
		$ran_zcat2 = 1;
		sleep(3);
	    }

	    if($samp_pair{$sample} eq "PE"){
		if(!-e "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.done"){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT2_$sample\_v2_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r2_files ">$output/intFiles/$sample/$sample\_v2_R2.fastq"`;
		    `/bin/touch $output/progress/$pre\_$uID\_ZCAT2_$sample\_v2_R2.done`;
		    push @zcat2_jids, "$pre\_$uID\_ZCAT2_$sample\_v2_R2";
		    $ran_zcat2 = 1;
		    sleep(3);
		}
	    }	

	    my $zcat2j = join(",", @zcat2_jids);
	    if($chimerascan){
		### NOTE: CHIMERASCAN ONLY WORKS FOR PE READS
		if($samp_pair{$sample} eq "PE"){
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
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $CHIMERASCAN/chimerascan_run.py -p 6 --quals solexa --multihits=10 --filter-false-pos=$Bin/data/hg19_bodymap_false_positive_chimeras.txt $CHIMERASCAN_INDEX $output/intFiles/$sample/$sample\_v2_R1.fastq $output/intFiles/$sample/$sample\_v2_R2.fastq $output/fusion/chimerascan/$sample/`;
			`/bin/touch $output/progress/$pre\_$uID\_CHIMERASCAN_$sample.done`;
			push @fusion_jids, "$pre\_$uID\_CHIMERASCAN_$sample";
			$ran_fusion = 1;
		    }
		    push @fusions, "--chimerascan $output/fusion/chimerascan/$sample/chimeras.bedpe";
		}
	    }

	    if($star_fusion){
		`/bin/mkdir -m 775 -p $output/fusion/star`;
		`/bin/mkdir -m 775 -p $output/fusion/star/$sample`;
		
		my $inReads = "$output/intFiles/$sample/$sample\_v2_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " $output/intFiles/$sample/$sample\_v2_R2.fastq";
		}

		if(!-e "$output/progress/$pre\_$uID\_STAR_CHIMERA_$sample.done" || $ran_zcat2){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STAR_CHIMERA_$sample", job_hold => "$zcat2j", cpu => "6", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_STAR_CHIMERA_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 6 --outFileNamePrefix $output/fusion/star/$sample/$sample\_STAR_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --chimSegmentMin 20`;
		    `/bin/touch $output/progress/$pre\_$uID\_STAR_CHIMERA_$sample.done`;
		    push @fusion_jids, "$pre\_$uID\_STAR_CHIMERA_$sample";
		    $ran_fusion = 1;
		}
		push @fusions, "--star $output/fusion/star/$sample/$sample\_STAR_Chimeric.out.junction";
	    }

	    if($mapsplice){
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice`;
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice/$sample`;

		my $inReads = "-1 $output/intFiles/$sample/$sample\_v2_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " -2 $output/intFiles/$sample/$sample\_v2_R2.fastq"
		}

		if(!-e "$output/progress/$pre\_$uID\_MAPSPLICE_$sample.done" || $ran_zcat2){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MAPSPLICE_$sample", job_hold => "$zcat2j", cpu => "6", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_MAPSPLICE_$sample.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $MAPSPLICE/mapsplice.py -p 6 --bam --fusion-non-canonical -c $chrSplits -x $BOWTIE_INDEX -o $output/fusion/mapsplice/$sample $inReads --gene-gtf $GTF`;
		    `/bin/touch $output/progress/$pre\_$uID\_MAPSPLICE_$sample.done`;
		    push @fusion_jids, "$pre\_$uID\_MAPSPLICE_$sample";
		    $ran_fusion = 1;
		}
		push @fusions, "--mapsplice $output/fusion/mapsplice/$sample/fusions_well_annotated.txt";
	    }

	    if($defuse){
		### NOTE: DEFUSE ONLY WORKS FOR PE READS
		my $defuse_config = '';
		if($species =~ /human|hg19/i){
		    $defuse_config = "$DEFUSE/scripts/config_homo_sapiens.txt";
		}
		elsif($species =~ /mouse|mm10/i){
		    $defuse_config = "$DEFUSE/scripts/config_mus_musculus.txt";
		}
		
		if($samp_pair{$sample} eq "PE"){
		    `/bin/mkdir -m 775 -p $output/fusion/defuse`;
		    `/bin/mkdir -m 775 -p $output/fusion/defuse/$sample`;
				
		    if(!-e "$output/progress/$pre\_$uID\_DEFUSE_$sample.done" || $ran_zcat2){
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DEFUSE_$sample", job_hold => "$zcat2j", cpu => "12", mem => "150", cluster_out => "$output/progress/$pre\_$uID\_DEFUSE_$sample.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $DEFUSE/scripts/defuse.pl --config $defuse_config --output $output/fusion/defuse/$sample --parallel 12 --1fastq $output/intFiles/$sample/$sample\_v2_R1.fastq --2fastq $output/intFiles/$sample/$sample\_v2_R2.fastq`;
			`/bin/touch $output/progress/$pre\_$uID\_DEFUSE_$sample.done`;
			push @fusion_jids, "$pre\_$uID\_DEFUSE_$sample";
			$ran_fusion = 1;
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
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $FUSIONCATCHER/bin/fusioncatcher -d $FUSIONCATCHER/data/ensembl_v77d $inReads -o $output/fusion/fusioncatcher/$sample -p 6 --skip-update-check --config=$FUSIONCATCHER/bin/configuration.cfg`;
		    `/bin/touch $output/progress/$pre\_$uID\_FC_$sample.done`;
		    push @fusion_jids, "$pre\_$uID\_FC_$sample";
		    $ran_fusion = 1;
		}
		push @fusions, "--fusioncatcher $output/fusion/fusioncatcher/$sample/final-list_candidate-fusion-genes.txt";
	    }

	    my $mergeFusions = join(" ", @fusions);
	    my $fusionj = join(",", @fusion_jids);
	    if(!-e "$output/progress/$pre\_$uID\_MERGE_FUSION_$sample.done" || $ran_fusion){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_FUSION_$sample", job_hold => "$fusionj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_FUSION_$sample.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/MergeFusion $mergeFusions --out $output/fusion/$pre\_merged_fusions_$sample\.txt --normalize_gene $Bin/data/hugo_data_073013.tsv`;
		`/bin/touch $output/progress/$pre\_$uID\_MERGE_FUSION_$sample.done`;
	    }
	}
	else{
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSKIPPING FUSION CALLING FOR SAMPLE $sample BECAUSE IT ISN'T SUPPORTED FOR $species; CURRENTLY ONLY SUPPORT FOR hg19";
	}
    }

    if($transcript){
	my $r1_gz_files1 = join(" ", @R1);
	my $r2_gz_files1 = join(" ", @R2);
	my $r1_gz_files2 = join(",", @R1);
	my $r2_gz_files2 = join(",", @R2);
	
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
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_BOWTIE2_$sample", cpu => "12", mem => "60", cluster_out => "$output/progress/$pre\_$uID\_BOWTIE2_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $BOWTIE2/bowtie2 --all --maxins 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed --threads 12 -x $TRANS_INDEX_DEDUP $inReads -S $output/intFiles/bowtie2/$sample/$sample\_bowtie2.sam`;
	    `/bin/touch $output/progress/$pre\_$uID\_BOWTIE2_$sample.done`;
	    $bowtie2j = "$pre\_$uID\_BOWTIE2_$sample";
	    $ran_bowtie2 = 1;
	}
	
	if(!-e "$output/progress/$pre\_$uID\_EXPRESS_$sample.done" || $ran_bowtie2){
	    sleep(3);
	    `/bin/mkdir -m 775 -p $output/transcript`;
	    `/bin/mkdir -m 775 -p $output/transcript/express`;
	    `/bin/mkdir -m 775 -p $output/transcript/express/$sample`;
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_EXPRESS_$sample", job_hold => "$bowtie2j", cpu => "5", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_EXPRESS_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $EXPRESS/express --output-dir $output/transcript/express/$sample --no-update-check $TRANS_FASTA_DEDUP $output/intFiles/bowtie2/$sample/$sample\_bowtie2.sam`;
	    `/bin/touch $output/progress/$pre\_$uID\_EXPRESS_$sample.done`;
	}
	
	
	my $ran_zcat3 = 0;
	my @zcat3_jids = ();
	my $kinReads = "$output/intFiles/$sample/$sample\_v3_R1.fastq";
	if(!-e "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.done"){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT3_$sample\_v3_R1", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r1_gz_files1 ">$output/intFiles/$sample/$sample\_v3_R1.fastq"`;
	    `/bin/touch $output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R1.done`;
	    push @zcat3_jids, "$pre\_$uID\_ZCAT3_$sample\_v3_R1";
	    $ran_zcat3 = 1;
	}
	
	if($samp_pair{$sample} eq "PE"){
	    $kinReads .= " $output/intFiles/$sample/$sample\_v3_R2.fastq";
	    if(!-e "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.done"){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT3_$sample\_v3_R2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/zcat $r2_gz_files1 ">$output/intFiles/$sample/$sample\_v3_R2.fastq"`;
		`/bin/touch $output/progress/$pre\_$uID\_ZCAT3_$sample\_v3_R2.done`;
		push @zcat3_jids, "$pre\_$uID\_ZCAT3_$sample\_v3_R2";
		$ran_zcat3 = 1;
	    }
	}
	    
	my $zcat3j = join(",", @zcat3_jids);
	if(!-e "$output/progress/$pre\_$uID\_KALLISTO_$sample.done" || $ran_zcat3){
	    sleep(3);
	    `/bin/mkdir -m 775 -p $output/transcript`;
	    `/bin/mkdir -m 775 -p $output/transcript/kallisto`;
	    `/bin/mkdir -m 775 -p $output/transcript/kallisto/$sample`;
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_KALLISTO_$sample", job_hold => "$zcat3j", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_KALLISTO_$sample.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $KALLISTO/kallisto quant -i $KALLISTO_INDEX -o $output/transcript/kallisto/$sample -b 100 $kinReads`;
	    `/bin/touch $output/progress/$pre\_$uID\_KALLISTO_$sample.done`;
	}
    }
}

my $ran_shmatrix = 0;
my $shmatrixj = '';
if($star){
    if($htseq){
	my $starhtseqj = join(",", @starhtseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.done" || $ran_starhtseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_HTSEQ_STAR", job_hold => "$starhtseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_gene .htseq_count $curDir/$output/counts_gene/$pre\_htseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_HTSEQ_STAR.done`;
	    $shmatrixj = "$pre\_$uID\_MATRIX_HTSEQ_STAR";
	    $ran_shmatrix = 1;
	}
    }
    if($dexseq){
	my $stardexseqj = join(",", @stardexseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_DEX_STAR.done" || $ran_stardexseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_DEX_STAR", job_hold => "$stardexseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_DEX_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_exon .dexseq_count $curDir/$output/counts_exon/$pre\_dexseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_DEX_STAR.done`;
	}
    }

    my $crmfiles = join(" ", @crm);
    my $scrmj = join(",", @starcrm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_STAR.done" || $ran_starcrm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_STAR", job_hold => "$scrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_STAR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/mergePicardMetrics.pl $crmfiles ">$output/metrics/$pre\_CollectRnaSeqMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_STAR.done`;
    }

    my $asmfiles = join(" ", @asm);
    my $sasmj = join(",", @starasm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM_STAR.done" || $ran_starasm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ASM_STAR", job_hold => "$sasmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ASM_STAR.log");
	my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/mergePicardMetrics.pl $asmfiles ">$output/metrics/$pre\_AlignmentSummaryMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM_STAR.done`;
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
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_gene/tophat2 .htseq_count $curDir/$output/counts_gene/tophat2/$pre\_htseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_HTSEQ_TOPHAT.done`;
	    $thmatrixj = "$pre\_$uID\_MATRIX_HTSEQ_TOPHAT";
	    $ran_thmatrix = 1;
	}
    }

    if($dexseq){
	my $tophatdexseqj = join(",", @tophatdexseq_jids);
	if(!-e "$output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.done" || $ran_tophatdexseq){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MATRIX_DEX_TOPHAT", job_hold => "$tophatdexseqj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_exon/tophat2 .dexseq_count $curDir/$output/counts_exon/tophat2/$pre\_dexseq_all_samples.txt $geneNameConversion`;
	    `/bin/touch $output/progress/$pre\_$uID\_MATRIX_DEX_TOPHAT.done`;
	}
    }

    my $crmfiles_tophat = join(" ", @crm_tophat);
    my $tcrmj = join(",", @tophatcrm_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.done" || $ran_tophatcrm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CRM_TOPHAT", job_hold => "$tcrmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/mergePicardMetrics.pl $crmfiles_tophat ">$output/metrics/tophat2/$pre\_CollectRnaSeqMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CRM_TOPHAT.done`;
    }

    my $asmfiles_tophat = join(" ", @asm_tophat);
    my $tasmj = join(",", @tophatasm_jids);

    if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.done" || $ran_tophatasm){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ASM_TOPHAT", job_hold => "$tasmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/mergePicardMetrics.pl $asmfiles_tophat ">$output/metrics/tophat2$pre\_AlignmentSummaryMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM_TOPHAT.done`;
    }
}

### SAMPLE COMMAND
####qsub -N Proj_4226_DESeq ~/bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript ~/RNAseqPipe/trunk/bin/RunDE.R "\"proj.id='4226'\" \"output.dir='/ifs/data/byrne/rnaseq/Proj_4226'\" \"counts.file='Proj_4226_ALL_samples.htseq.count'\" \"key.file='/ifs/data/byrne/rnaseq/Proj_4226/sampleKey.txt'\" \"comps=c('hi - lo')\""


if($deseq){    
    if($star){
	`/bin/mkdir -m 775 -p $output/differentialExpression_gene`;
	`/bin/mkdir -m 775 -p $output/clustering`;
	`/bin/mkdir -m 775 -p $output/gsa`;
	###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq_STAR -hold_jid $pre\_$uID\_MATRIX_HTSEQ_STAR -pe alloc 1 -l virtual_free=1G $Bin/qCMD $R/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"species='$species'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"diff.exp.dir='$curDir/$output/differentialExpression_gene'\\\"\" \"\\\"counts.file='$curDir/$output/counts_gene/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"counts.dir='$curDir/$output/counts_gene'\\\"\" \"\\\"clustering.dir='$curDir/$output/clustering'\\\"\" \"\\\"gsa.dir='$curDir/$output/gsa'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
	if(!-e "$output/progress/$pre\_$uID\_DESeq_STAR.done" || $ran_shmatrix){
            my $reps = '';
            if($no_replicates){
                $reps = '-no_replicates';
            }
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_STAR", job_hold => "$shmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DESeq_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/run_DESeq_wrapper.pl -pre $pre -diff_out $curDir/$output/differentialExpression_gene -count_out $curDir/$output/counts_gene -cluster_out $curDir/$output/clustering -gsa_out $curDir/$output/gsa -config $config -bin $Bin -species $species -counts $curDir/$output/counts_gene/$pre\_htseq_all_samples.txt -samplekey $samplekey -comparisons $comparisons $reps`;
	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_STAR.done`;
	}
    }

    if($tophat){
	`/bin/mkdir -m 775 -p $output/differentialExpression_gene/tophat2`;
	`/bin/mkdir -m 775 -p $output/clustering/tophat2`;
	`/bin/mkdir -m 775 -p $output/gsa/tophat2`;

	if(!-e "$output/progress/$pre\_$uID\_DESeq_TOPHAT.done" || $ran_thmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DESeq_TOPHAT", job_hold => "$thmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_DESeq_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/run_DESeq_wrapper.pl -pre $pre -diff_out $curDir/$output/differentialExpression_gene/tophat2 -count_out $curDir/$output/counts_gene/tophat2 -cluster_out $curDir/$output/clustering/tophat2 -gsa_out $curDir/$output/gsa/tophat2 -config $config -bin $Bin -species $species -counts $curDir/$output/counts_gene/tophat2/$pre\_htseq_all_samples.txt -samplekey $samplekey -comparisons $comparisons`;
	    `/bin/touch $output/progress/$pre\_$uID\_DESeq_TOPHAT.done`;
	}
    }
}
else{
    if($star){
	`/bin/mkdir -m 775 -p $output/clustering`;
	if(!-e "$output/progress/$pre\_$uID\_CLUSTERING_STAR.done" || $ran_shmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTERING_STAR", job_hold => "$shmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_CLUSTERING_STAR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/run_DESeq_wrapper.pl -cluster_out $curDir/$output/clustering -config $config -bin $Bin -counts $curDir/$output/counts_gene/$pre\_htseq_all_samples.txt -clusterOnly`;
	    `/bin/touch $output/progress/$pre\_$uID\_CLUSTERING_STAR.done`;
	}
    }

    if($tophat){
	`/bin/mkdir -m 775 -p $output/clustering/tophat2`;
	if(!-e "$output/progress/pre\_$uID\_CLUSTERING_TOPHAT.done" || $ran_thmatrix){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "pre\_$uID\_CLUSTERING_TOPHAT", job_hold => "$thmatrixj", cpu => "1", mem => "1", cluster_out => "$output/progress/pre\_$uID\_CLUSTERING_TOPHAT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/run_DESeq_wrapper.pl -cluster_out $curDir/$output/clustering/tophat2 -config $config -bin $Bin -counts $curDir/$output/counts_gene/tophat2/$pre\_htseq_all_samples.txt -clusterOnly`;	    
	}
    }
}


if($r1adaptor){
    my $cagj = join(",", @cag_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_CAS.done" || $ran_cag){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", job_hold => "$cagj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/mergeCutAdaptStats.py . '*CUTADAPT_STATS.txt' $output/metrics/$pre\_CutAdaptStats.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CAS.done`;
    }
}

close LOG;


sub verifyConfig{
    my $paths = shift;

    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
	chomp;
	
	my @conf = split(/\s+/, $_);

	if($conf[0] =~ /picard/i){
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
	}
	elsif($conf[0] =~ /star/i){
	    if(!-e "$conf[1]/STAR"){
		die "CAN'T FIND STAR IN $conf[1] $!";
	    }
	    $STAR = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
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
	    $BOWTIE2 = $conf[1];
	}
	elsif($conf[0] =~ /^bowtie$/i){
	    ### need bowtie in path for chimerascan to run
	    if(!-e "$conf[1]/bowtie"){
		die "CAN'T FIND bowtie IN $conf[1] $!";
	    }
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";	    
	}
	elsif($conf[0] =~ /java/i){
	    if(!-e "$conf[1]/java"){
		die "CAN'T FIND java IN $conf[1] $!";
	    }
	    $JAVA = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
	}
	elsif($conf[0] =~ /perl/i){
	    if(!-e "$conf[1]/perl"){
		die "CAN'T FIND perl IN $conf[1] $!";
	    }
	    $PERL = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
	}
 	elsif($conf[0] =~ /python/i){
	    if(!-e "$conf[1]/python"){
		die "CAN'T FIND python IN $conf[1] $!";
	    }
	    $PYTHON = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
	}
	elsif($conf[0] =~ /^r$/i){
	    if(!-e "$conf[1]/R"){
		die "CAN'T FIND R IN $conf[1] $!";
	    }
	    $R = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
	}
	elsif($conf[0] =~ /cutadapt/i){
	    if(!-e "$conf[1]/cutadapt"){
		die "CAN'T FIND cutadapt IN $conf[1] $!";
	    }
	    $CUTADAPT = $conf[1];
	}
    }
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

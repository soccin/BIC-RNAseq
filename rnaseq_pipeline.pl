#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

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


my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1, $lncrna, $lincrna_BROAD, $output);

$pre = 'TEMP';
$output = "results";
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
            'species=s' => \$species,
            'lncrna' => \$lncrna,
 	    'output|out|o=s' => \$output,
            'lincrna_BROAD' => \$lincrna_BROAD) or exit(1);


if(!$map || !$species || !$config || $help){
    print <<HELP;

    USAGE: ./rnaseq_pipeline.pl -map MAP -species SPECIES -config CONFIG -pre PRE -samplekey SAMPLEKEY -comparisons COMPARISONS
	* MAP: file listing sample mapping information for processing (REQUIRED)
	* SPECIES: only hg19 and mm9 and human-mouse hybrid (hybrid) currently supported (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B (if -deseq, REQUIRED)
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B (if -deseq, REQUIRED)
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless -pass1 specified; tophat2 (-tophat); if no aligner specifed, will default to STAR
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons); fusion callers chimerascan (-chimerascan), rna star (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs
	* OUTPUT: output results directory (default: results)
        * OPTIONS: lncRNA analysis (-lncrna) runs all analyses based on lncRNA GTF (hg19 only); 
HELP
exit;
}

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($species !~ /human|hg19|mouse|mm9|hybrid|zebrafish|zv9/i){
    die "Species must be human (hg19) or mouse (mm9) or human-mouse hybrid (hybrid) or zebrafish (zv9)";
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
if($samplekey){
    $commandLine .= " -samplekey $samplekey";
}
if($comparisons){
    $commandLine .= " -comparisons $comparisons";
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
if($species){
    $commandLine .= " -species $species";
}
if($lncrna){
    $commandLine .= " -lncrna";
}

my $numArgs = $#ARGV + 1;
foreach my $argnum (0 .. $#ARGV) {
    $commandLine .= " $ARGV[$argnum]";
}

my $GTF = '';
my $DEXSEQ_GTF = '';
my $CHIMERASCAN_INDEX = '';
my $geneNameConversion = '';
my $starDB = '';
my $chrSplits = '';
my $BOWTIE_INDEX = '';
my $BOWTIE2_INDEX = '';
my $TRANS_INDEX = '';
my $REF_SEQ = '';
my $RIBOSOMAL_INTERVALS='';
my $REF_FLAT = '';

if($species =~ /human|hg19/i){
    $species = 'hg19';
    $REF_SEQ = '/ifs/data/bio/assemblies/H.sapiens/hg19/hg19.fasta';
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $CHIMERASCAN_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/chimerascan/';
    $BOWTIE_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie/hg19_bowtie';
    $BOWTIE2_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie2/hg19_bowtie2';
    $chrSplits = '/ifs/data/bio/assemblies/H.sapiens/hg19/chromosomes';
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_hg19.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__hg19.txt.gz";

    if($lncrna){
        $GTF = "$Bin/data/lncipedia.gtf";
        $starDB = '/ifs/data/bio/assemblies/H.sapiens/hg19/star/LNCipedia';
        $geneNameConversion = '';
        $TRANS_INDEX = '';
    } else {
        $GTF = "$Bin/data/gencode.v18.annotation.gtf";
        $starDB = '/ifs/data/bio/assemblies/H.sapiens/hg19/star';
        $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
        $TRANS_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie2/transcriptome/gencode18';
    }
}
elsif($species =~ /mouse|mm9/i){
    $species = 'mm9';
    $REF_SEQ = '/ifs/data/bio/assemblies/M.musculus/mm9/mm9.fasta';
    $GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.gtf";
    $DEXSEQ_GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf";
    $starDB = '/ifs/data/bio/assemblies/M.musculus/mm9/star';
    $geneNameConversion = "$Bin/data/mm9Ensembl67IDToGeneName.txt";
    $BOWTIE_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie/mm9_bowtie';
    $BOWTIE2_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie2/mm9_bowtie2';
    $chrSplits = '/ifs/data/bio/assemblies/M.musculus/mm9/chromosomes';
    $TRANS_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie2/transcriptome/ensembl';
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_MM9_assemblies.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__mm9.txt.gz";
}
elsif($species =~ /human-mouse|mouse-human|hybrid/i){
    $species = 'hybrid';
    $REF_SEQ = '/ifs/data/mpirun/genomes/human/hg19_mm9_hybrid/hg19_mm9_hybrid.fasta';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $starDB = '/ifs/data/mpirun/genomes/human/hg19_mm9_hybrid/star';
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
    $RIBOSOMAL_INTERVALS = "$Bin/data/ribosomal_hg19.interval_file";
    $REF_FLAT = "$Bin/data/refFlat__hg19.txt.gz";
}
elsif($species =~ /zebrafish|zv9/i){
    $species = 'zv9';
    $REF_SEQ = '/ifs/data/bio/assemblies/D.rerio/zv9/zv9.fasta';
    $GTF = "$Bin/data/zv9.gtf";
    $starDB = '/ifs/data/bio/assemblies/D.rerio/zv9/star';
    $chrSplits = '/ifs/data/bio/assemblies/D.rerio/zv9/chromosomes';
    $geneNameConversion = "Bin/data/zv9EnsemblIDtoGeneName.txt";
}


my $CUFFLINKS = '/opt/bin';
my $HTSEQ = '/opt/bin';
my $DEXSEQ = '/opt/bin';
my $PICARD = '/opt/bin';
my $CHIMERASCAN = '/opt/bin';
my $MAPSPLICE = '/opt/bin';
my $STAR = '/opt/bin';
my $DEFUSE = '/opt/bin';
my $FUSIONCATCHER = '/opt/bin';
my $TOPHAT = '/opt/bin';

### THIS VALUE IS USED IN CREATING JUNCTION IS 2ND PASS OF STAR
my $MIN_READ_LENGTH = 50;

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

`/bin/mkdir -m 775 -p $output`; 
`/bin/mkdir -m 775 -p $output/intFiles`; 
`/bin/mkdir -m 775 -p $output/progress`;
`/bin/mkdir -m 775 -p $output/alignments`;
`/bin/mkdir -m 775 -p $output/metrics`;

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

    ###if(-e "$data[1]/$data[0]/$data[2]/FINISHED_$data[1]\_$data[0]\_$data[2]"){
	###print LOG "SKIPPING $_ because $data[1]/$data[0]/$data[2]/FINISHED_$data[1]\_$data[0]\_$data[2] exists and therefore considered processed";
	###next;
    ###}

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
foreach my $sample (keys %samp_libs_run){
    my @R1 = ();
    my @R2 = ();
    my $readsFlag = 0;
    foreach my $lib (keys %{$samp_libs_run{$sample}}){
	foreach my $run (keys %{$samp_libs_run{$sample}{$lib}}){	    
	    open(READS, "$output/intFiles/$sample/$lib/$run/files_$sample\_$lib\_$run") || die "Can't open $output/intFiles/$sample/$lib/$run/files_$sample\_$lib\_$run $!";
	    while(<READS>){
		chomp;
		    
		my @readPair = split(/\s+/, $_);
		
		### NOTE: JUST CHECKING TO SEE IF IT EXISTS
		###       HOWEVER DOES NOT GURANTEE THAT IT'S NON-EMPTY
		if(-e "$output/intFiles/$sample/$lib/$run/$readPair[0]" && -e "$output/intFiles/$sample/$lib/$run/$readPair[1]"){
		    push @R1, "$output/intFiles/$sample/$lib/$run/$readPair[0]";
		    push @R2, "$output/intFiles/$sample/$lib/$run/$readPair[1]";
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
    if($tophat){
	`/bin/mkdir -m 775 -p $output/intFiles/tophat2`;
	`/bin/mkdir -m 775 -p $output/intFiles/tophat2/$sample`;
	`/bin/mkdir -m 775 -p $output/alignments/tophat2`;
	`/bin/mkdir -m 775 -p $output/metrics/tophat2`;

	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}
	
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_TOPHAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD $TOPHAT/tophat2 -p 6 --zpacker /opt/pigz-2.1.6/pigz  -r 70 --mate-std-dev 90 --GTF $GTF --transcriptome-index=$TRANS_INDEX -o $output/intFiles/tophat2/$sample $BOWTIE2_INDEX $inReads`;

	sleep(5);

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_REORDER_$sample -hold_jid $pre\_$uID\_TOPHAT_$sample -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar ReorderSam I=$output/intFiles/tophat2/$sample/accepted_hits.bam O=$output/intFiles/tophat2/$sample/accepted_hits_RO.bam REFERENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true`;

	sleep(5);

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_TOPHAT_AORRG_$sample -hold_jid $pre\_$uID\_REORDER_$sample -pe alloc 3 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$output/intFiles/tophat2/$sample/accepted_hits_RO.bam O=$output/alignments/tophat2/$pre\_$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true RGID=$sample\_1 RGLB=_1 RGPL=Illumina RGPU=$sample\_1 RGSM=$sample`;

	sleep(5);

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_TOPHAT_CRM -hold_jid $pre\_$uID\_TOPHAT_AORRG_$sample -pe alloc 1 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/alignments/tophat2/$pre\_$sample\.bam O=$output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/tophat2/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=NONE METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	push @crm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_TOPHAT_ASM -hold_jid $pre\_$uID\_TOPHAT_AORRG_$sample -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/alignments/tophat2/$pre\_$sample\.bam OUTPUT=$output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	push @asm_tophat, "-metrics $output/intFiles/tophat2/$pre\_$sample\_AlignmentSummaryMetrics.txt";
    }

    my $starOut = '';    
    if($star){
	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_1PASS_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_1PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
	if($pass1){
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SP_$sample -hold_jid $pre\_$uID\_STAR_1PASS_$sample -pe alloc 1 -l virtual_free=2G $Bin/qCMD $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_1PASS_Aligned.out.sam`;
	    
	    $starOut = "$output/intFiles/$sample/$sample\_STAR_1PASS_Aligned.out.sam_filtered.sam";
	}
	else{
	    `/bin/mkdir -m 775 -p $output/intFiles/$sample`;
	    `/bin/mkdir -m 775 -p $output/intFiles/$sample/star2passGG`;
	    sleep(5);
	
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_GG2_$sample -hold_jid $pre\_$uID\_STAR_1PASS_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --runMode genomeGenerate --genomeDir $output/intFiles/$sample/star2passGG --genomeFastaFiles $REF_SEQ --sjdbFileChrStartEnd $output/intFiles/$sample/$sample\_STAR_1PASS_SJ.out.tab --sjdbOverhang $MIN_READ_LENGTH --runThreadN 12`;
	    
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_2PASS_$sample -hold_jid $pre\_$uID\_STAR_GG2_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $output/intFiles/$sample/star2passGG --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $output/intFiles/$sample/$sample\_STAR_2PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
	    
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SP_$sample -hold_jid $pre\_$uID\_STAR_2PASS_$sample -pe alloc 1 -l virtual_free=2G $Bin/qCMD $Bin/starProcessing.pl $output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam`;
	    $starOut = "$output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.sam";
	}

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_MERGE_$sample -hold_jid $pre\_$uID\_SP_$sample -pe alloc 12 -l virtual_free=7G -q lau.q,lcg.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles I=$starOut O=$output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true`;

	sleep(5);

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_AORRG_$sample -hold_jid $pre\_$uID\_STAR_MERGE_$sample -pe alloc 3 -l virtual_free=3G -q lau.q,lcg.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$output/intFiles/$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.bam O=$output/alignments/$pre\_$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true RGID=$sample\_1 RGLB=_1 RGPL=Illumina RGPU=$sample\_1 RGSM=$sample`;

	sleep(5);
	
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_CRM -hold_jid $pre\_$uID\_STAR_AORRG_$sample -pe alloc 1 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx3g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectRnaSeqMetrics I=$output/alignments/$pre\_$sample\.bam O=$output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt CHART_OUTPUT=$output/metrics/$pre\_$sample\_CollectRnaSeqMetrics_chart.pdf REF_FLAT=$REF_FLAT RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS STRAND_SPECIFICITY=NONE METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE`;
	push @crm, "-metrics $output/intFiles/$pre\_$sample\_CollectRnaSeqMetrics.txt";

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_ASM -hold_jid $pre\_$uID\_STAR_AORRG_$sample -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$output/alignments/$pre\_$sample\.bam OUTPUT=$output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	push @asm, "-metrics $output/intFiles/$pre\_$sample\_AlignmentSummaryMetrics.txt";
    }

    if($detectFusions){
	`/bin/mkdir -m 775 -p $output/fusion`;
	if($species =~ /hg19|human/i){	
	    my @fusions = ();
	    my $r1_files = join(" ", @R1);
	    my $r2_files = join(" ", @R2);

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/zcat $r1_files ">$output/intFiles/$sample/$sample\_R1.fastq"`;
	    if($samp_pair{$sample} eq "PE"){
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/zcat $r2_files ">$output/intFiles/$sample/$sample\_R2.fastq"`;
	    }
	
	    if($chimerascan){
		### NOTE: CHIMERASCAN ONLY WORKS FOR PE READS
		if($samp_pair{$sample} eq "PE"){
		    `/bin/mkdir -m 775 -p $output/fusion/chimerascan`;
		    `/bin/mkdir -m 775 -p $output/fusion/chimerascan/$sample`;
		
		    ### NOTE: CHIMERASCAN FAILS WHEN A READ PAIR IS OF DIFFERENT LENGTHS
		    ###       e.g. WHEN WE CLIP AND TRIM VARIABLE LENGTHS FROM
		    ###       SO HAVE TO USE UNPROCESSED READS
		    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CHIMERASCAN_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=1G $Bin/qCMD /opt/bin/python $CHIMERASCAN/chimerascan_run.py -p 6 --quals solexa --multihits=10 --filter-false-pos=$Bin/data/hg19_bodymap_false_positive_chimeras.txt $CHIMERASCAN_INDEX $output/intFiles/$sample/$sample\_R1.fastq $output/intFiles/$sample/$sample\_R2.fastq $output/fusion/chimerascan/$sample/`;
		
		    push @fusions, "--chimerascan $output/fusion/chimerascan/$sample/chimeras.bedpe";
		}
	    }

	    if($star_fusion){
		`/bin/mkdir -m 775 -p $output/fusion/star`;
		`/bin/mkdir -m 775 -p $output/fusion/star/$sample`;
		
		my $inReads = "$output/intFiles/$sample/$sample\_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " $output/intFiles/$sample/$sample\_R2.fastq";
		}

		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_CHIMERA_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 6 --outFileNamePrefix $output/fusion/star/$sample/$sample\_STAR_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --chimSegmentMin 20`;

		push @fusions, "--star $output/fusion/star/$sample/$sample\_STAR_Chimeric.out.junction";
	    }

	    if($mapsplice){
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice`;
		`/bin/mkdir -m 775 -p $output/fusion/mapsplice/$sample`;

		my $inReads = "-1 $output/intFiles/$sample/$sample\_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= " -2 $output/intFiles/$sample/$sample\_R2.fastq"
		}

		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MAPSPLICE_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD /opt/bin/python $MAPSPLICE/mapsplice.py -p 6 --bam --fusion-non-canonical -c $chrSplits -x $BOWTIE_INDEX -o $output/fusion/mapsplice/$sample $inReads --gene-gtf $GTF`;

		push @fusions, "--mapsplice $output/fusion/mapsplice/$sample/fusions_well_annotated.txt";
	    }

	    if($defuse){
		### NOTE: DEFUSE ONLY WORKS FOR PE READS
		if($samp_pair{$sample} eq "PE"){
		    `/bin/mkdir -m 775 -p $output/fusion/defuse`;
		    `/bin/mkdir -m 775 -p $output/fusion/defuse/$sample`;
				
		    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEFUSE_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD $DEFUSE/scripts/defuse.pl --config $DEFUSE/scripts/config.txt --output $output/fusion/defuse/$sample --parallel 6 --1fastq $output/intFiles/$sample/$sample\_R1.fastq --2fastq $output/intFiles/$sample/$sample\_R2.fastq`;
		
		    push @fusions, "--defuse $output/fusion/defuse/$sample/results.filtered.tsv";
		}
	    }

	    if($fusioncatcher){
		`/bin/mkdir -m 775 -p $output/fusion/fusioncatcher`;
		`/bin/mkdir -m 775 -p $output/fusion/fusioncatcher/$sample`;

		my $inReads = "-i $output/intFiles/$sample/$sample\_R1.fastq";	
		if($samp_pair{$sample} eq "PE"){
		    $inReads .= ",$output/intFiles/$sample/$sample\_R2.fastq"
		}

		### NOTE: DUE TO PYTHON MODULE COMPATIBILITY ISSUES ON NODES
		###       HAVE TO USE DIFFERENT VERSION OF PYTHON
		`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_FC_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD_FC $FUSIONCATCHER/bin/fusioncatcher -d $FUSIONCATCHER/data/ensembl_v73 $inReads -o $output/fusion/fusioncatcher/$sample -p 6`;

		push @fusions, "--fusioncatcher $output/fusion/fusioncatcher/$sample/final-list_candidate-fusion-genes.txt";

	    }

	    my $mergeFusions = join(" ", @fusions);
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MERGE_FUSION_$sample -hold_jid $pre\_$uID\_CHIMERASCAN_$sample,$pre\_$uID\_STAR_CHIMERA_$sample,$pre\_$uID\_MAPSPLICE_$sample,$pre\_$uID\_DEFUSE_$sample,$pre\_$uID\_FC_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/MergeFusion $mergeFusions --out $output/fusion/$pre\_merged_fusions_$sample\.txt --normalize_gene $Bin/data/hugo_data_073013.tsv`;
	}
	else{
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSKIPPING FUSION CALLING FOR SAMPLE $sample BECAUSE IT ISN'T SUPPORTED FOR $species; CURRENTLY ONLY SUPPORT FOR hg19";
	}
    }
    if($cufflinks){
	`/bin/mkdir -m 775 -p $output/cufflinks`;
	if($star){
	    `/bin/mkdir -m 775 -p $output/cufflinks/$sample`;
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUFFLINKS_STAR -hold_jid $pre\_$uID\_STAR_MERGE_$sample -pe alloc 5 -l virtual_free=2G $Bin/qCMD $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o $output/cufflinks/$sample $output/alignments/$pre\_$sample\.bam`;
	}
	
	if($tophat){
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2`;
	    `/bin/mkdir -m 775 -p $output/cufflinks/tophat2/$sample`;
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUFFLINKS_TOPHAT -hold_jid $pre\_$uID\_TOPHAT_$sample -pe alloc 5 -l virtual_free=2G $Bin/qCMD $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o $output/cufflinks/tophat2/$sample $output/intFiles/tophat2/$sample/accepted_hits.bam`;
	}
    }
    
    if($htseq || $dexseq){
	if($star){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_QNS_STAR_$sample -hold_jid $pre\_$uID\_SP_$sample -pe alloc 12 -l virtual_free=7G -q lau.q,lcg.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles I=$starOut O=$output/intFiles/$sample/$sample\_STAR_queryname_sorted.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID USE_THREADING=true`;
	}
	
	if($tophat){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_QNS_TOPHAT_$sample -hold_jid $pre\_$uID\_TOPHAT_$sample -pe alloc 12 -l virtual_free=7G -q lau.q,lcg.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles INPUT=$output/intFiles/tophat2/$sample/accepted_hits.bam OUTPUT=$output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam SORT_ORDER=queryname TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT`;
	}
	
	sleep(5);
	
	if($htseq){
	    `/bin/mkdir -m 775 -p $output/counts_gene`;
	    if($star){		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HT_STAR -hold_jid $pre\_$uID\_QNS_STAR_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD "$HTSEQ/htseq-count -m intersection-strict -s no -t exon $output/intFiles/$sample/$sample\_STAR_queryname_sorted.sam $GTF > $output/counts_gene/$sample.htseq_count"`;
	    }
	    
	    if($tophat){
		`/bin/mkdir -m 775 -p $output/counts_gene/tophat2`;
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HT_TOPHAT -hold_jid $pre\_$uID\_QNS_TOPHAT_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD "$HTSEQ/htseq-count -m intersection-strict -s no -t exon $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $GTF > $output/counts_gene/tophat2/$sample.htseq_count"`;
	    }
	}
	
	if($dexseq){
	    `/bin/mkdir -m 775 -p $output/counts_exon`;
	    if($star){		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEX_STAR -hold_jid $pre\_$uID\_QNS_STAR_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $output/intFiles/$sample/$sample\_STAR_queryname_sorted.sam $output/counts_exon/$sample.dexseq_count`;
	    }
	    
	    if($tophat){
		`/bin/mkdir -m 775 -p $output/counts_exon/tophat2`;
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEX_TOPHAT -hold_jid $pre\_$uID\_QNS_TOPHAT_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $output/intFiles/$sample/$sample\_TOPHAT2_accepted_hits_queryname_sorted.sam $output/counts_exon/tophat2/$sample.dexseq_count`;
	    }
	}
    }
}

sleep(5);
if($htseq){
    if($star){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_HTSEQ_STAR -hold_jid $pre\_$uID\_HT_STAR -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_gene '*.htseq_count' $curDir/$output/counts_gene/$pre\_htseq_all_samples.txt $geneNameConversion`;
    }

    if($tophat){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_HTSEQ_TOPHAT -hold_jid $pre\_$uID\_HT_TOPHAT -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_gene/tophat2 '*.htseq_count' $curDir/$output/counts_gene/tophat2/$pre\_htseq_all_samples.txt $geneNameConversion`;
    }
}
if($dexseq){
    if($star){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_DEX_STAR -hold_jid $pre\_$uID\_DEX_STAR -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_exon '*.dexseq_count' $curDir/$output/counts_exon/$pre\_dexseq_all_samples.txt $geneNameConversion`;
    }

    if($tophat){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_DEX_TOPHAT -hold_jid $pre\_$uID\_DEX_TOPHAT -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/$output/counts_exon/tophat2 '*.dexseq_count' $curDir/$output/counts_exon/tophat2/$pre\_dexseq_all_samples.txt $geneNameConversion`;
    }
}

if($star){
    my $crmfiles = join(" ", @crm);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_CRM -hold_jid $pre\_$uID\_STAR_CRM -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/mergePicardMetrics.pl $crmfiles ">$output/metrics/$pre\_CollectRnaSeqMetrics.txt"`;

    my $asmfiles = join(" ", @asm);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_ASM -hold_jid $pre\_$uID\_STAR_ASM -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/mergePicardMetrics.pl $asmfiles ">$output/metrics/$pre\_AlignmentSummaryMetrics.txt"`;
}
if($tophat){
    my $crmfiles_tophat = join(" ", @crm_tophat);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_CRM_TOPHAT -hold_jid $pre\_$uID\_TOPHAT_CRM -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/mergePicardMetrics.pl $crmfiles_tophat ">$output/metrics/tophat2/$pre\_CollectRnaSeqMetrics.txt"`;

    my $asmfiles_tophat = join(" ", @asm_tophat);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_ASM_TOPHAT -hold_jid $pre\_$uID\_TOPHAT_ASM -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/mergePicardMetrics.pl $asmfiles_tophat ">$output/metrics/tophat2$pre\_AlignmentSummaryMetrics.txt"`;
}

### SAMPLE COMMAND
####qsub -N Proj_4226_DESeq ~/bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript ~/RNAseqPipe/trunk/bin/RunDE.R "\"proj.id='4226'\" \"output.dir='/ifs/data/byrne/rnaseq/Proj_4226'\" \"counts.file='Proj_4226_ALL_samples.htseq.count'\" \"key.file='/ifs/data/byrne/rnaseq/Proj_4226/sampleKey.txt'\" \"comps=c('hi - lo')\""


if($deseq){    
    open(COMP, "$comparisons") || die "Can't open comparisons file $comparisons $!";
    my @comps = ();
    while(<COMP>){
	chomp;
	
	my @com = split(/\s+/, $_);
	push @comps, "'$com[0] - $com[1]'";
    }
    close COMP;
    
    my $cmpStr = join(",", @comps);
    
    if($star){
	`/bin/mkdir -m 775 -p $output/differentialExpression_gene`;
	`/bin/mkdir -m 775 -p $output/clustering`;
	`/bin/mkdir -m 775 -p $output/gsa`;
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq_STAR -hold_jid $pre\_$uID\_MATRIX_HTSEQ_STAR -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"species='$species'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"diff.exp.dir='$curDir/$output/differentialExpression_gene'\\\"\" \"\\\"counts.file='$curDir/$output/counts_gene/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"counts.dir='$curDir/$output/counts_gene'\\\"\" \"\\\"clustering.dir='$curDir/$output/clustering'\\\"\" \"\\\"gsa.dir='$curDir/$output/gsa'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
    }
    if($tophat){
	`/bin/mkdir -m 775 -p $output/differentialExpression_gene/tophat2`;
	`/bin/mkdir -m 775 -p $output/clustering/tophat2`;
	`/bin/mkdir -m 775 -p $output/gsa/tophat2`;
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq_TOPHAT -hold_jid $pre\_$uID\_MATRIX_HTSEQ_TOPHAT -pe alloc 1 -l virtual_free=1G $Bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"species='$species'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"diff.exp.dir='$curDir/$output/differentialExpression_gene/tophat2'\\\"\" \"\\\"counts.file='$curDir/$output/counts_gene/tophat2/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"counts.dir='$curDir/$output/counts_gene/tophat2'\\\"\" \"\\\"clustering.dir='$curDir/$output/clustering/tophat2'\\\"\" \"\\\"gsa.dir='$curDir/$output/gsa/tophat2'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
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

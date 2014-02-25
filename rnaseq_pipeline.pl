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

my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq);

$pre = 'TEMP';
GetOptions ('map=s' => \$map,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'samplekey=s' => \$samplekey,
	    'comparisons=s' => \$comparisons,
	    'help' => \$help,
	    'cufflinks' => \$cufflinks,
	    'dexseq' => \$dexseq,
	    'htseq' => \$htseq,
	    'deseq' => \$deseq,
	    'chimerascan' => \$chimerascan,
            'species=s' => \$species);

if(!$map || !$species || $help){
    print <<HELP;

    USAGE: ./rnaseq_pipeline/exome_pipeline.pl -map MAP -species SPECIES -pre PRE -config CONFIG
	* MAP: file listing sample information for processing (REQUIRED)
	* SPECIES: only hg19 and mm9 and human-mouse hybrid (hybrid) currently supported (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (default: /opt/bin)
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B
	* ANALYSES SUPPORTED: cufflinks (-cufflinks), htseq (=htseq), dexseq (-dexseq), deseq (-deseq; must specify samplekey and comparisons), chimerascan (-chimerascan)

HELP
exit;
}

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($species !~ /human|hg19|mouse|mm9|hybrid/i){
    die "Species must be human (hg19) or mouse (mm9) or human-mouse hybrid (hybrid)";
}

my $GTF = '';
my $DEXSEQ_GTF = '';
my $CHIMERASCAN_INDEX = '';
my $geneNameConversion = '';
if($species =~ /human|hg19/i){
    $species = 'hg19';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $CHIMERASCAN_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/chimerascan/';
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
}
elsif($species =~ /mouse|mm9/i){
    $species = 'mm9';
    $GTF = "$Bin/data//Mus_musculus.NCBIM37.67_ENSEMBL.gtf";
    $DEXSEQ_GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf";
    $geneNameConversion = "$Bin/data/mm9Ensembl67IDToGeneName.txt";
}
elsif($species =~ /human-mouse|mouse-human|hybrid/i){
    $species = 'hybrid';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
}


my $CUFFLINKS = '/opt/bin';
my $HTSEQ = '/opt/bin';
my $DEXSEQ = '/opt/bin';
my $PICARD = '/opt/bin';
my $CHIMERASCAN = '/opt/bin';

### Check that all programs are available
&verifyConfig($config);

if($deseq){
    if(!-e $comparisons || !-e $samplekey){
	print "MUST PROVIDE SOMPARISONS AND SAMPLEKEY FILES";
	die;
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

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my %samp_libs_run = ();
my $slr_count = 0;

open(LOG, ">$cd\_rnaseq_pipeline.log") or die "can't write to output log";
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
	print LOG "WARNING: $data[1]\t$data[0]\t$data[2] isn't unique;";
	$data[2] = "$data[2]\_$slr_count";
	print LOG "writing instead to $data[1]\/$data[0]\/$data[2]\n";
    }

    ###if(-e "$data[1]/$data[0]/$data[2]/FINISHED_$data[1]\_$data[0]\_$data[2]"){
	###print LOG "SKIPPING $_ because $data[1]/$data[0]/$data[2]/FINISHED_$data[1]\_$data[0]\_$data[2] exists and therefore considered processed";
	###next;
    ###}

    if(!-d "$data[1]"){
	`mkdir $data[1]`;
    }
    if(!-d "$data[1]/$data[0]"){
	`mkdir $data[1]/$data[0]`;
    }
    if(!-d "$data[1]/$data[0]/$data[2]"){
	`mkdir $data[1]/$data[0]/$data[2]`;
    }

    $samp_libs_run{$data[1]}{$data[0]}{$data[2]} = 1;

    `ln -s $data[3]/* $data[1]/$data[0]/$data[2]/`;

    chdir "$data[1]/$data[0]/$data[2]";

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
		    print OUT "$file\t$file_R2\n";
		}
		elsif($data[4] =~ /se/i){
		    print OUT "$file\n";
		}
	    }
	}
    }
    close OUT;
   

    if($htseq || $dexseq || $cufflinks){
	if($data[4] =~ /pe/i){
	    `$Bin/rnaseq_illumina_PE.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "RGID=$data[1]\_$data[0]\_$data[2] RGPL=Illumina RGPU=$data[1]\_$data[0]\_$data[2] RGLB=$data[1]\_$data[0] RGSM=$data[1]" -species $species -config $config > files_$data[1]\_$data[0]\_$data[2]\_solexa_PE.log 2>&1`;    
	    
	    ###`/common/sge/bin/lx24-amd64/qsub $Bin/qCMD $Bin/rnaseq_illumina_PE.pl -file files -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "RGID=$data[1]\_$data[0]\_$data[2] RGPL=Illumina RGPU=$data[1]\_$data[0]\_$data[2] RGLB=$data[1]\_$data[0] RGSM=$data[1]" -species $species -config $config`;
	}
	elsif($data[4] =~ /se/i){
	    `$Bin/rnaseq_illumina_SE.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "RGID=$data[1]\_$data[0]\_$data[2] RGPL=Illumina RGPU=$data[1]\_$data[0]\_$data[2] RGLB=$data[1]\_$data[0] RGSM=$data[1]" -species $species -config $config > files_$data[1]\_$data[0]\_$data[2]\_solexa_SE.log 2>&1`;
	}
	
    }
    ###`/bin/touch FINISHED_$data[1]\_$data[0]\_$data[2]`;
    chdir $curDir;
}

if($htseq || $dexseq || $cufflinks){
    my $qsubADDRGCheck = `/common/sge/bin/lx24-amd64/qstat -j $pre\_$uID\_ADDRG`;
    while($qsubADDRGCheck){
	
	print LOG "$pre\_$uID\_ADDRG is still running; will check again in an hour...\n";
	sleep(3600);
	
	$qsubADDRGCheck = `/common/sge/bin/lx24-amd64/qstat -j $pre\_$uID\_ADDRG`;
    }
}

foreach my $sample (keys %samp_libs_run){
    my @filteredBams = ();
    my @R1 = ();
    my @R2 = ();
    my $readsFlag = 0;
    foreach my $lib (keys %{$samp_libs_run{$sample}}){
	foreach my $run (keys %{$samp_libs_run{$sample}{$lib}}){
	    if($htseq || $dexseq || $cufflinks){
		if(!-e "$sample/$lib/$run/$sample\_$lib\_$run\_FILTERED.bam"){
		    print LOG  "Can't locate $sample/$lib/$run/$sample\_$lib\_$run\_FILTERED.bam $!";
		    die;
		}
	    }
	    push @filteredBams, "I=$sample/$lib/$run/$sample\_$lib\_$run\_FILTERED.bam";
	    
	    if($chimerascan){
		if(!-d "chimerascan"){
		    `mkdir chimerascan`;
		}

		open(READS, "$sample/$lib/$run/files_$sample\_$lib\_$run") || die "Can't open $sample/$lib/$run/files_$sample\_$lib\_$run $!";
		while(<READS>){
		    chomp;
		    
		    my @readPair = split(/\s+/, $_);
		    if(-e "$sample/$lib/$run/$readPair[0]" && -e "$sample/$lib/$run/$readPair[1]"){
			push @R1, "$sample/$lib/$run/$readPair[0]";
			push @R2, "$sample/$lib/$run/$readPair[1]";
		    }
		    else{
			$readsFlag = 1;
			print "SKIPPING CHIMERASCAN BECAUSE CAN'T LOCATE $sample/$lib/$run/$sample\_$lib\_$run\_R1.fastq && $sample/$lib/$run/$sample\_$lib\_$run\_R2.fastq\n";
		    }
		}
		close READS;
	    }
	}
    }
    
    if($htseq || $dexseq || $cufflinks){
	my $fin = join(" ", @filteredBams);
	if(scalar(@filteredBams) == 1){
	    my @fpath = split(/\//, $filteredBams[0]);
	    my $fshift = shift @fpath;
	    my $fBam = join("/", @fpath);
	    my $fBai = $fBam;
	    $fBai =~ s/\.bam$/\.bai/;

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_$sample -hold_jid $pre\_$uID\_ADDRG -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/ln -s $fBam $sample/$sample\.bam`;
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_$sample -hold_jid $pre\_$uID\_ADDRG -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/ln -s $fBai $sample/$sample\.bai`;
	}
	else{
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_$sample -hold_jid $pre\_$uID\_ADDRG -pe alloc 24 -l virtual_free=3G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar $fin O=$sample/$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=5000000`;
	}
    }

    if($chimerascan && $readsFlag == 0){
	if($species !~ /hg19|human/i){
	    print "CHIMERASCAN ISN'T SUPPORTED FOR $species; CURRENTLY ONLY SUPPORT FOR hg19";
	    next;
	}
	
	my $r1_files = join(" ", @R1);
	my $r2_files = join(" ", @R2);
	
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/cat $r1_files ">$sample/$sample\_R1.fastq.gz"`;
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/cat $r2_files ">$sample/$sample\_R2.fastq.gz"`;
	
	if(!-d "chimerascan/$sample"){
	    `mkdir chimerascan/$sample`;
	}
	
	### NOTE: CHIMERASCAN FAILS WHEN A READ PAIR IS OF DIFFERENT LENGTHS
	###       e.g. WHEN WE CLIP AND TRIM VARIABLE LENGTHS FROM
	###       SO HAVE TO USE UNPROCESSED READS
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CHIMERASCAN_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=7G $Bin/qCMD /opt/bin/python $CHIMERASCAN/chimerascan_run.py -p 6 --quals solexa --multihits=10 --filter-false-pos=/ifs/data/mpirun/data/rnaSeq/hg19_bodymap_false_positive_chimeras.txt $CHIMERASCAN_INDEX $sample/$sample\_R1.fastq.gz $sample/$sample\_R2.fastq.gz chimerascan/$sample/`;
    }
    
    if($cufflinks){
	if(!-d "cufflinks"){
	    `/bin/mkdir cufflinks`;
	}
	if(!-d "cufflinks/$sample"){
	    `/bin/mkdir cufflinks/$sample`;
	}
	
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUFFLINKS -hold_jid $pre\_$uID\_MERGE_$sample -pe alloc 5 -l virtual_free=2G -q lau.q $Bin/qCMD $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o cufflinks/$sample $sample/$sample\.bam`;
    }
    
    if($htseq || $dexseq){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_QNS_$sample -hold_jid $pre\_$uID\_MERGE_$sample -pe alloc 12 -l virtual_free=7G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar INPUT=$sample/$sample\.bam OUTPUT=$sample/$sample\_queryname_sorted.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID USE_THREADING=true`;
    }
    
    sleep(5);
    
    if($htseq){
	if(!-d "htseq"){
	    `mkdir htseq`;
	}
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HT -hold_jid $pre\_$uID\_QNS_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD "$HTSEQ/htseq-count -m intersection-strict -s no -t exon $sample/$sample\_queryname_sorted.sam $GTF > htseq/$sample.htseq_count"`;
    }
    if($dexseq){
	if(!-d "dexseq"){
	    `mkdir dexseq`;
	}
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEX -hold_jid $pre\_$uID\_QNS_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $sample/$sample\_queryname_sorted.sam dexseq/$sample.dexseq_count`;
    }
}

if($htseq){
    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_HTSEQ -hold_jid $pre\_$uID\_HT -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/htseq '*.htseq_count' $curDir/htseq/$pre\_htseq_all_samples.txt $geneNameConversion`;
}
if($dexseq){
    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_DEX -hold_jid $pre\_$uID\_DEX -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/dexseq '*.dexseq_count' $curDir/dexseq/$pre\_dexseq_all_samples.txt $geneNameConversion`;
}

### SAMPLE COMMAND
####qsub -N Proj_4226_DESeq ~/bin/qCMD /opt/R-2.15.0/bin/Rscript ~/RNAseqPipe/trunk/bin/RunDE.R "\"proj.id='4226'\" \"output.dir='/ifs/data/byrne/rnaseq/Proj_4226'\" \"counts.file='Proj_4226_ALL_samples.htseq.count'\" \"key.file='/ifs/data/byrne/rnaseq/Proj_4226/sampleKey.txt'\" \"comps=c('hi - lo')\""


if($deseq){
    if(!-d "DESeq"){
	`mkdir DESeq`;
    }
    
    open(COMP, "$comparisons") || die "Can't open comparisons file $comparisons $!";
    my @comps = ();
    while(<COMP>){
	chomp;
	
	my @com = split(/\s+/, $_);
	push @comps, "'$com[0] - $com[1]'";
    }
    close COMP;
    
    my $cmpStr = join(",", @comps);
    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq -hold_jid $pre\_$uID\_MATRIX_HTSEQ -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/R-2.15.0/bin/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"output.dir='$curDir/DESeq'\\\"\" \"\\\"counts.file='$curDir/htseq/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
}

sub verifyConfig{
    my $paths = shift;

    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
	chomp;
	
	my @conf = split(/\s+/, $_);

	if($conf[0] =~ /picard/i){
	    if(!-e "$conf[1]/MergeSamFiles.jar" || !-e "$conf[1]/AddOrReplaceReadGroups.jar"){
		die "CAN'T FIND MergeSamFiles.jar and/or AddOrReplaceReadGroups.jar IN $conf[1] $!";
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
    }
    close CONFIG;
}

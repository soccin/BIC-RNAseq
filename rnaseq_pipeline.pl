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


my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1);

$pre = 'TEMP';
GetOptions ('map=s' => \$map,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'samplekey=s' => \$samplekey,
	    'comparisons=s' => \$comparisons,
	    'help' => \$help,
	    'star' => \$star,
	    'pass1' => \$pass1,
	    'tophat' => \$tophat,
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
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless -pass1 specified; tophat2 (-tophat)
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons); fusion callers chimerascan (-chimerascan), rna star (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs

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
my $starDB = '';
my $chrSplits = '';
my $BOWTIE_INDEX = '';
my $BOWTIE2_INDEX = '';
my $TRANS_INDEX = '';
my $REF_SEQ = '';

if($species =~ /human|hg19/i){
    $species = 'hg19';
    $REF_SEQ = '/ifs/data/bio/assemblies/H.sapiens/hg19/hg19.fasta';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $CHIMERASCAN_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/chimerascan/';
    $starDB = '/ifs/data/bio/assemblies/H.sapiens/hg19/star';
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
    $BOWTIE_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie/hg19_bowtie';
    $BOWTIE2_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie2/hg19_bowtie2';
    $chrSplits = '/ifs/data/bio/assemblies/H.sapiens/hg19/chromosomes';
    $TRANS_INDEX = '/ifs/data/bio/assemblies/H.sapiens/hg19/bowtie2/transcriptome/gencode18';
}
elsif($species =~ /mouse|mm9/i){
    $species = 'mm9';
    $REF_SEQ = '/ifs/data/bio/assemblies/M.musculus/mm9/mm9.fasta';
    $GTF = "$Bin/data//Mus_musculus.NCBIM37.67_ENSEMBL.gtf";
    $DEXSEQ_GTF = "$Bin/data/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf";
    $starDB = '/ifs/data/bio/assemblies/M.musculus/mm9/star';
    $geneNameConversion = "$Bin/data/mm9Ensembl67IDToGeneName.txt";
    $BOWTIE_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie/mm9_bowtie';
    $BOWTIE2_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie2/mm9_bowtie2';
    $chrSplits = '/ifs/data/bio/assemblies/M.musculus/mm9/chromosomes';
    $TRANS_INDEX = '/ifs/data/bio/assemblies/M.musculus/mm9/bowtie2/transcriptome/ensembl';
}
elsif($species =~ /human-mouse|mouse-human|hybrid/i){
    $species = 'hybrid';
    $REF_SEQ = '/ifs/data/mpirun/genomes/human/hg19_mm9_hybrid/hg19_mm9_hybrid.fasta';
    $GTF = "$Bin/data/gencode.v18.annotation.gtf";
    $DEXSEQ_GTF = "$Bin/data/gencode.v18.annotation_dexseq.gtf";
    $starDB = '/ifs/data/mpirun/genomes/human/hg19_mm9_hybrid/star';
    $geneNameConversion = "$Bin/data/gencode18IDToGeneName.txt";
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
my %samp_pair = ();

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

foreach my $sample (keys %samp_libs_run){
    my @R1 = ();
    my @R2 = ();
    my $readsFlag = 0;
    foreach my $lib (keys %{$samp_libs_run{$sample}}){
	foreach my $run (keys %{$samp_libs_run{$sample}{$lib}}){	    
	    open(READS, "$sample/$lib/$run/files_$sample\_$lib\_$run") || die "Can't open $sample/$lib/$run/files_$sample\_$lib\_$run $!";
	    while(<READS>){
		chomp;
		    
		my @readPair = split(/\s+/, $_);
		
		### NOTE: JUST CHECKING TO SEE IF IT EXISTS
		###       HOWEVER DOES NOT GURANTEE THAT IT'S NON-EMPTY
		if(-e "$sample/$lib/$run/$readPair[0]" && -e "$sample/$lib/$run/$readPair[1]"){
		    push @R1, "$sample/$lib/$run/$readPair[0]";
		    push @R2, "$sample/$lib/$run/$readPair[1]";
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
	if(!-d "tophat2/$sample"){
	    `/bin/mkdir -p tophat2/$sample`;
	}

	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_TOPHAT2_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD $TOPHAT/tophat2 -p 6 --zpacker /opt/pigz-2.1.6/pigz  -r 70 --mate-std-dev 90 --GTF $GTF --transcriptome-index=$TRANS_INDEX -o tophat2/$sample $BOWTIE2_INDEX $inReads`;
    }

    my $starOut = '';    
    if($star){
	my $inReads = "$r1_gz_files";	
	if($samp_pair{$sample} eq "PE"){
	    $inReads .= " $r2_gz_files";
	}

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_1PASS_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $sample/$sample\_STAR_1PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
	if($pass1){
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SP_$sample -hold_jid $pre\_$uID\_STAR_1PASS_$sample -pe alloc 1 -l virtual_free=2G -q lau.q $Bin/qCMD $Bin/starProcessing.pl $sample/$sample\_STAR_1PASS_Aligned.out.sam`;
	    
	    $starOut = "$sample/$sample\_STAR_1PASS_Aligned.out.sam_filtered.sam";
	}
	else{
	    `/bin/mkdir -p $sample/star2passGG`;
	    sleep(5);
	
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_GG2_$sample -hold_jid $pre\_$uID\_STAR_1PASS_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --runMode genomeGenerate --genomeDir $sample/star2passGG --genomeFastaFiles $REF_SEQ --sjdbFileChrStartEnd $sample/$sample\_STAR_1PASS_SJ.out.tab --sjdbOverhang $MIN_READ_LENGTH --runThreadN 12`;
	    
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_2PASS_$sample -hold_jid $pre\_$uID\_STAR_GG2_$sample -pe alloc 12 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $sample/star2passGG --readFilesIn $inReads --runThreadN 12 --outFileNamePrefix $sample/$sample\_STAR_2PASS_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --readFilesCommand zcat`;
	    
	    sleep(5);
	    
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SP_$sample -hold_jid $pre\_$uID\_STAR_2PASS_$sample -pe alloc 1 -l virtual_free=2G -q lau.q $Bin/qCMD $Bin/starProcessing.pl $sample/$sample\_STAR_2PASS_Aligned.out.sam`;
	    $starOut = "$sample/$sample\_STAR_2PASS_Aligned.out.sam_filtered.sam";
	}

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_MERGE_$sample -hold_jid $pre\_$uID\_SP_$sample -pe alloc 12 -l virtual_free=7G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar I=$starOut O=$sample/$sample\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true`;
    }

    if($detectFusions){
	if($species !~ /hg19|human/i){
	    print "FUSION CALLING ISN'T SUPPORTED FOR $species; CURRENTLY ONLY SUPPORT FOR hg19";
	    next;
	}
	
	my @fusions = ();
	my $r1_files = join(" ", @R1);
	my $r2_files = join(" ", @R2);

	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/zcat $r1_files ">$sample/$sample\_R1.fastq"`;
	if($samp_pair{$sample} eq "PE"){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CAT_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /bin/zcat $r2_files ">$sample/$sample\_R2.fastq"`;
	}
	
	if($chimerascan){
	    ### NOTE: CHIMERASCAN ONLY WORKS FOR PE READS
	    if($samp_pair{$sample} eq "PE"){
		if(!-d "fusion/chimerascan/$sample"){
		    `/bin/mkdir -p fusion/chimerascan/$sample`;
		}
		
		### NOTE: CHIMERASCAN FAILS WHEN A READ PAIR IS OF DIFFERENT LENGTHS
		###       e.g. WHEN WE CLIP AND TRIM VARIABLE LENGTHS FROM
		###       SO HAVE TO USE UNPROCESSED READS
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CHIMERASCAN_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=1G $Bin/qCMD /opt/bin/python $CHIMERASCAN/chimerascan_run.py -p 6 --quals solexa --multihits=10 --filter-false-pos=$Bin/data/hg19_bodymap_false_positive_chimeras.txt $CHIMERASCAN_INDEX $sample/$sample\_R1.fastq $sample/$sample\_R2.fastq fusion/chimerascan/$sample/`;
		
		push @fusions, "--chimerascan fusion/chimerascan/$sample/chimeras.bedpe";
	    }
	}

	if($star_fusion){
	    if(!-d "fusion/star/$sample"){
		`/bin/mkdir -p fusion/star/$sample`;
	    }

	    my $inReads = "$sample/$sample\_R1.fastq";	
	    if($samp_pair{$sample} eq "PE"){
		$inReads .= " $sample/$sample\_R2.fastq";
	    }

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_STAR_CHIMERA_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=3G $Bin/qCMD $STAR/STAR --genomeDir $starDB --readFilesIn $inReads --runThreadN 6 --outFileNamePrefix fusion/star/$sample/$sample\_STAR_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMattributes All --outSAMunmapped Within --chimSegmentMin 20`;

	    push @fusions, "--star fusion/star/$sample/$sample\_STAR_Chimeric.out.junction";
	}

	if($mapsplice){
	    if(!-d "fusion/mapsplice/$sample"){
		`/bin/mkdir -p fusion/mapsplice/$sample`;
	    }

	    my $inReads = "-1 $sample/$sample\_R1.fastq";	
	    if($samp_pair{$sample} eq "PE"){
		$inReads .= " -2 $sample/$sample\_R2.fastq"
	    }

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MAPSPLICE_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD /opt/bin/python $MAPSPLICE/mapsplice.py -p 6 --bam --fusion-non-canonical -c $chrSplits -x $BOWTIE_INDEX -o fusion/mapsplice/$sample $inReads --gene-gtf $GTF`;

	    push @fusions, "--mapsplice fusion/mapsplice/$sample/fusions_well_annotated.txt";
	}

	if($defuse){
	    ### NOTE: DEFUSE ONLY WORKS FOR PE READS
	    if($samp_pair{$sample} eq "PE"){
		if(!-d "fusion/defuse/$sample"){
		    `/bin/mkdir -p fusion/defuse/$sample`;
		}
				
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEFUSE_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD $DEFUSE/scripts/defuse.pl --config $DEFUSE/scripts/config.txt --output fusion/defuse/$sample --parallel 6 --1fastq $sample/$sample\_R1.fastq --2fastq $sample/$sample\_R2.fastq`;
		
		push @fusions, "--defuse fusion/defuse/$sample/results.filtered.tsv";
	    }
	}

	if($fusioncatcher){
	    if(!-d "fusion/fusioncatcher/$sample"){
		`/bin/mkdir -p fusion/fusioncatcher/$sample`;
	    }

	    my $inReads = "-i $sample/$sample\_R1.fastq";	
	    if($samp_pair{$sample} eq "PE"){
		$inReads .= ",$sample/$sample\_R2.fastq"
	    }

	    ### NOTE: DUE TO PYTHON MODULE COMPATIBILITY ISSUES ON NODES
	    ###       HAVE TO USE DIFFERENT VERSION OF PYTHON
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_FC_$sample -hold_jid $pre\_$uID\_CAT_$sample -pe alloc 6 -l virtual_free=2G $Bin/qCMD_FC $FUSIONCATCHER/bin/fusioncatcher -d $FUSIONCATCHER/data/ensembl_v73 $inReads -o fusion/fusioncatcher/$sample -p 6`;

	    push @fusions, "--fusioncatcher fusion/fusioncatcher/$sample/final-list_candidate-fusion-genes.txt";

	}

	my $mergeFusions = join(" ", @fusions);
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MERGE_FUSION_$sample -hold_jid $pre\_$uID\_CHIMERASCAN_$sample,$pre\_$uID\_STAR_CHIMERA_$sample,$pre\_$uID\_MAPSPLICE_$sample,$pre\_$uID\_DEFUSE_$sample,$pre\_$uID\_FC_$sample -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/MergeFusion $mergeFusions --out $pre\_merged_fusions_$sample\.txt --normalize_gene $Bin/data/hugo_data_073013.tsv`;
    }
    
    if($cufflinks){
	if($star){
	    if(!-d "cufflinks/$sample"){
		`/bin/mkdir -p cufflinks/$sample`;
	    }
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUFFLINKS_STAR -hold_jid $pre\_$uID\_STAR_MERGE_$sample -pe alloc 5 -l virtual_free=2G -q lau.q $Bin/qCMD $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o cufflinks/$sample $sample/$sample\.bam`;
	}
	
	if($tophat){
	    if(!-d "cufflinks_tophat2/$sample"){
		`/bin/mkdir -p cufflinks_tophat2/$sample`;
	    }
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUFFLINKS_TOPHAT2 -hold_jid $pre\_$uID\_TOPHAT2_$sample -pe alloc 5 -l virtual_free=2G -q lau.q $Bin/qCMD $CUFFLINKS/cufflinks -q -p 12 --no-update-check -N -G $GTF -o cufflinks/$sample $sample/tophat2/accepted_hits.bam`;
	}
    }
    
    if($htseq || $dexseq){
	if($star){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_QNS_STAR_$sample -hold_jid $pre\_$uID\_SP_$sample -pe alloc 12 -l virtual_free=7G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar I=$starOut O=$sample/$sample\_queryname_sorted.sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true`;
	}
	
	if($tophat){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_QNS_TOPHAT2_$sample -hold_jid $pre\_$uID\_TOPHAT2_$sample -pe alloc 12 -l virtual_free=7G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar INPUT=$sample/tophat2/accepted_hits.bam OUTPUT=$sample/tophat2/accepted_hits_queryname_sorted.sam SORT_ORDER=queryname TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT`;    
	}
	
	sleep(5);
	
	if($htseq){
	    if($star){
		if(!-d "htseq"){
		    `/bin/mkdir -p htseq`;
		}
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HT_STAR -hold_jid $pre\_$uID\_QNS_STAR_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD "$HTSEQ/htseq-count -m intersection-strict -s no -t exon $sample/$sample\_queryname_sorted.sam $GTF > htseq/$sample.htseq_count"`;
	    }
	    
	    if($tophat){
		if(!-d "htseq_tophat2"){
		    `/bin/mkdir -p htseq_tophat2`;
		}
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HT_TOPHAT2 -hold_jid $pre\_$uID\_QNS_TOPHAT2_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD "$HTSEQ/htseq-count -m intersection-strict -s no -t exon $sample/tophat2/accepted_hits_queryname_sorted.sam $GTF > htseq_tophat2/$sample.htseq_count"`;
	    }
	}
	
	if($dexseq){
	    if($star){
		if(!-d "dexseq"){
		    `/bin/mkdir -p dexseq`;
		}
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEX_STAR -hold_jid $pre\_$uID\_QNS_STAR_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $sample/$sample\_queryname_sorted.sam dexseq/$sample.dexseq_count`;
	    }
	    
	    if($tophat){
		if(!-d "dexseq_tophat2"){
		    `/bin/mkdir -p dexseq_tophat2`;
		}
		
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_DEX_TOPHAT2 -hold_jid $pre\_$uID\_QNS_TOPHAT2_$sample -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $DEXSEQ/dexseq_count.py -s no $DEXSEQ_GTF $sample/tophat2/accepted_hits_queryname_sorted.sam dexseq_tophat2/$sample.dexseq_count`;
	    }
	}
    }
}

if($htseq){
    if($star){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_HTSEQ_STAR -hold_jid $pre\_$uID\_HT_STAR -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/htseq '*.htseq_count' $curDir/htseq/$pre\_htseq_all_samples.txt $geneNameConversion`;
    }

    if($tophat){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_HTSEQ_TOPHAT2 -hold_jid $pre\_$uID\_HT_TOPHAT2 -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/htseq_tophat2 '*.htseq_count' $curDir/htseq_tophat2/$pre\_htseq_all_samples.txt $geneNameConversion`;
    }
}
if($dexseq){
    if($star){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_DEX_STAR -hold_jid $pre\_$uID\_DEX_STAR -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/dexseq '*.dexseq_count' $curDir/dexseq/$pre\_dexseq_all_samples.txt $geneNameConversion`;
    }

    if($tophat){
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_MATRIX_DEX_TOPHAT2 -hold_jid $pre\_$uID\_DEX_TOPHAT2 -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/bin/python $Bin/rnaseq_count_matrix.py $curDir/dexseq_tophat2 '*.dexseq_count' $curDir/dexseq_tophat2/$pre\_dexseq_all_samples.txt $geneNameConversion`;
    }
}
### SAMPLE COMMAND
####qsub -N Proj_4226_DESeq ~/bin/qCMD /opt/R-2.15.0/bin/Rscript ~/RNAseqPipe/trunk/bin/RunDE.R "\"proj.id='4226'\" \"output.dir='/ifs/data/byrne/rnaseq/Proj_4226'\" \"counts.file='Proj_4226_ALL_samples.htseq.count'\" \"key.file='/ifs/data/byrne/rnaseq/Proj_4226/sampleKey.txt'\" \"comps=c('hi - lo')\""


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
	if(!-d "DESeq"){
	    `/bin/mkdir -p DESeq`;
	}
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq_STAR -hold_jid $pre\_$uID\_MATRIX_HTSEQ_STAR -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/R-2.15.0/bin/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"output.dir='$curDir/DESeq'\\\"\" \"\\\"counts.file='$curDir/htseq/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
    }
    if($tophat){
	if(!-d "DESeq_tophat2"){
	    `/bin/mkdir -p DESeq_tophat2`;
	}
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_DESeq_TOPHAT2 -hold_jid $pre\_$uID\_MATRIX_HTSEQ_TOPHAT2 -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/R-2.15.0/bin/Rscript $Bin/RunDE.R \"\\\"bin='$Bin'\\\"\" \"\\\"proj.id='$pre'\\\"\" \"\\\"output.dir='$curDir/DESeq_tophat2'\\\"\" \"\\\"counts.file='$curDir/htseq_tophat2/$pre\_htseq_all_samples.txt'\\\"\" \"\\\"key.file='$samplekey'\\\"\" \"\\\"comps=c($cmpStr)\\\"\"`;
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

#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use Schedule;
use Cluster;
use File::Basename;

my ($config, $sample, $bam, $bed, $intdir, $outdir, $progdir, $scheduler, $readlen, $layout, $forceall, $sync, $help);

my $priority_project = "ngs";
my $priority_group = "Pipeline";
my $pre = "TEMP";

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

GetOptions ('pre=s' => \$pre,
            'config=s' => \$config,
            'bam=s' => \$bam,
            'bed=s' => \$bed,
            'sample=s' => \$sample,
            'readlen=i' => \$readlen,
            'layout=s' => \$layout,
            'outdir=s' => \$outdir,
            'intdir=s' => \$intdir,
            'forceall' => \$forceall,
            'sync' => \$sync,
            'progdir=s' => \$progdir,
            'scheduler=s' => \$scheduler,
            'priority_project=s' => \$priority_project,
            'priority_group=s' => \$priority_group,
            'help' => \$help) or exit(1);

if(!$config || !$sample || !$bam || !$outdir || !$bed || !$intdir || !$progdir || !$scheduler){
print <<HELP;

    USAGE: rseqc.pl -config CONFIG -bam BAM -sample SAMPLE -intdir INTDIR -outdir OUTDIR -progdir PROGDIR -scheduler SCHEDULER
        *           CONFIG:  file containing configuration, specifically path to python
        *              BAM:  bam file contianing final alignments for one sample
        *           SAMPLE:  sample name (e.g., s_SAMPLE_1)
        *          READLEN:  read length
        *           LAYOUT:  experiment layout (SE or PE)
        *              BED:  BED file used for junction/splicing information
        *           INTDIR:  intermediate file directory
        *           OUTDIR:  final output directory (metrics/images)
        *          PROGDIR:  progress directory
        *        SCHEDULER:  LSF or SGE
        *         FORCEALL:  force all modules to run even if they have been run before
        *             SYNC:  print string of all jobs submitted for syncing purposes
        * PRIORITY_PROJECT:  sge notion of priority assigned to projects (default: ngs)
        *   PRIORITY_GROUP:  lsf notion of priority assigned to groups (default: Pipeline);
HELP
exit;
}

#added by Julia 07/22/2019 ###
my $SINGULARITY = '';
my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";
#added by Julia 07/22/2019 end ###
my $PYTHON = '';

open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /^python$/i){
        if(!-e "$conf[1]/python"){
        #    die "CAN'T FIND PYTHON IN $conf[1] $!";
        }
        $PYTHON = $conf[1];
    }
    elsif($conf[0] =~ /^singularity$/i){
        if(!-e "$conf[1]/singularity"){
        #    die "CAN'T FIND singularity IN $conf[1] $!";
        }
        $SINGULARITY = $conf[1];
    }
    elsif($conf[0] =~ /singularity_r/i){
        if(!-e "$conf[1]/R"){
        #    die "CAN'T FIND R IN $conf[1] $!";
        }
        #### config now contains new R, but this needs to run with old R
        $singularityenv_prepend_path = "$conf[1]:$singularityenv_prepend_path";
    }
}

my $root_bin = dirname($Bin);
my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$root_bin/rnaseq_pipeline_singularity_prod.simg");
$singularityParams = Schedule::singularityParams(%sinParams);
$singularityBind = Schedule::singularityBind($scheduler);

$ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
$ENV{'SINGULARITY_BINDPATH'} = $singularityBind;

$ENV{'LD_LIBRARY_PATH'} = "/opt/common/CentOS_6/gcc/gcc-4.9.3/lib64:$ENV{'LD_LIBRARY_PATH'}";
close CONFIG;

my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

my $ran_bs = 0;
my $ran_mp = 0;
my $ran_rdp = 0;
my $ran_rdst = 0;
my $ran_rgc = 0;
my $ran_rn = 0;
my $ran_rq = 0;
my $ran_ip = 0;
my $ran_ja = 0;
my $ran_js = 0;
my $ran_ie = 0;
my $ran_id = 0;
my $ran_dp = 0;
my $ran_cp = 0;

my @rseqc_jids = ();
my $file_pre = "rseqc_$sample";

#my $curDir = `pwd`;
#chomp $curDir;
#my $cd = $curDir;
#$cd =~ s/\//_/g;
#open(LOG, ">$cd\_$sample\_rseqc.log") or die "can't write to output log";
#my $alignments_found = `grep $sample Alignment_counts.txt | cut -f 5`;
#if($alignments_found == 0){
#    my @currentTime = &getTime();
#    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tNO ALIGNMENTS FOUND FOR SAMPLE $sample. NO RSEQC STATS.\n"; 
#    exit(0);
#}

### bam stats
if(!-e "$progdir/$pre\_$uID\_RSEQC_BS_$sample.done" || $forceall){
    print "Running bam_stat.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_BS_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_BS_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/bam_stat.py -i $bam ">$intdir/$file_pre\_bam_stat.txt"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_BS_$sample.done`;
    
    push @rseqc_jids, "$pre\_$uID\_RSEQC_BS_$sample";
    $ran_bs = 1;
}

### mismatch profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_MP_$sample.done" || $forceall){
    print "Running mismatch_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_MP_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_MP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/mismatch_profile.py -i $bam -o "$intdir/$file_pre" -l $readlen`; 
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_MP_$sample.done`;
    $ran_mp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_MP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_MP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_MP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.mismatch_profile.pdf && mv $intdir/$file_pre.mismatch_profile.pdf $outdir || true"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_MP_MV_$sample";
}

### read duplication
if(!-e "$progdir/$pre\_$uID\_RSEQC_RDP_$sample.done" || $forceall){
    print "Running read_duplication.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDP_$sample", cpu => "1", mem => "60", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);   
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/read_duplication.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RDP_$sample.done`;
    $ran_rdp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RDP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.DupRate_plot.pdf && mv $intdir/$file_pre.DupRate_plot.pdf $outdir"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RDP_MV_$sample";
}

### read distribution
if(!-e "$progdir/$pre\_$uID\_RSEQC_RDST_$sample.done" || $forceall){
    print "Running read_distribution.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDST_$sample", cpu => "1", mem => "45", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDST_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/read_distribution.py -i $bam -r $bed ">$intdir/$file_pre\_read_distribution.txt"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RDST_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RDST_$sample";
    $ran_rdst = 1;
}

### read GC
if(!-e "$progdir/$pre\_$uID\_RSEQC_RGC_$sample.done" || $forceall){
    print "Running read_GC.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RGC_$sample", cpu => "1", mem => "50", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RGC_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/read_GC.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RGC_$sample.done`;
    $ran_rgc = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RGC_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RGC_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RGC_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.GC_plot.pdf &&  mv $intdir/$file_pre.GC_plot.pdf $outdir"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RGC_MV_$sample";
}

### read NVC
if(!-e "$progdir/$pre\_$uID\_RSEQC_NVC_$sample.done" || $forceall){
    print "Running read_NVC.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_NVC_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_NVC_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/read_NVC.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_NVC_$sample.done`;
    $ran_rn = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_NVC_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_NVC_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_NVC_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.NVC_plot.pdf && mv $intdir/$file_pre.NVC_plot.pdf $outdir"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_NVC_MV_$sample";
}

### read quality
if(!-e "$progdir/$pre\_$uID\_RSEQC_RQ_$sample.done" || $forceall){
    print "Running read_quality.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RQ_$sample", cpu => "1", mem => "320", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RQ_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    #`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/read_quality.py -i $bam -o "$intdir/$file_pre"`;
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} -W 30000 $singularityParams $PYTHON/read_quality.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RQ_$sample.done`;
    $ran_rq = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RQ_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RQ_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RQ_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams mv $intdir/$file_pre.qual.*.pdf $outdir || true"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RQ_MV_$sample";
}

### insertion profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_IP_$sample.done" || $forceall){
    print "Running insertion_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IP_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/insertion_profile.py -i $bam -o "$intdir/$file_pre" -s $layout` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_IP_$sample.done`;
    $ran_ip = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_IP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams mv $intdir/$file_pre.insertion_profile*.pdf $outdir || true"`; 
    push @rseqc_jids, "$pre\_$uID\_RSEQC_IP_MV_$sample";
}

### deletion profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_DP_$sample.done" || $forceall){
    print "Running deletion_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_DP_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_DP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/deletion_profile.py -i $bam -o "$intdir/$file_pre" -l 50` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_DP_$sample.done`;
    $ran_dp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_DP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_DP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_DP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams mv $intdir/$file_pre.deletion_profile*.pdf $outdir || true"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_DP_MV_$sample";
}

### clipping profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_CP_$sample.done" || $forceall){
    print "Running clipping_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_CP_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_CP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/clipping_profile.py -i $bam -o "$intdir/$file_pre" -s $layout` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_CP_$sample.done`;
    $ran_cp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_CP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_CP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_CP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams mv $intdir/$file_pre.clipping_profile*.pdf $outdir || true"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_CP_MV_$sample";
}

### junction annotation
if(!-e "$progdir/$pre\_$uID\_RSEQC_JA_$sample.done" || $forceall){
    print "Running junction_annotation.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JA_$sample", cpu => "1", mem => "10", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JA_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/junction_annotation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_JA_$sample.done`;
    $ran_ja = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JA_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_JA_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JA_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams mv $intdir/$file_pre.splice_*.pdf $outdir || true"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_JA_MV_$sample";
}

### junction saturation
if(!-e "$progdir/$pre\_$uID\_RSEQC_JS_$sample.done" || $forceall){
    print "Running junction_saturation.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JS_$sample", cpu => "1", mem => "80", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JS_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/junction_saturation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_JS_$sample.done`;
    $ran_js = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JS_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_JS_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_Js_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.junctionSaturation_plot.pdf && mv $intdir/$file_pre.junctionSaturation_plot.pdf $outdir"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_JS_MV_$sample";
}

### infer experiment
if(!-e "$progdir/$pre\_$uID\_RSEQC_IE_$sample.done" || $forceall){
    print "Running infer_experiment.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IE_$sample", cpu => "1", mem => "10", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IE_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/infer_experiment.py -i $bam ">$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_IE_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_IE_$sample";
    $ran_ie = 1;
}

### inner distance
if(!-e "$progdir/$pre\_$uID\_RSEQC_ID_$sample.done" || $forceall){
    print "Running inner_distance.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_ID_$sample", cpu => "1", mem => "30", cluster_out => "$progdir/$pre\_$uID\_RSEQC_ID_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/inner_distance.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_ID_$sample.done`;
    $ran_id = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_ID_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_ID_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_ID_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams test -e $intdir/$file_pre.inner_distance_plot.pdf && mv $intdir/$file_pre.inner_distance_plot.pdf $outdir"`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_ID_MV_$sample";
 
}

if($sync){
    foreach my $rseqc_jid (@rseqc_jids){
        `$Bin/jobSync $scheduler $rseqc_jid`;
    } 
}

`sleep 10`;
exit(0);


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


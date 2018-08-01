#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use Schedule;
use Cluster;

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

my $PYTHON = '';
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /^python$/i){
        if(!-e "$conf[1]/python"){
            die "CAN'T FIND PYTHON IN $conf[1] $!";
        }
        $PYTHON = $conf[1];
    }
}
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

### bam stats
if(!-e "$progdir/$pre\_$uID\_RSEQC_BS_$sample.done" || $forceall){
    print "Running bam_stat.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_BS_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_BS_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/bam_stat.py -i $bam ">$intdir/$file_pre\_bam_stat.txt"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_BS_$sample.done`;
    
    push @rseqc_jids, "$pre\_$uID\_RSEQC_BS_$sample";
    $ran_bs = 1;
}

### mismatch profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_MP_$sample.done" || $forceall){
    print "Running mismatch_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_MP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_MP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/mismatch_profile.py -i $bam -o "$intdir/$file_pre" -l $readlen`; 
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_MP_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_MP_$sample";
    $ran_mp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_MP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_MP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_MP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.mismatch_profile.pdf && mv $intdir/$file_pre.mismatch_profile.pdf $outdir\"`;
}

### read duplication
if(!-e "$progdir/$pre\_$uID\_RSEQC_RDP_$sample.done" || $forceall){
    print "Running read_duplication.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/read_duplication.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RDP_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RDP_$sample";
    $ran_rdp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RDP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.DupRate_plot.pdf && mv $intdir/$file_pre.DupRate_plot.pdf $outdir\"`;
}

### read distribution
if(!-e "$progdir/$pre\_$uID\_RSEQC_RDST_$sample.done" || $forceall){
    print "Running read_distribution.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RDST_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RDST_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/read_distribution.py -i $bam -r $bed ">$intdir/$file_pre\_read_distribution.txt"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RDST_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RDST_$sample";
    $ran_rdst = 1;
}

### read GC
if(!-e "$progdir/$pre\_$uID\_RSEQC_RGC_$sample.done" || $forceall){
    print "Running read_GC.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RGC_$sample", cpu => "1", mem => "4", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RGC_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/read_GC.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RGC_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RGC_$sample";
    $ran_rgc = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RGC_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RGC_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RGC_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.GC_plot.pdf &&  mv $intdir/$file_pre.GC_plot.pdf $outdir\"`;
}

### read NVC
if(!-e "$progdir/$pre\_$uID\_RSEQC_NVC_$sample.done" || $forceall){
    print "Running read_NVC.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_NVC_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_NVC_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);    
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/read_NVC.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_NVC_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_NVC_$sample";
    $ran_rn = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_NVC_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_NVC_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_NVC_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.NVC_plot.pdf && mv $intdir/$file_pre.NVC_plot.pdf $outdir\"`;
}

### read quality
if(!-e "$progdir/$pre\_$uID\_RSEQC_RQ_$sample.done" || $forceall){
    print "Running read_quality.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RQ_$sample", cpu => "1", mem => "300", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RQ_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/read_quality.py -i $bam -o "$intdir/$file_pre"`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_RQ_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_RQ_$sample";
    $ran_rq = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_RQ_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_RQ_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_RQ_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams 'for f in $intdir/$file_pre.qual.*.pdf\; do  [ -e \"\$f\" ] && mv \"\$f\" $outdir\; done\;'`;
}

### insertion profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_IP_$sample.done" || $forceall){
    print "Running insertion_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/insertion_profile.py -i $bam -o "$intdir/$file_pre" -s $layout` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_IP_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_IP_$sample";
    $ran_ip = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_IP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams 'for f in $intdir/$file_pre.insertion_profile*.pdf\; do [ -e \"\$f\" ] && mv \"\$f\" $outdir\; done\;'`;
}

### deletion profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_DP_$sample.done" || $forceall){
    print "Running deletion_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_DP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_DP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/deletion_profile.py -i $bam -o "$intdir/$file_pre" -l 50` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_DP_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_DP_$sample";
    $ran_dp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_DP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_DP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_DP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.deletion_profile.pdf && mv $intdir/$file_pre.deletion_profile.pdf $outdir\"`;
}

### clipping profile
if(!-e "$progdir/$pre\_$uID\_RSEQC_CP_$sample.done" || $forceall){
    print "Running clipping_profile.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_CP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_CP_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/clipping_profile.py -i $bam -o "$intdir/$file_pre" -s $layout` ;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_CP_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_CP_$sample";
    $ran_cp = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_CP_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_CP_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_CP_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams 'for f in $intdir/$file_pre.clipping_profile*.pdf\; do [ -e \"\$f\" ] && mv \"\$f\" $outdir\; done\;'`;
}

### junction annotation
if(!-e "$progdir/$pre\_$uID\_RSEQC_JA_$sample.done" || $forceall){
    print "Running junction_annotation.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JA_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JA_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/junction_annotation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_JA_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_JA_$sample";
    $ran_ja = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JA_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_JA_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JA_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams 'for f in $intdir/$file_pre.splice_*.pdf\; do [ -e \"\$f\" ] &&  mv \"\$f\" $outdir\; done\;'`;
}

### junction saturation
if(!-e "$progdir/$pre\_$uID\_RSEQC_JS_$sample.done" || $forceall){
    print "Running junction_saturation.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JS_$sample", cpu => "1", mem => "8", cluster_out => "$progdir/$pre\_$uID\_RSEQC_JS_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/junction_saturation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_JS_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_JS_$sample";
    $ran_js = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_JS_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_JS_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_Js_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.junctionSaturation_plot.pdf && mv $intdir/$file_pre.junctionSaturation_plot.pdf $outdir\"`;
}

### infer experiment
if(!-e "$progdir/$pre\_$uID\_RSEQC_IE_$sample.done" || $forceall){
    print "Running infer_experiment.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_IE_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_IE_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/infer_experiment.py -i $bam ">$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_IE_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_IE_$sample";
    $ran_ie = 1;
}

### inner distance
if(!-e "$progdir/$pre\_$uID\_RSEQC_ID_$sample.done" || $forceall){
    print "Running inner_distance.py\n";
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_ID_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_ID_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/inner_distance.py -i $bam -o "$intdir/$file_pre" -r $bed`;
    `/bin/touch $progdir/$pre\_$uID\_RSEQC_ID_$sample.done`;
    push @rseqc_jids, "$pre\_$uID\_RSEQC_ID_$sample";
    $ran_id = 1;

    ## move only the output pdf file to delivered results
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSEQC_ID_MV_$sample", job_hold => "$pre\_$uID\_RSEQC_ID_$sample", cpu => "1", mem => "1", cluster_out => "$progdir/$pre\_$uID\_RSEQC_ID_MV_$sample.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams \"test -e $intdir/$file_pre.inner_distance_plot.pdf && mv $intdir/$file_pre.inner_distance_plot.pdf $outdir\"`;
}

if($sync){
    foreach my $rseqc_jid (@rseqc_jids){
        `$Bin/jobSync $scheduler $rseqc_jid`;
    } 
}

`sleep 10`;
exit(0);

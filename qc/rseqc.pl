#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
#use Schedule;
#use Cluster;
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
            die "CAN'T FIND PYTHON IN $conf[1] $!";
        }
        $PYTHON = $conf[1];
    }
    elsif($conf[0] =~ /singularity/i){
        if(!-e "$conf[1]/singularity"){
            die "CAN'T FIND singularity IN $conf[1] $!";
        }
        $SINGULARITY = $conf[1];
    }
    elsif($conf[0] =~ /^r$/i){
        if(!-e "$conf[1]/R"){
            die "CAN'T FIND R IN $conf[1] $!";
        }
        $singularityenv_prepend_path .= ":$conf[1]";
    }
}

my $root_bin = dirname($Bin);
##added by Julia on 07/22/2019 ###
#my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$root_bin/rnaseq_pipeline_singularity_prod.simg");
#$singularityParams = Schedule::singularityParams(%sinParams);
#$singularityBind = Schedule::singularityBind($scheduler);

#$ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
#$ENV{'SINGULARITY_BINDPATH'} = $singularityBind;
##added by Julia on 07/22/2019 ended## 
close CONFIG;


#my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
#my $additionalParams = Schedule::additionalParams(%addParams);

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

my @outFiles = ();

print "Running bam_stat.py\n";
`$PYTHON/bam_stat.py -i $bam ">$intdir/$file_pre\_bam_stat.txt"`;
   
print "Running mismatch_profile.py\n";
`$PYTHON/mismatch_profile.py -i $bam -o "$intdir/$file_pre" -l $readlen`;
`test -e "$intdir/$file_pre.mismatch_profile.pdf" && mv "$intdir/$file_pre.mismatch_profile.pdf" $outdir`;

print "Running read_duplication.py\n";
`read_duplication.py -i $bam -o "$intdir/$file_pre"`;
`test -e "$intdir/$file_pre.DupRate_plot.pdf" && mv "$intdir/$file_pre.DupRate_plot.pdf" $outdir`;

print "Running read_distribution.py\n";
`$PYTHON/read_distribution.py -i $bam -r $bed ">$intdir/$file_pre\_read_distribution.txt"`;

print "Running read_GC.py\n";
`$PYTHON/read_GC.py -i $bam -o "$intdir/$file_pre"`;
`test -e "$intdir/$file_pre.GC_plot.pdf" &&  mv "$intdir/$file_pre.GC_plot.pdf" $outdir`;

print "Running read_NVC.py\n";
`$PYTHON/read_NVC.py -i $bam -o "$intdir/$file_pre"`;
`test -e "$intdir/$file_pre.NVC_plot.pdf" && mv "$intdir/$file_pre.NVC_plot.pdf" $outdir`;

print "Running read_quality.py\n";
`$PYTHON/read_quality.py -i $bam -o "$intdir/$file_pre"`;
@outFiles = glob "$intdir/$file_pre.qual.*.pdf";
foreach my $fl (@outFiles) {
    `mv $fl $outdir`;
}

print "Running insertion_profile.py\n";
`$PYTHON/insertion_profile.py -i $bam -o "$intdir/$file_pre" -s $layout`;
@outFiles = glob "$intdir/$file_pre.insertion_profile*.pdf";
foreach my $fl (@outFiles) {
    `mv $fl $outdir`;
}

print "Running deletion_profile.py\n";
`$PYTHON/deletion_profile.py -i $bam -o "$intdir/$file_pre" -l 50` ;
@outFiles = glob "$intdir/$file_pre.deletion_profile*.pdf";
foreach my $fl (@outFiles) {
    `mv $fl $outdir`;
}

print "Running clipping_profile.py\n";
`$PYTHON/clipping_profile.py -i $bam -o "$intdir/$file_pre" -s $layout` ;
@outFiles = "$intdir/$file_pre.clipping_profile*.pdf";
foreach my $fl (@outFiles) {
    `mv $fl $outdir`;
}

print "Running junction_annotation.py\n";
`$PYTHON/junction_annotation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
@outFiles = "$intdir/$file_pre.splice_*.pdf";
foreach my $fl (@outFiles) {
    `mv $fl $outdir`;
}

print "Running junction_saturation.py\n";
`$PYTHON/junction_saturation.py -i $bam -o "$intdir/$file_pre" -r $bed`;
`test -e "$intdir/$file_pre.junctionSaturation_plot.pdf" && mv "$intdir/$file_pre.junctionSaturation_plot.pdf" $outdir`;

print "Running infer_experiment.py\n";
`$PYTHON/infer_experiment.py -i $bam ">$intdir/$file_pre" -r $bed`;

print "Running inner_distance.py\n";
`$PYTHON/inner_distance.py -i $bam -o "$intdir/$file_pre" -r $bed`;
`test -e "$intdir/$file_pre.inner_distance_plot.pdf" && mv "$intdir/$file_pre.inner_distance_plot.pdf" $outdir`;

`sleep 10`;
exit(0);

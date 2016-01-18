#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);


my ($config, $bin, $species, $counts, $samplekey, $comparisons, $clusterOnly, $diff_out, $count_out, $cluster_out, $gsa_out, $help, $no_replicates);

my $pre = 'TEMP';
GetOptions ('pre=s' => \$pre,
	    'config=s' => \$config,
	    'bin=s' => \$bin,
	    'species=s' => \$species,
	    'counts=s' => \$counts,
	    'samplekey=s' => \$samplekey,
	    'comparisons=s' => \$comparisons,
	    'diff_out=s' => \$diff_out,
	    'count_out=s' => \$count_out,
	    'cluster_out=s' => \$cluster_out,
	    'gsa_out=s' => \$gsa_out,
 	    'clusterOnly' => \$clusterOnly,
            'no_replicates' => \$no_replicates,
	    'help' => \$help) or exit(1);

my $R = '';
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /^r$/i){
	if(!-e "$conf[1]/R"){
	    die "CAN'T FIND R IN $conf[1] $!";
	}
	$R = $conf[1];
    }
}
close CONFIG;

if($clusterOnly){
    `$R/Rscript $bin/RunDE.R "bin='$bin'" "counts.file='$counts'" "clustering.dir='$cluster_out'"`;
}
else{
    open(COMP, "$comparisons") || die "Can't open comparisons file $comparisons $!";
    my @comps = ();
    while(<COMP>){
	chomp;
	
	my @com = split(/\s+/, $_);
	push @comps, "'$com[0] - $com[1]'";
    }
    close COMP;
    
    my $cmpStr = join(",", @comps);
    my $reps = "";
    if($no_replicates){
        $reps = "no.replicates=TRUE";
    }

    my $run_gsa = "gsa.dir='$gsa_out'";
    if(!$gsa_out){
	$run_gsa = "GSA=FALSE";
    }
    
    print "command: $R/Rscript $bin/RunDE.R \"bin='$bin'\" \"species='$species'\" \"proj.id='$pre'\" \"diff.exp.dir='$diff_out'\" \"counts.file='$counts'\" \"counts.dir='$count_out'\" \"clustering.dir='$cluster_out'\" \"$run_gsa\" \"key.file='$samplekey'\" \"comps=c($cmpStr)\" \"$reps\"\n";
    `$R/Rscript $bin/RunDE.R "bin='$bin'" "species='$species'" "proj.id='$pre'" "diff.exp.dir='$diff_out'" "counts.file='$counts'" "counts.dir='$count_out'" "clustering.dir='$cluster_out'" \"$run_gsa\" "key.file='$samplekey'" "comps=c($cmpStr)" "$reps"`;
}

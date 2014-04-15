#!/usr/bin/perl

use strict;

### input is sam file from star
### removes 1.) alignments with isize > (2^31)-1 = 2,147,483,647 
### because this causes error in picard size isize is an int
###         2) alignments where read length does not match qual length; seems to be a problem in star
###         3) where there aren't enough fields (have that as 11 right now; through quals)

my $file = shift;
open(IN, "$file");
open(OUT, ">$file\_filtered.sam");
while(my $line = <IN>){
    chomp $line;


    if($line =~ /^\@/){
	if($line !~ /\@HD/){
	    print OUT "$line\n";
	}
	next;
    }

    my @data = split(/\t/, $line);
    if((abs($data[8]) >= 2**31 - 1) || (length($data[9]) != length($data[10])) || (scalar(@data) < 11)){
	next;
    }

    print OUT "$line\n";

}
close IN;
close OUT;  

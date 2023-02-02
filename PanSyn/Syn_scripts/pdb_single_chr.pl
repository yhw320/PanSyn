#!/usr/bin/perl -w
#Author: Hongwei Yu
#Date: 2022-05-31
#Modified:
#Description:
#use strict;

use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","c=s","o=s","h|help");            
if ( !(defined $opts{i}  and defined $opts{c} and defined $opts{o}) ) {
    die "********************************************************
		-i	Full path to the *.pdb file
		-c	The letter corresponding to a specified single chromosome, where Chr1=a, Chr2=b etc. (example:a)
		-o	Full path to output file, the suffix should be pdb (example: chra.pdb)
*********************************************************\n";
}
if (defined $opts{h} or defined $opts{help}  ) {
    die "********************************************************
		-i	Full path to the *.pdb file
		-c	The letter corresponding to a specified single chromosome, where Chr1=a, Chr2=b etc. (example:a)
		-o	Full path to output file, the suffix should be pdb (example: chra.pdb)
*********************************************************\n";
}

open I,"<$opts{i}";
open O,">$opts{o}";
while(my $a=<I>){
	chomp $a;
	if ($a!~/^TITLE/ and $a!~/^REMARK/) {
		if ($a=~/^MODEL(\s+)1/) {
			print O "$a\n";
		}
		else{
			if ($a=~/chr(\s)$opts{c}/) {
				print O "$a\n";
			}
		}
		if ($a=~/^MODEL(\s+)2/) {
			exit;
		}
	}
}
close I;
close O;



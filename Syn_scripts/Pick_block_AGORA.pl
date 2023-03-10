#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","g=s","b=s","o=s");
if(!(defined $opts{i} and defined $opts{g} and defined $opts{b} and defined $opts{o})){

	die"**********************************************\n
	options:
		-i	Full path to the [ancGenome.*.list] file (generated by AGORA software)
		-g	Full path to the reference species' GFF file
		-b	Full path to the [block.id] file
	    -o	Full path to the output file
		*********************************************\n";
}
my %h_gff=();
open I,"<$opts{g}";
while ($a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h_gff{$it[1]}=$it[0]."\t".$it[2]."\t".$it[3];
}
close I;

my %h_id=();
open I,"<$opts{b}";
while ($a=<I>){
	chomp $a;
	$h_id{$a}=0;
}
close I;

open I,"<$opts{i}";
open O,">$opts{o}";
while ($a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	if (exists $h_id{$it[0]}) {
		my @it2=split/ /,$it[4];
		foreach my $i (@it2) {
			if (exists $h_gff{$i}) {
				my @it3=split/\t/,$h_gff{$i};
				print O "$it[0]\t$it3[0]\t$it3[1]\t$it3[2]\n";
			}
		}
	}
}
close I;close O;
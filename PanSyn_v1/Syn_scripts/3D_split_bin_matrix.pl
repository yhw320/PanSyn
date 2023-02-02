#!/usr/bin/perl
#Author: Hongwei Yu
#Modified:
#Description:
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","b=s","a=s","o=s","h|help");            
if ( !(defined $opts{i} and defined $opts{b} and defined $opts{a} and defined $opts{o}) ) {
    die "********************************************************
	-i	Full path to the matrix file from HicPro (example:500k_iced.matrix)
	-b	Full path to the bed file from HicPro (example:500k_abs.bed)
	-a	Chromosome name of interest (example:chr1) or all chromosomes (input:genome)
	-o	Full path to the output file
*********************************************************\n";
}
if ($opts{a} eq "genome") {
	my %h;my %h1;my %h2;
	open I,"<$opts{i}";
	while(my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		$h{$it[0]}{$it[1]}=$it[2];
		$h{$it[1]}{$it[0]}=$it[2];
		$h1{$it[0]}=0;
		$h2{$it[1]}=0;
	}
	close I;
	my $num=0;
	foreach my $i (keys %h2) {
		$num=$num+1;
	}
	$zero=0;
	open O,">$opts{o}/$opts{a}_matrix.txt";
	foreach my $i1  (sort {$a<=>$b} keys %h1) {
		my $nn=0;
		foreach my $i2 (sort {$a<=>$b} keys %h2) {
			$nn=$nn+1;
			if ($nn<$num) {
				if (exists $h{$i1} and exists $h{$i1}{$i2}) {
					print O "$h{$i1}{$i2}\t";
				}
				else{
					print O "$zero\t";
				}
			}
			else{
				if (exists $h{$i1} and exists $h{$i1}{$i2}) {
					print O "$h{$i1}{$i2}\n";
				}
				else{
					print O "$zero\n";
				}
			}
		}
	}
	close O;
}
else{
	open I,"<$opts{b}";
	while(my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		if ($it[0] eq $opts{a}) {
			$h_chr_num{$it[3]}=0;
		}
	}
	close I;
	
	my %h;my %h1;my %h2;
	open I,"<$opts{i}";
	while(my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		if (exists $h_chr_num{$it[0]} and exists $h_chr_num{$it[1]}) {
			$h{$it[0]}{$it[1]}=$it[2];
			$h{$it[1]}{$it[0]}=$it[2];
			$h1{$it[0]}=0;
			$h2{$it[1]}=0;
		}
	}
	close I;
	my $num=0;
	foreach my $i (keys %h2) {
		$num=$num+1;
	}
	$zero=0;
	open O,">$opts{o}/$opts{a}_matrix.txt";
	foreach my $i1  (sort {$a<=>$b} keys %h1) {
		my $nn=0;
		foreach my $i2 (sort {$a<=>$b} keys %h2) {
			$nn=$nn+1;
			if ($nn<$num) {
				if (exists $h{$i1} and exists $h{$i1}{$i2}) {
					print O "$h{$i1}{$i2}\t";
				}
				else{
					print O "$zero\t";
				}
			}
			else{
				if (exists $h{$i1} and exists $h{$i1}{$i2}) {
					print O "$h{$i1}{$i2}\n";
				}
				else{
					print O "$zero\n";
				}
			}
		}
	}
	close O;
}


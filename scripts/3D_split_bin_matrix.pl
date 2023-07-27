#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","b=s","a=s","o=s","h|help");            
if ( !(defined $opts{i} and defined $opts{b} and defined $opts{a} and defined $opts{o}) ) {
    die "********************************************************
	-i	Full path to the [*.matrix] file (example:500k_iced.matrix)
	-b	Full path to the [*.bed] file (example:500k_abs.bed)
	-a	Enter the name of the chromosome of interest (e.g. chr1) or type [genome] to include all chromosomes.
	-o	Full path to the [outputDir] directory cotaining output files
	optional:
	-h|-help Print this help page
*********************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [*.matrix] file (example:500k_iced.matrix)
	-b	Full path to the [*.bed] file (example:500k_abs.bed)
	-a	Enter the name of the chromosome of interest (e.g. chr1) or type [genome] to include all chromosomes.
	-o	Full path to the [outputDir] directory cotaining output files
	optional:
	-h|-help Print this help page
		*************************************************\n";
}
my $slash;
####
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####

my $zero; my %h_chr_num=(); 

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


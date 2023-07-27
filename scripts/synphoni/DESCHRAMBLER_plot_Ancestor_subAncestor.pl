#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"ic=s","im=s","c=s","o8=s","h|help");
if(!(defined $opts{ic} and defined $opts{im} and defined $opts{c} and defined $opts{o8})){

	die"**********************************************\n
	options:
		-ic	Full path to the [chrom_sizes] directory
		-im Full path to the [maps] directory
		-c	Full path to the [Order_of_evolution.txt] file
		-o8	Full path to the new directory containing output files
		Optional:
		-h|-help Print this help page
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
		-ic	Full path to the [chrom_sizes] directory
		-im Full path to the [maps] directory
		-c	Full path to the [Order_of_evolution.txt] file
		-o8	Full path to the new directory containing output files
		Optional:
		-h|-help Print this help page
		*********************************************\n";
}

my $slash;
###
if ($opts{ic} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{ic} =~ s/$slash$//;
}
####
my $slash;
###
if ($opts{im} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{im} =~ s/$slash$//;
}
####
my $slash;
###
if ($opts{o8} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o8} =~ s/$slash$//;
}
####

$dir1=$opts{ic};
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^(\S+)-APCF_size.txt$/ ) {
		$file=$dir1."/$aa";

		if (-d "$opts{o8}/$1") {
			print "The $opts{o8}/$1 directory already exists!\nThe old directory has been deleted.\n";
			system "rm -r $opts{o8}/$1";
		}

		system"mkdir $opts{o8}/$1";
		open O,">$opts{o8}/$1/$1.chromosome" or die "Could not open outfile.\n";
        open IN,"<","$file" or die;
        while ($a=<IN>){
			chomp $a;
			my @items=split/\t/,$a;
			if ($items[0] ne "total") {
				print O "$items[0]\t0\t$items[1]\n";
			}
		}
	}
	close O;close IN;
}
close I;
##################################################################
$dir1=$opts{im};
opendir I,$dir1;
my $houzhui_name;my $hg19;
while ($aa=readdir I) {
	if ( $aa=~/APCF_(\S+)\.merged\.map$/ ) {
		$houzhui_name="APCF_".$1."\.merged\.map";
		$hg19=$1;
	}
}
close I;
##################################################################
my $n=0;

if (-d "$opts{o8}/reference_bed") {
	print "The $opts{o8}/reference_bed directory already exists!\nThe old directory has been deleted.\n";
	system "rm -r $opts{o8}/reference_bed";
}


system "mkdir $opts{o8}/reference_bed";
open II,"<","$opts{c}" or die;
while ($a=<II>){
	chomp $a;
	$n=$n+1;
	if ($n==1) {
		open IN,"<","$opts{im}/$a-$houzhui_name" or die;
		open O,">","$opts{o8}/reference_bed/$a\_reference.bed" or die;
		while ($b=<IN>){
			chomp $b;
			if ($b=~/^APCF\.(\S+):/) {
				$h{'ha'}=$a."-".$1;
			}
			if ($b=~$hg19) {
				my @it=split/ /,$b;
				if ($it[0]=~/\.(\S+):(\S+)-(\S+)$/) {
					my $na1=$1;my $na2=$2;my $na3=$3;
					print O "$na1\t$na2\t$na3\t$h{'ha'}\n";
				}
			}
		}
		close O;
		open IN2,"<","$opts{o8}/$a/$a.chromosome" or die;
		open O2,">","$opts{o8}/$a/$a\_plot.bed" or die;
		while ($c=<IN2>){
			chomp $c;
			my @it=split/\t/,$c;
			$color=$a."-".$it[0];
			$jishu{$color}=0;
			print O2 "$it[0]\t$c\t$color\n";
		}
		close IN;close O2;close IN2;
	}
	else{
		open IN,"<","$opts{im}/$a-$houzhui_name" or die;
		open O,">","$opts{o8}/reference_bed/$a\_reference.bed" or die;
		while ($b=<IN>){
			chomp $b;
			if ($b=~/^APCF\.(\S+):(\S+)-(\S+)/) {
				$start=$2;$end=$3;
				$h{'ha'}=$a."-".$1."_".$start."-".$end;
			}
			if ($b=~$hg19) {
				my @it=split/ /,$b;
				if ($it[0]=~/\.(\S+):(\S+)-(\S+)$/) {
					my $na1=$1;my $na2=$2;my $na3=$3;
					print O "$na1\t$na2\t$na3\t$h{'ha'}\n";
				}
			}
		}
		close IN;close O;
	}
}
close II;


###########second####################
$n=0;
open II,"<","$opts{c}" or die;
while ($a=<II>){
	chomp $a;
	$n=$n+1;
	if ($n==1) {
		$hsap_file=$a."_reference.bed";
		$hsap_name=$a;
	}
	else {
		$other_file=$a."_reference.bed";
		system "bedtools intersect -a $opts{o8}/reference_bed/$hsap_file -b $opts{o8}/reference_bed/$other_file -wo > $opts{o8}/$hsap_name\_$a\_common.bed";
		system "sort -k 1,1 -k 2n,2 $opts{o8}/$hsap_name\_$a\_common.bed >$opts{o8}/$hsap_name\_$a\_common.bed.sorted";
		open IN,"<","$opts{o8}/$hsap_name\_$a\_common.bed.sorted" or die;
		open O,">","$opts{o8}/1.bed" or die;
		while ($b=<IN>){
			chomp $b;
			my @it1=split/\t/,$b;
			$cha=$it1[2]-$it1[1];
			print O "$it1[3]\t$cha\t$it1[7]\n";
		}
		close IN;close O;
		my %h1=();my %h_sum=();
		open IN,"<","$opts{o8}/1.bed" or die;
		open O,">","$opts{o8}/2.bed" or die;
		while ($b=<IN>){
			chomp $b;
			my @it1=split/\t/,$b;
			if (!exists $h_sum{$it1[2]}) {
				$h_sum{$it1[2]}=$it1[1];
			}
			else{
				$h_sum{$it1[2]}=$h_sum{$it1[2]}+$it1[1];
			}
			if (!exists $h1{$it1[2]}) {
				$n1=1;
				$h1{$it1[2]}{$n1}=$it1[0]."\t".$it1[1];
			}
			else{
				$n1=$n1+1;
				$h1{$it1[2]}{$n1}=$it1[0]."\t".$it1[1];
			}
		}
		foreach my $i1 (keys %h1) {
			foreach my $i2 (sort keys %{$h1{$i1}}) {
				my @it3=split/\t/,$h1{$i1}{$i2};
				$bizhi=$it3[1]/$h_sum{$i1};
				print O "$it3[0]\t$i1\t$bizhi\n";
			}
		}
		close IN;close O;
		my %h_start=();
		open IN,"<","$opts{o8}/2.bed" or die;
		open O,">","$opts{o8}/$a/$a\_plot.bed" or die;
		while ($b=<IN>){
			chomp $b;
			my @it=split/\t/,$b;
			my @it2=split/_/,$it[1];
			my @it3=split/-/,$it2[1];
			my @it4=split/-/,$it2[0];
			if (!exists $h_start{$it[1]}) {
				$start=$it3[0];
				if ($it[2]==1) {
					$end=$it3[1];
				}
				else{
					$he=$it3[1]-$it3[0];
					$end=$start+$he*$it[2];
				}
				print O "$it4[1]\t$it4[1]\t$start\t$end\t$it[0]\n";
				$h_start{$it[1]}=$end;
			}
			else{
				$start=$h_start{$it[1]};
				$he=$it3[1]-$it3[0];
				$end=$start+$he*$it[2];
				print O "$it4[1]\t$it4[1]\t$start\t$end\t$it[0]\n";
				$h_start{$it[1]}=$end;
			}
		}
		close IN;close O;
		system "rm $opts{o8}/$hsap_name\_$a\_common.bed";
	}
}
close II;
	
#system "rm $opts{o8}/1.bed";
#system "rm $opts{o8}/2.bed";

open O,">$opts{o8}/color.txt" or die "Could not open outfile.\n";
foreach my $i (keys %jishu) {
	print O "$i\n";
}
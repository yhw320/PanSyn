#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","list=s","o=s","h|help");
if(!(defined $opts{i} and defined $opts{list} and defined $opts{o})){

	die"**********************************************\n
	-i	Full path to the [SynNetBuild*] directory
	-list	Full path to the [Species_list.txt] file
	-o	Full path to the [outputDir_S3] directory
	Optional:
	-h|-help Print this help page
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
	die"**********************************************\n
	-i	Full path to the [SynNetBuild*] directory
	-list	Full path to the [Species_list.txt] file
	-o	Full path to the [outputDir_S3] directory
	Optional:
	-h|-help Print this help page
		*********************************************\n";
}
my $slash;
####
if ($opts{i} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{i} =~ s/$slash$//;
}
####
####
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####

my %h_gene_num=(); my $sum_name;my $sum_name1;my $sum_name2;my $ji;my $sum;my $n;my $all; my %h_final=();my %h_or2=();my %h_order=();my $nn;my $n1;my $n2;my $zhi;my %h_dui=();my $i;my $w;my $w2;
my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).gff$/ and $c!~/_/) {
		my $spe_name=$1;
		my $file=$dir1."/$c";
		my $n=0;
		open B1,"<$file";
		while (my $a=<B1>){
			chomp $a;
			$n=$n+1;
		}
		close B1;
		$h_gene_num{$spe_name}=$n;
	}
}
close P;
############################################
my %h_tandem=();
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).tandem.collinear$/) {
		my $spe_name=$1;
		my $file=$dir1."/$c";
		open B1,"<$file";
		while (my $a=<B1>){
			if ($a!~/^Anchor1/) {
				chomp $a;
				my @it=split/[,|\t]+/,$a;
				foreach my $i (@it) {
					$h_tandem{$spe_name}{$i}=0;
				}
			}
		}
		close B1;
	}
}
close P;   
############################################
my %h_pair=();
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).collinearity$/ and $c!~/_/) {
		my $spe_name=$1;
		my $file=$dir1."/$c";
		open B1,"<$file";
		while (my $a=<B1>){
			if ($a!~/^#/) {
				chomp $a;
				my @it=split/\t/,$a;
				$sum_name=$spe_name."\t".$spe_name;
				$h_pair{$sum_name}{$it[1]}=0;
				$h_pair{$sum_name}{$it[2]}=0;
			}
		}
		close B1;
	}
	if ( $c=~/^(\S+)_(\S+).collinearity$/) {
		my $spe_name1=$1;
		my $spe_name2=$2;
		my $file=$dir1."/$c";
		open B1,"<$file";
		while (my $a=<B1>){
			if ($a!~/^#/) {
				chomp $a;
				my @it=split/\t/,$a;
				$sum_name1=$spe_name1."\t".$spe_name2;
				$sum_name2=$spe_name2."\t".$spe_name1;
				$ji=$spe_name1."_";
				if ($it[1]=~/$ji/) {
					$h_pair{$sum_name1}{$it[1]}=0;
					$h_pair{$sum_name2}{$it[2]}=0;
				}
				else{
					$h_pair{$sum_name1}{$it[2]}=0;
					$h_pair{$sum_name2}{$it[1]}=0;
				}
			}
		}
		close B1;
	}
}
close P;

#######################################
foreach my $i1 (keys %h_pair) {
	$n=0;$sum=0;
	my @it=split/\t/,$i1;
	my %h=();
	foreach my $i2 (keys %{$h_pair{$i1}}) {
		$h{$i2}=0;
	}
	foreach my $i3 (keys %h_tandem) {
		if ($i3 eq $it[0]) {
			foreach my $i4 (keys %{$h_tandem{$i3}}) {
				$h{$i4}=0;
			}
		}
	}
	foreach my $i5 (keys %h) {
		$sum=$sum+1;
	}
	$all=$h_gene_num{$it[0]};
	$zhi=$sum/$all;
	$h_final{$i1}=$zhi;
}
################################33
open O,">$opts{o}/Syntenic_percentages.txt";
foreach my $i (keys %h_final) {
	print O "$i\t$h_final{$i}\n"; 
}
close O;


##################################
open I,"<$opts{list}";
$nn=1;
while (my $a=<I>){
	chomp $a;
	$h_order{$a}=$nn;
	$h_or2{$nn}=$a;
	$nn=$nn+1;
}
close I;
$n=0;
foreach my $i (keys %h_order) {
	$n=$n+1;
}
##################################3
open I,"<$opts{o}/Syntenic_percentages.txt";
while (my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$n1=$h_order{$it[0]};
	$n2=$h_order{$it[1]};
	$h_dui{$n1}{$n2}=$it[2];
}
close I;

open O,">$opts{o}/Syntenic_percentages.matrix";
print O "Species\t";
for ($i=1;$i<$n;$i=$i+1) {
	print O"$h_or2{$i}\t";
}
print O "$h_or2{$n}\n";

for ($w=1;$w<=$n;$w=$w+1) {
	print O "$h_or2{$w}\t";
	for ($w2=1;$w2<$n;$w2=$w2+1) {
		if (exists $h_dui{$w} and exists $h_dui{$w}{$w2}) {
			print O "$h_dui{$w}{$w2}\t";
		}
		else{
			print O "0\t";
		}
	}
	if (exists $h_dui{$w} and exists $h_dui{$w}{$n}) {
		print O "$h_dui{$w}{$n}\n";
	}
	else{
		print O "0\n";
	}
}
close O;


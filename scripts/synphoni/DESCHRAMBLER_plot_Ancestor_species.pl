#!/usr/bin/perl
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"ic=s","im=s","o7=s","h|help");
if(!(defined $opts{ic} and defined $opts{im} and defined $opts{o7})){

	die"**********************************************\n
	options:
		-ic	Full path to the [chrom_sizes] directory
		-im Full path to the [maps] directory
		-o7	Full path to the new directory containing output files
		Optional:
		-h|-help Print this help page
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
		-ic	Full path to the [chrom_sizes] directory
		-im Full path to the [maps] directory
		-o7	Full path to the new directory containing output files
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
if ($opts{o7} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o7} =~ s/$slash$//;
}
####

$dir1=$opts{ic};
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^(\S+)_chrom.sizes$/ ) {
		$file=$dir1."/$aa";

		if (-d "$opts{o7}/$1") {
			print "The $opts{o7}/$1 directory already exists!\nThe old directory has been deleted.\n";
			system "rm -r $opts{o7}/$1";
		}

		system"mkdir $opts{o7}/$1";
		open O,">$opts{o7}/$1/$1.chromosome" or die "Could not open outfile.\n";
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
%jishu=();
$dir1=$opts{im};
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^APCF_(\S+).merged.map$/ ) {
		$file=$dir1."/$aa";
		$name=$1;
		open O,">$opts{o7}/$1/$1\_annotation_pos.txt" or die "Could not open outfile.\n";
        open IN,"<","$file" or die;
        while ($a=<IN>){
			chomp $a;
			if ($a=~/^APCF\.(\S+):/) {
				$h{'ha'}=$1;
				$jishu{$1}=0;
			}
			if ($a=~$name) {
				my @it=split/ /,$a;
				if ($it[0]=~/\.(\S+):(\S+)-(\S+)$/) {
					my $na1=$1;my $na2=$2;my $na3=$3;
					print O "$na1\t$na1\t$na2\t$na3\t$h{'ha'}\n";
				}
			}
		}
	}
	close O;close IN;
}

close I;

open O,">$opts{o7}/color.txt" or die "Could not open outfile.\n";
foreach my $i (keys %jishu) {
	print O "$i\n";
}
#!/usr/bin/perl
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i1=s","i2=s","o=s");
if(!(defined $opts{i1} and defined $opts{i2} and defined $opts{o})){

	die"**********************************************\n
	options:
		-i1	Full path to the chrom_sizes fold
		-i2 Full path to the maps fold
		-o	Full path to the folder containing output files
		*********************************************\n";
}
$dir1=$opts{i1};
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^(\S+)_chrom.sizes$/ ) {
		$file=$dir1."/$aa";
		system"mkdir $opts{o}/$1";
		open O,">$opts{o}/$1/$1.chromosome" or die "Could not open outfile.\n";
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
$dir1=$opts{i2};
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^APCF_(\S+).merged.map$/ ) {
		$file=$dir1."/$aa";
		$name=$1;
		open O,">$opts{o}/$1/$1\_annotation_pos.txt" or die "Could not open outfile.\n";
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

open O,">$opts{o}/color.txt" or die "Could not open outfile.\n";
foreach my $i (keys %jishu) {
	print O "$i\n";
}
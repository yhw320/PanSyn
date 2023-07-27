#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i2=s","t=s","color=s","o6=s","width=s","height=s","h|help");
if (!( defined $opts{i2} and defined $opts{t} and defined $opts{color} and defined $opts{o6})) {
		die "************************************************\n
	-i2	Full path to the new [inputDir2] directory containing input files
	-t	Full path to the [tree.nwk] file
	-color	Full path to the Ancestor'chromosome color file
	-o6	Full path to the new [outputDir6] directory storing the output files
	-Optional:
	-width	Specify the width of the output image (default: 22)
	-height	Specify the height of the output image (default: 4)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i2	Full path to the new [inputDir2] directory containing input files
	-t	Full path to the [tree.nwk] file
	-color	Full path to the Ancestor'chromosome color file
	-o6	Full path to the new [outputDir6] directory storing the output files
	-Optional:
	-width	Specify the width of the output image (default: 22)
	-height	Specify the height of the output image (default: 4)
	-h|-help Print this help page
		*************************************************\n";
}
my $slash;
###
if ($opts{i2} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{i2} =~ s/$slash$//;
}
####
###
if ($opts{o6} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o6} =~ s/$slash$//;
}
####

if (!(defined $opts{width})) {
	$opts{width}=22;
}
if (!(defined $opts{height})) {
	$opts{height}=4;
}

my $dir1=$opts{i2};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_chr_breakage_fusion.result$/ ) {
		my $file=$dir1."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,"> $opts{o6}/$c";
	    while (my $a=<I>) {
		    chomp $a;
			if ($a=~/Ancestor/ and $a=~/Chr/) {
		    }
		    else {
				if ($a=~/_part/) {
					$a=~s/_part//g;
				}
				if ($a=~/_all/) {
					$a=~s/_all//g;
				}
				print O "$a\n";
			}
		}
    close I;close O;
	}
}
close P;

system "python3 infer_chr_fusion.py --tree_dir \"$opts{t}\" --work_space \"$opts{o6}\" ";
##########
open I,"< $opts{color}";
open O,"> $opts{o6}/color.txt";
while (my $a=<I>) {
	chomp $a;
	my @it=split/\t/,$a;
	print O "$it[0]_block\t$it[1]\n";
}
close I;close O;

#########

$dir2="$opts{o6}/ancestor";
opendir P,$dir2;
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_anc.result$/ ) {
		$out_name=$1;
		my $file=$dir2."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,"> $opts{o6}/$c-R";
		my %jishu=();
	    while (my $a=<I>) {
		    chomp $a;
			my @it=split/\t/,$a;
			my @it2=split/,/,$it[1];
			foreach my $i2 (@it2) {
				if (!exists $jishu{$it[0]}) {
					$jishu{$it[0]}=1;
				}
				else{
					$jishu{$it[0]}=$jishu{$it[0]}+1;
				}
			}
		}
		close I;
		open I,"<$file" or die("Could not open $file\n"); 
		print  O "Chr\tNumber\tAncestor\n";
	    while (my $a=<I>) {
		    chomp $a;
			my @it=split/\t/,$a;
			my @it2=split/,/,$it[1];
			foreach my $i3 (@it2) {
				$zhi=1/$jishu{$it[0]};
				print  O "$it[0]\t$zhi\t$i3\_block\n";
			}
		}
		close I;
		close O;
		system "Rscript Chr_breakage_fusion2.R $opts{o6} $opts{o6}/color.txt $opts{o6}/$c-R $out_name $opts{width} $opts{height}";
	}
}
close P;




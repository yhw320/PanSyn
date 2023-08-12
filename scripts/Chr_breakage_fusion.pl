#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","t=s","color=s","o=s","width=s","height=s","h|help");
if (!( defined $opts{i} and defined $opts{t} and defined $opts{color} and defined $opts{o})) {
		die "************************************************\n
	-i	Full path to the [inputDir_S23] directory
	-t	Full path to the [tree.tre] file
	-color	Full path to the [*_ancestor_chr_color.txt] file
	-o	Full path to the [outputDir_S23] directory
	-Optional:
	-width	The width of the output image (default: 22)
	-height	The height of the output image (default: 4)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [inputDir_S23] directory
	-t	Full path to the [tree.tre] file
	-color	Full path to the [*_ancestor_chr_color.txt] file
	-o	Full path to the [outputDir_S23] directory
	-Optional:
	-width	The width of the output image (default: 22)
	-height	The height of the output image (default: 4)
	-h|-help Print this help page
		*************************************************\n";
}
my $slash;
###
if ($opts{i} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{i} =~ s/$slash$//;
}
####
###
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####

if (!(defined $opts{width})) {
	$opts{width}=22;
}
if (!(defined $opts{height})) {
	$opts{height}=4;
}

my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_chr_breakage_fusion.result$/ ) {
		my $file=$dir1."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,"> $opts{o}/$c";
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

system "infer_chr_fusion --tree_dir \"$opts{t}\" --work_space \"$opts{o}\" ";
##########
open I,"< $opts{color}";
open O,"> $opts{o}/color.txt";
while (my $a=<I>) {
	chomp $a;
	my @it=split/\t/,$a;
	print O "$it[0]_block\t$it[1]\n";
}
close I;close O;

#########

$dir2="$opts{o}/ancestor";
opendir P,$dir2;
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_anc.result$/ ) {
		$out_name=$1;
		my $file=$dir2."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,"> $opts{o}/$c-R";
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
		system "Chr_breakage_fusion2 $opts{o} $opts{o}/color.txt $opts{o}/$c-R $out_name $opts{width} $opts{height}";
	}
}
close P;




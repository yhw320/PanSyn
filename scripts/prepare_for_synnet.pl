#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"g=s","p=s","o=s","h|help");
if (!( defined $opts{g} and defined $opts{p} and defined $opts{o})) {
		die "************************************************\n
	-g	Full path to the [gffs_S3] directory
	-p	Full path to the [peps_S3] directory
	-o	Full path to the [outputDir_S3] directory
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-g	Full path to the [gffs_S3] directory
	-p	Full path to the [peps_S3] directory
	-o	Full path to the [outputDir_S3] directory
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}

my $slash;
####
if ($opts{g} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{g} =~ s/$slash$//;
}
####
####
if ($opts{p} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{p} =~ s/$slash$//;
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


if (-d "$opts{o}/prepare_data") {
	print "The $opts{o}/prepare_data directory already exists!\nThe old directory has been deleted.\n";
	system "rm -r $opts{o}/prepare_data";
}

system "mkdir $opts{o}/prepare_data";

my $dir1;my $name;my $name2;

$dir1=$opts{p};
opendir P,$dir1;
while (my $c=readdir P) {
	if ($c=~/^(\S+).pep$/) {
		$name=$1;
		my $file=$dir1."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,">$opts{o}/prepare_data/$name.pep"; 
		while (my $a=<I>) {
			chomp $a;
			if ($a=~/>(\S+)/) {
				$name2=$1;
				print O ">$name\_$name2\n";
			}
			else{print O "$a\n";}
		}
		close I;close O;
	}
}
close P;

$dir1=$opts{g};
opendir P,$dir1;
while (my $c=readdir P) {
	if ($c=~/^(\S+)_simplified.gff$/) {
		$name=$1;
		my $file=$dir1."/$c";
		open I,"<$file" or die("Could not open $file\n"); 
		open O,">$opts{o}/prepare_data/$name.bed"; 
		while (my $a=<I>) {
			chomp $a;
			my @it=split/\t/,$a;
			print O "$name\_$it[0]\t$name\_$it[1]\t$it[2]\t$it[3]\n";
		}
		close I;close O;
	}
}
close P;




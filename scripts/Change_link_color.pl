#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","c=s","o=s","h|help");
if(!(defined $opts{i} and defined $opts{c} and defined $opts{o})){
		die "************************************************\n
	-i	Full path to the [All_species_microsyn_genes.links] file
	-c	Full path to the [chr_link_color.txt] file
	-o	Full path to the [outputDir_S1] directory
	Optional:
	-h|-help Print this help page
	*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [All_species_microsyn_genes.links] file
	-c	Full path to the [chr_link_color.txt] file
	-o	Full path to the [outputDir_S1] directory
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}

####
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    my $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####


open I,"<$opts{c}";
while (my $a=<I>){
	chomp $a;
	my @items=split/\t/,$a;
	$h_color{$items[0]}=$items[1];
}
close I;

open I,"<$opts{i}";
open O,">$opts{o}/genes.links";
while (my $a=<I>){
	chomp $a;
	my @items=split/\t/,$a;
	if (exists $h_color{$items[0]}) {
		print O "$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]\t$items[5]\tcolor=$h_color{$items[0]}\n";
	}
}
close I;
close O;

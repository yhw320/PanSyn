#!/usr/bin/perl
use Getopt::Long;
use warnings;
my %opts;

GetOptions(\%opts,"i=s");

open I,"<$opts{i}";
while (my $line =<I>){
	chomp $line;
	next if $line =~ /^\s*$/;  # 忽略空行
	system($line);
}
close I;

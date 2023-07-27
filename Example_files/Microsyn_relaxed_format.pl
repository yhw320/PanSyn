#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"p=s","g=s","o=s","h|help");
if (!( defined $opts{p} and defined $opts{g} and defined $opts{o})) {
		die "************************************************\n
	-p	Full path to the [gffs] directory containing [*_simplified.gff] files
	-g	Full path to the [peps] directory containing [*.pep] files
	-o	Full path to the [outputDir] directory containing output files
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-p	Full path to the [gffs] directory containing [*_simplified.gff] files
	-g	Full path to the [peps] directory containing [*.pep] files
	-o	Full path to the [outputDir] directory containing output files
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}

####################################

my $slash;
####
if ($opts{p} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{p} =~ s/$slash$//;
}
####

####
if ($opts{g} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{g} =~ s/$slash$//;
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

if (-d "$opts{o}/peps") {
	system "rm -r $opts{o}/peps";
} 
system "mkdir $opts{o}/peps";
opendir P,$opts{p};
while (my $c=readdir P) {
	if ( $c=~/^(\S+).pep$/ ) {
		my $file=$opts{p}."/$c";
		system "sed 's/^>/>$1\_/' $file >$opts{o}/peps/$1.pep";
	}
}
close P;

if (-d "$opts{o}/gffs") {
	system "rm -r $opts{o}/gffs";
} 
system "mkdir $opts{o}/gffs";
opendir P,$opts{g};
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_simplified.gff$/ ) {
		my $file=$opts{g}."/$c";
		my $name=$1;
		open IN,"<$file" or die("Could not open $file\n"); 
		open O,">$opts{o}/gffs/$1.chrom" or die("Could not open $opts{o}/gffs/$1.chrom\n"); 
		while (my $a=<IN>) {
			chomp $a;
			my @it=split/\t/,$a;
			print O "$name\t$name\_$it[1]\t$it[0]\t$it[4]\t$it[2]\t$it[3]\n";
		}
		close IN;close O;
	}
}
close P;

#system "orthofinder -f $opts{o}/peps -t $opts{t}";

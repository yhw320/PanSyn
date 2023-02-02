#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"p=s","g=s","o=s","h|help");
if (!( defined $opts{p} and defined $opts{g} and defined $opts{o})) {
		die "************************************************\n
	-p	Full path to the [gffs] folder containing gff files
	-g	Full path to the [peps] folder containing pep files
	-o	Full path to the [outputDir] folder containing output files
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-p	Full path to the [gffs] folder containing gff files
	-g	Full path to the [peps] folder containing pep files
	-o	Full path to the [outputDir] folder containing output files
	-h|-help Print this help page
		*************************************************\n";
}

####################################

system "mkdir $opts{o}/peps";
opendir P,$opts{p};
while (my $c=readdir P) {
	if ( $c=~/^(\S+).pep$/ ) {
		my $file=$opts{p}."/$c";
		system "sed 's/^>/>$1\_/' $file >$opts{o}/peps/$1.pep";
	}
}
close P;

system "mkdir $opts{o}/gffs";
opendir P,$opts{g};
while (my $c=readdir P) {
	if ( $c=~/^(\S+).gff$/ ) {
		my $file=$opts{g}."/$c";
		$name=$1;
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

#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"g=s","p=s","o=s","h|help");
if (!( defined $opts{g} and defined $opts{p} and defined $opts{o})) {
		die "************************************************\n
	-g	Full path to the folder containing gff data files
	-p	Full path to the folder containing pep data files
	-o	Full path to the folder storing the output results
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-g	Full path to the folder containing gff data files
	-p	Full path to the folder containing pep data files
	-o	Full path to the folder storing the output results
	-h|-help Print this help page
		*************************************************\n";
}
system "mkdir $opts{o}/prepare_data";


my $dir1=$opts{p};
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

my $dir1=$opts{g};
opendir P,$dir1;
while (my $c=readdir P) {
	if ($c=~/^(\S+).gff$/) {
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




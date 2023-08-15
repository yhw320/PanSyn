#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"p1=s","p2=s","o=s","a=s","t=s","e=s","m=s","h|help");
if (!( defined $opts{p1} and defined $opts{p2} and defined $opts{a} and defined $opts{o})) {
		die "************************************************\n
	-p1	Full path to pep file of the reference species
	-p2	Full path to pep file of the interested species
	-o	Full path to the [outputDir_S9B] directory
	-a	Alignment software (diamond or blast)
	Optional:
	-t	Number of threads (default: 12)
	-e	Sequence alignment E-value (default: 1e-5 for blast, 1e-3 for diamond)
	-m	The number of best non-self alignment hits (default:1)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-p1	Full path to pep file of the reference species
	-p2	Full path to pep file of the interested species
	-o	Full path to the [outputDir_S9B] directory
	-a	Alignment software (diamond or blast)
	Optional:
	-t	Protein sequence alignment threads (default:12)
	-e	Protein sequence alignment evalue (default:1e-5/0.001)
	-m	The number of best non-self alignment hits (default:1)
	-h|-help	Print this help page
		*************************************************\n";
}

my $slash;

####
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{e})) {
	if ($opts{a} eq "blast") {
		$opts{e}=1e-5;
	}
	if ($opts{a} eq "diamond") {
		$opts{e}=0.001;
	}
}
if (!(defined $opts{m})) {
	$opts{m}=1;
}

######################
my $dir = "$opts{o}"; 

opendir my $dh, $dir or die "can't open $dir: $!";

my @files = readdir($dh);       
closedir($dh);   

if (@files == 2) {
} 
else {
	system "rm -r $opts{o}/*";
}

#######################

if ($opts{a} eq "blast") {

	my $dirphr = "$opts{p2}.phr";
	my $dirpin = "$opts{p2}.pin";
	my $dirpsq = "$opts{p2}.psq";
	if (-e $dirphr and -e $dirpin and -e $dirpsq) {
	}
	else{
		system "makeblastdb -in $opts{p2} -dbtype prot";
	}
	system "blastp -query $opts{p1} -db $opts{p2} -outfmt 6 -out $opts{o}/1_2.blastp -evalue $opts{e} -num_threads $opts{t} -max_target_seqs $opts{m}";				
	open I,"< $opts{o}/1_2.blastp";
	open O,"> $opts{o}/1.pairs";
	while (my $a=<I>) {
		chomp $a;
		my @itms=split/\t/,$a;
		print O "$itms[0]\t$itms[1]\n";
	}
	close O;

	system "sort $opts{o}/1.pairs | uniq >$opts{o}/1_2.pairs";
}
if ($opts{a} eq "diamond") {

	my $dirphr = "$opts{p2}.database.dmnd";
	if (-e $dirphr) {
	}
	else{
		system "diamond makedb --in $opts{p2} -d $opts{p2}.database";
	}
	system "diamond blastp --query $opts{p1} --db $opts{p2} --outfmt 6 --out $opts{o}/1_2.blastp  --evalue $opts{e} --threads $opts{t} --max-target-seqs $opts{m}";

	open I,"< $opts{o}/1_2.blastp";
	open O,"> $opts{o}/1.pairs";
	while (my $a=<I>) {
		chomp $a;
		my @itms=split/\t/,$a;
		print O "$itms[0]\t$itms[1]\n";
	}
	close O;

	system "sort $opts{o}/1.pairs | uniq >$opts{o}/1_2.pairs";
}
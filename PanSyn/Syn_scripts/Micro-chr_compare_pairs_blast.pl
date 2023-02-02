#!/usr/bin/perl -w
use Getopt::Long;
my %opts;
GetOptions(\%opts,"p1=s","p2=s","o=s","t=s","e=s","m=s","h|help");
if (!( defined $opts{p1} and defined $opts{p2} and defined $opts{o})) {
		die "************************************************\n
	-p1	Full path to the reference genome' pep file
	-p2	Full path to the query genome' pep file
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	Blast alignment threads (default:12)
	-e	Blast alignment evalue (default:1e-5)
	-m	The number of best non-self BLASTP hits (default:1)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-p1	Full path to the reference genome' pep file
	-p2	Full path to the query genome' pep file
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	Blast alignment threads (default:12)
	-e	Blast alignment evalue (default:1e-5)
	-m	The number of best non-self BLASTP hits (default:1)
	-h|-help	Print this help page
		*************************************************\n";
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{e})) {
	$opts{e}=1e-5;
}
if (!(defined $opts{m})) {
	$opts{m}=1;
}

system "makeblastdb -in $opts{p2} -dbtype prot";

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
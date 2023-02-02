#!/usr/bin/perl-w
#Author:Hongwei Yu
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","m=s","t=s","e=s","o=s","s=s","h|help");
if(!(defined $opts{i} and defined $opts{o})){
	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing required data files
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	Blast alignment threads (default:12)
	-e	Blast alignment evalue (default:1e-5)
	-s	The number of best non-self BLASTP hits (default:5)
	-m	Full path where the MCScanX software resides (/software/MCScanX-master/MCScanX)
		*********************************************\n";
}

if (defined $opts{h} or defined $opts{help}) {
	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing required data files
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	Blast alignment threads (default:12)
	-e	Blast alignment evalue (default:1e-5)
	-s	The number of best non-self BLASTP hits (default:5)
	-m	Full path where the MCScanX software resides (/software/MCScanX-master/MCScanX)
		*********************************************\n";
}

if (!(defined $opts{e})) {
	$opts{e}=1e-5;
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{s})) {
	$opts{s}=5;
}

if (!(defined $opts{m})) {
	$opts{m}="MCScanX";
}
$n=1;
my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).pep$/ ) {
		my $file=$dir1."/$c";
		system "makeblastdb -in $file -dbtype prot";
		$h{$n}=$1;
		$n=$n+1;
	}
}
close P;
$num=keys %h;
system "mkdir $opts{o}/microsyn_genes_links";
system "mkdir $opts{o}/MCScanX_result";
foreach my $i (sort {$a<=>$b} keys %h) {
	for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
		if ($i!=$num) {
			system "mkdir $opts{o}/$h{$i}_$h{$n1}";
			system "blastp -query $opts{i}/$h{$i}.pep -db $opts{i}/$h{$n1}.pep -outfmt 6 -out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast -evalue $opts{e} -num_threads $opts{t} -max_target_seqs $opts{s}";
			system "cat $opts{i}/$h{$i}.gff $opts{i}/$h{$n1}.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
			chdir "$opts{o}/$h{$i}_$h{$n1}";
			system "$opts{m} $h{$i}_$h{$n1}";
			system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/MCScanX_result";
			posi2($h{$i},$h{$n1}); 
		}
	}
}
sub posi2{
	%h1=();%h2=();
	my $spe1=$_[0];	my $spe2=$_[1];
	open B1,"<$opts{i}/$spe1.gff";
		while ($a=<B1>){
		chomp $a;
		my @items1=split/\t/,$a;
		$h1{$items1[1]}=$items1[0]."\t".$items1[2]."\t".$items1[3];
	}
	close B1;
	open B2,"<$opts{i}/$spe2.gff";
	while ($a2=<B2>){
		chomp $a2;
		my @items2=split/\t/,$a2;
		$h2{$items2[1]}=$items2[0]."\t".$items2[2]."\t".$items2[3];
	}
	close B2;

	open O,">$opts{o}/microsyn_genes_links/$spe1\_$spe2\_microsyn_genes.links";
	open I,"<$opts{o}/MCScanX_result/$spe1\_$spe2.collinearity";
	while ($b=<I>){
		chomp $b;
		if ($b!~/^#/) {
			my @items3=split/\t/,$b;
			if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
				print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=chr6\n";
			}
			if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
				print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=chr6\n";
			}
		}
	}
	close I;
	close O;
}
system "cat $opts{o}/microsyn_genes_links/*_microsyn_genes.links > $opts{o}/microsyn_genes_links/All_species_microsyn_genes.links";




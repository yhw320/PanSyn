#!/usr/bin/perl-w
#Author:Hongwei Yu
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","m=s","t=s","ca=s","ci=s","e=s","o=s","h|help");
if(!(defined $opts{i} and defined $opts{o})){

	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder storing output files
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-ca	The color of the lines of macro-synteny conserved genes in the resulting graph (default: chr5)
	-ci	The color of the lines of micro-synteny genes in the resulting graph (default: chr24)
	-m	Full path where the MCScanX software command resides (/mnt/inspurfs/home/liyl/yuhw/software/MCScanX-master/MCScanX)
	-h|-help Print this help page
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder storing output files
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-ca	The color of the lines of macro-synteny conserved genes in the resulting graph (default: chr5)
	-ci	The color of the lines of micro-synteny genes in the resulting graph (default: chr24)
	-m	Full path where the MCScanX software command resides (/mnt/inspurfs/home/liyl/yuhw/software/MCScanX-master/MCScanX)
	-h|-help Print this help page
		*************************************************\n";
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{e})) {
	$opts{e}=1e-3;
}
if (!(defined $opts{ca})) {
	$opts{ca}=chr5;
}
if (!(defined $opts{ci})) {
	$opts{ci}=chr24;
}
if (!(defined $opts{m})) {
	$opts{m}="MCScanX";
}

$n=1;
my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ($c!~/_/) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $file=$dir1."/$c";
			system "diamond makedb --in $file -d $file.database";
			$h{$n}=$1;
			$n=$n+1;
		}
	}
}
close P;

$num=keys%h;
system "mkdir $opts{o}/microsyn_genes.links";
system "mkdir $opts{o}/collinearity";
foreach my $i (sort {$a<=>$b} keys %h) {
	for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
		if ($i!=$num) {
			system "mkdir $opts{o}/$h{$i}_$h{$n1}";
	        system "diamond blastp --query $opts{i}/$h{$i}.pep --db $opts{i}/$h{$n1}.pep.database --outfmt 6 --out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs 5";
			system "cat $opts{i}/$h{$i}.gff $opts{i}/$h{$n1}.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
			chdir "$opts{o}/$h{$i}_$h{$n1}";
			system "$opts{m} $h{$i}_$h{$n1}";
			system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/collinearity";
			color($h{$i},$h{$n1}); 
			posi2($h{$i},$h{$n1}); 
		}
	}
}
system "cat $opts{o}/microsyn_genes.links/*_microsyn_genes.links > $opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
sub color{
	my $spe1=$_[0];	my $spe2=$_[1];%h_gene1=();%h_gene2=();
	open GENE1,"<$opts{i}/$spe1\_diagonal_gene_pairs.pep";
	while (my $a=<GENE1>){
		chomp $a;
		if ($a=~/>(\S+)/) {
			$h_gene1{$1}=1;
		}
	}
	close GENE1;
	open GENE2,"<$opts{i}/$spe2\_diagonal_gene_pairs.pep";
	while (my $b=<GENE2>){
		chomp $b;
		if ($b=~/>(\S+)/) {
			$h_gene2{$1}=1;
		}
	}
	close GENE2;
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
	open O,">$opts{o}/microsyn_genes.links/$spe1\_$spe2\_microsyn_genes.links";
	open I,"<$opts{o}/collinearity/$spe1\_$spe2.collinearity";
	while ($b=<I>){
		chomp $b;
		if ($b!~/^#/) {
			my @items3=split/\t/,$b;
			if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
				if (exists $h_gene1{$items3[1]} and exists $h_gene2{$items3[2]}) {
					print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ca}\n";
				}
				else{
					if (exists $h_gene1{$items3[1]} or exists $h_gene2{$items3[2]}) {
						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ci}\n";
					}
					else{
						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ci}\n";
					}
				}
			}
			if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
				if (exists $h_gene2{$items3[1]} and exists $h_gene1{$items3[2]}) {
					print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ca}\n";
				}
				else{
					if (exists $h_gene2{$items3[1]} or exists $h_gene1{$items3[2]}) {
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ci}\n";
					}
					else{
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ci}\n";
					}
				}
			}
		}
	}
	close I;
	close O;
}
#micro_highlight
open I,"<$opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
open O,">$opts{o}/microsyn_genes.links/1";
while (my $a=<I>){
	chomp $a;
	my @yu=split/\t/,$a;
	print O "$yu[0]\t$yu[1]\t$yu[2]\n";
	print O "$yu[3]\t$yu[4]\t$yu[5]\n";
}
close I;close O;
system "sort -k 1,1 -k 2n,2 $opts{o}/microsyn_genes.links/1| uniq >$opts{o}/microsyn_genes.links/micro_highlight.txt";
system "rm $opts{o}/microsyn_genes.links/1";

#macro_highlight
open O,">$opts{o}/microsyn_genes.links/macro_highlight.txt";
my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+)_diagonal_gene_pairs.pep$/ ) {
		$name=$1;
		my $file=$dir1."/$c";
		open I,"<$file";
		open I2,"<$opts{i}/$name.gff";
		while (my $a=<I2>){
			chomp $a;
			my @ui=split/\t/,$a;
			$gene_chr_p{$ui[1]}=$ui[0]."\t".$ui[2]."\t".$ui[3];
		}
		while (my $a=<I>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				if (exists $gene_chr_p{$1}) {
					print O "$gene_chr_p{$1}\n";
				}
			}
		}
		close I;
		close I2;
	}
}
close P;
close O;

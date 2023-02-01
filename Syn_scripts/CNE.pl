#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","b=s","r=s","p=s","o=s","e=s","l=s","t=s","height=s","width=s","z=s","s=s");
if(!(defined $opts{i} and defined $opts{b} and defined $opts{r} and defined $opts{o} and defined $opts{p} and defined $opts{s})){

	die"**********************************************\n
		-i	Full path to the [inputDir] folder containing input files
		-b	Full path to the [Species_cluster.bed] file
		-r	Reference species by name abbreviations (example: HSap)
		-o	Full path to the [outputDir] folder containing output files
		-p	Full path to the [*_gene.pos] file
		-s	Full path to the [syn_scripts] file PanSyn provided
		Optional:
		-z	Orientation of gene cluster display (1 or 2)
		-e	Threshold of similarity between sequences (0-1].(default: 0.95)
		-l	Minimum length of CNE (default: 50)
		-t	Number of threads to use (default: 12)
		-height	The height of the resulting plot (default: 7)
		-width	The width of the resulting plot (default: 7)
	*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
	die "************************************************\n
		-i	Full path to the [inputDir] folder containing input files
		-b	Full path to the [Species_cluster.bed] file
		-r	Reference species by name abbreviations (example: HSap)
		-o	Full path to the [outputDir] folder containing output files
		-p	Full path to the [*_gene.pos] file
		-s	Full path to the [syn_scripts] file PanSyn provided
		Optional:
		-z	Orientation of gene cluster display (1 or 2)
		-e	Threshold of similarity between sequences (0-1].(default: 0.95)
		-l	Minimum length of CNE (default: 50)
		-t	Number of threads to use (default: 12)
		-height	The height of the resulting plot (default: 7)
		-width	The width of the resulting plot (default: 7)
		*************************************************\n";
}

if (!(defined $opts{e})) {
	$opts{e}=0.95;
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{l})) {
	$opts{l}=50;
}
if (!(defined $opts{height})) {
	$opts{height}=7;
}
if (!(defined $opts{width})) {
	$opts{width}=7;
}
if (!(defined $opts{z})) {
	$opts{z}=1;
}

my %h_cluster_pos=();
my %h_spe_chr=();
my %h_spename=();
my %h_reverse=();

open I,"<$opts{b}";
while (my $a=<I>){
	chomp $a;
	my @items=split/\t/,$a;
	$h_cluster_pos{$items[0]}=$items[2]."\t".$items[3];
	$h_spe_chr{$items[0]}=$items[1];
	$h_reverse{$items[0]}=$items[4];
}
close I;

my $dir1=$opts{i};
opendir P,$dir1;
my $num_spe=1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).fa$/ ) {
		if ($1 ne $opts{r}) {
			$a_b=$h_cluster_pos{$opts{r}};
			$c_d=$h_cluster_pos{$1};
			my @i_ab=split/\t/,$a_b;
			my @i_cd=split/\t/,$c_d;
			$h_spename{$1}=0;
			$num_spe=$num_spe+1;
			chdir "$opts{o}";
			system "cnef -r $opts{i}/$opts{r}.fa -q $opts{i}/$1.fa -e $opts{i}/$opts{r}.exon -f $opts{i}/$1.exon -y $h_spe_chr{$opts{r}} -z $h_spe_chr{$1} -a $i_ab[0] -b $i_ab[1] -c $i_cd[0] -d $i_cd[1] -t $opts{e} -l $opts{l} -o $opts{r}\_$1.cne -T $opts{t} -v $h_reverse{$1}";
		}
	}
}
close P;
############ start plot#################
system "mkdir $opts{o}/for_R";
open O,">$opts{o}/for_R/cne_pos.txt";
open O1,">$opts{o}/for_R/formatted_cne_pos.txt";
foreach my $i (keys %h_spename) {
	open I,"<$opts{o}/$opts{r}_$i.cne";
	while (my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		print O "$i\t$it[1]\t$it[2]\n";
		$color1="\"#008DF6\"";#blue
		print O1 "$h_spe_chr{$opts{r}}\t$it[1]\t$it[2]\t$color1\t$i\n";
	}
}
close O;
close O1;

if ($num_spe==2) {
	open I,"<$opts{p}";
	open O,">$opts{o}/for_R/gene_pos.txt";
	$color1="\"#FF7DBE\"";
	$gene="gene";
	while (my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		print O "$it[0]\t$it[1]\t$it[2]\t$color1\t$gene\n";
		$chr=$it[0];
	}
	close I;
	close O;
	system "cat $opts{o}/for_R/gene_pos.txt $opts{o}/for_R/formatted_cne_pos.txt >$opts{o}/for_R/all-for-R.txt";
}
if ($num_spe>2) {
	system "mkdir $opts{o}/for_common";
	foreach my $i (keys %h_spename) {
		open O,">$opts{o}/for_common/$i.pos";
		open I,"<$opts{o}/$opts{r}_$i.cne";
		while (my $a=<I>){
			chomp $a;
			my @it=split/\t/,$a;
			print O "$it[0]\t$it[1]\t$it[2]\n";
		}
		close O;
	}

	$n=0;
	my $dir2="$opts{o}/for_common";
	opendir P,$dir2;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pos$/ ) {
			$output="$opts{o}/for_common/$1.pos";
			$n=$n+1;
			if ($n==1) {
				$out_before=$output;
				$spe_name=$1;
			}
			else{
				system "bedtools intersect -a $out_before -b $output > $opts{o}/for_common/$spe_name\_$1.pos";
				$out_before="$opts{o}/for_common/$spe_name\_$1.pos";
				$common_file="$opts{o}/for_common/$spe_name\_$1.pos";
				$spe_name=$spe_name."_".$1;
			}
		}
	}
	close P;

	system "cp $common_file $opts{o}/allspe_common_cne.txt";

	open I,"<$opts{p}";
	open O,">$opts{o}/for_R/gene_pos.txt";
	$color1="\"#FF7DBE\"";
	$gene="gene";
	while (my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		print O "$it[0]\t$it[1]\t$it[2]\t$color1\t$gene\n";
		$chr=$it[0];
	}
	close I;
	$color2="\"#F2A672\"";
	$co="Common";
	open I,"<$opts{o}/allspe_common_cne.txt";
	while (my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		print O "$chr\t$it[1]\t$it[2]\t$color2\t$co\n";
	}
	close I;
	close O;
	system "cat $opts{o}/for_R/gene_pos.txt $opts{o}/for_R/formatted_cne_pos.txt >$opts{o}/for_R/all-for-R.txt";

	system "rm -r $opts{o}/for_common";
}

system "Rscript $opts{s}/CNE.R $opts{o}/for_R $opts{o}/for_R/all-for-R.txt $opts{z} $chr $opts{p} $opts{height} $opts{width}";

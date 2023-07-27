#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","o2=s","b=s","c=s","color=s","width=s","height=s","h|help");
if(!(defined $opts{i} and defined $opts{o2} and defined $opts{b} and defined $opts{c} and defined $opts{color})){
	die"**********************************************\n
	-i	The parameter refers to the same parameter [-i] used in the previously executed script [Microsyn_cluster.pl]
	-o2	Full path to the new [outputDir] directory containing output files	
	-b	Full path to the interested [*_microsyn_genes.blocks] directory generated by the previous script [Microsyn_cluster.pl]
	-c	Full path to the [A_single.cluster] file
	-color	Full path to the [genes_color.txt] file
	Optional:
	-width	The width of the output graph (default:8)
	-height	The height of the output graph (default:3)
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
	die"**********************************************\n
	-i	The parameter refers to the same parameter [-i] used in the previously executed script [Microsyn_cluster.pl]
	-o2	Full path to the new [outputDir] directory containing output files	
	-b	Full path to the interested [*_microsyn_genes.blocks] directory generated by the previous script [Microsyn_cluster.pl]
	-c	Full path to the [A_single.cluster] file
	-color	Full path to the [genes_color.txt] file
	Optional:
	-width	The width of the output graph (default:8)
	-height	The height of the output graph (default:3)
		*********************************************\n";
}
my $one_lie; my $two_lie;
my $slash; 
####
if ($opts{i} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{i} =~ s/$slash$//;
}
####

####
if ($opts{o2} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o2} =~ s/$slash$//;
}
####

####
if ($opts{b} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{b} =~ s/$slash$//;
}
####

if (!(defined $opts{width})) {
	$opts{width}=8;
}
if (!(defined $opts{height})) {
	$opts{height}=3;
}

######################

if (-d "$opts{o2}/visualize_cluster") {
	system "rm -r $opts{o2}/visualize_cluster";
	print "The previous [$opts{o2}/visualize_cluster] directory has been removed!\n";
} 

#######################

open I,"<$opts{c}" or die("Could not open $opts{c}.\n"); 
%h_PYes=();
while ($a=<I>){
	chomp $a;
	$h_PYes{$a}=1;
}
close I;
system "mkdir $opts{o2}/visualize_cluster";
print "The [$opts{o2}/visualize_cluster] directory has been created!\n";

$dir1="$opts{b}/conserved_clusters";
opendir I,$dir1;
%h_duiying=();
while ($aa=readdir I) {
	if ( $aa=~/^(\S+)_conserved.clusters$/ ) {
		$two_spe_name=$1;
		my @spe_name=split/_/,$two_spe_name;
		$spe1=$spe_name[0];
		$spe2=$spe_name[1];
		$file=$dir1."/$aa";
		open IN,"<","$file" or die;
		$n=0;
		while (my $b=<IN>) {
			chomp $b;
			if ($b=~/^#/) {
				if ($n==0) {
					$one_lie="haha";
					$two_lie="haha";
				}
				else{
					if (exists $h_PYes{$one_lie}) {
						$h_duiying{$spe2}=$two_lie;
						$h_duiying{$spe1}=$one_lie;
					}
					$one_lie="haha";
					$two_lie="haha";
				}
				$n=$n+1;
			}
			else {
				my @it2=split/\t/,$b;
				if ($one_lie eq "haha" and $two_lie eq "haha") {
					$one_lie=$it2[0];
					$two_lie=$it2[1];
				}
				else{
					$one_lie=$one_lie."\t".$it2[0];
					$two_lie=$two_lie."\t".$it2[1];
				}
			}
		}
		close IN;
	}
}
close I;
open O,">$opts{o2}/visualize_cluster/allspe_conserved.clusters";
foreach my $i (keys %h_duiying) {
	print O "$i\t$h_duiying{$i}\n";
}
close O;


###hox2.pl
open IN,"<$opts{o2}/visualize_cluster/allspe_conserved.clusters";
open O,">$opts{o2}/visualize_cluster/forR_microsyn_conserved.cluster";
print O "name\tstart\tend\tstrand\tcol\tfill\tGenome\tgene_type\n";

while (my $a=<IN>){
	chomp $a;
	my @it1=split/\t/,$a;
	open G,"<$opts{i}/$it1[0]_simplified.gff";
	%h_gene=();
	while (my $b=<G>){
		chomp $b;
		my @it2=split/\t/,$b;
		if ($it2[4] eq "-") {
			$ha=-1;
		}
		if ($it2[4] eq "+") {
			$ha=1;
		}
		$h_gene{$it2[1]}=$it2[0]."\t".$it2[2]."\t".$it2[3]."\t".$ha;
	}
	close G;
	system "sort -k1,1 -k 3n,3  $opts{i}/$it1[0]_simplified.gff >$opts{i}/$it1[0].gff_sort";
	open G2,"<$opts{i}/$it1[0].gff_sort";
	%h_chr=();%h_gene_234=();%h_gene_pos=();%h_gene_chr=();%h2_chr_jishu_gene=();
	while (my $c=<G2>){
		chomp $c;
		my @it4=split/\t/,$c;
		if (!exists $h_chr{$it4[0]}) {
			$n_jishu=1;
			$h_chr{$it4[0]}=0;
			$h_gene_234{$it4[1]}=$n_jishu;
			$h_gene_chr{$it4[1]}=$it4[0];
			$h2_chr_jishu_gene{$it4[0]}{$n_jishu}=$it4[1];
			$h_gene_pos{$it4[1]}=$it4[0]."\t".$it4[2]."\t".$it4[3];
		}
		else{
			$h_gene_chr{$it4[1]}=$it4[0];
			$n_jishu=$n_jishu+1;
			$h_gene_234{$it4[1]}=$n_jishu;
			$h2_chr_jishu_gene{$it4[0]}{$n_jishu}=$it4[1];
			$h_gene_pos{$it4[1]}=$it4[0]."\t".$it4[2]."\t".$it4[3];
		}
	}
	close G2;
	foreach my $i (1..$#it1) {
		if (exists $h_gene{$it1[$i]}) {
			my @it3=split/\t/,$h_gene{$it1[$i]};
			print O "$it1[$i]\t$it3[1]\t$it3[2]\t$it3[3]\t$i\t$i\t$it3[0]\tblocks\n";
		}
	}
	$jian=$#it1-1;
	$chr_name=$h_gene_chr{$it1[1]};
	foreach my $i (1..$jian) {
		$ijia=$i+1;
		$zhi=$h_gene_234{$it1[$ijia]}-$h_gene_234{$it1[$i]};
		if ($zhi>1) {
			$start=$h_gene_234{$it1[$i]}+1;
			for ($n3=$start;$n3<$h_gene_234{$it1[$ijia]};$n3=$n3+1) {
				$jiange_gene=$h2_chr_jishu_gene{$chr_name}{$n3};
				my @jiayou=split/\t/,$h_gene_pos{$jiange_gene};
				my @jiayou2=split/\t/,$h_gene{$jiange_gene};
				$ling=0;
				print O "Interval_gene\t$jiayou[1]\t$jiayou[2]\t$jiayou2[3]\t$ling\t$ling\t$jiayou[0]\tblocks\n";
			}
		}
		if ($zhi<-1) {
			$start=$h_gene_234{$it1[$ijia]}+1;
			for ($n3=$start;$n3<$h_gene_234{$it1[$i]};$n3=$n3+1) {
				$jiange_gene=$h2_chr_jishu_gene{$chr_name}{$n3};
				my @jiayou=split/\t/,$h_gene_pos{$jiange_gene};
				my @jiayou2=split/\t/,$h_gene{$jiange_gene};
				$ling=0;
				print O "Interval_gene\t$jiayou[1]\t$jiayou[2]\t$jiayou2[3]\t$ling\t$ling\t$jiayou[0]\tblocks\n";
			}
		}
	}
}
close IN;
close O;


open I1,"<$opts{o2}/visualize_cluster/forR_microsyn_conserved.cluster";
open O,">$opts{o2}/visualize_cluster/inter_forR_microsyn_conserved.cluster";
open I2,"<$opts{color}";
while ($a=<I2>){
	chomp $a;
	my @itt=split/\t/,$a;
	$h_color{$itt[0]}=$itt[1];
}
close I2;

while ($a=<I1>){
	chomp $a;
	my @itt=split/\t/,$a;
	if ($itt[0] eq "name") {
		print O "$a\n";
	}
	else{
		if (exists $h_color{$itt[4]}) {
			print O "$itt[0]\t$itt[1]\t$itt[2]\t$itt[3]\t$h_color{$itt[4]}\t$h_color{$itt[4]}\t$itt[6]\t$itt[7]\n";
		}
		else{
			print O "$itt[0]\t$itt[1]\t$itt[2]\t$itt[3]\tgrey\tgrey\t$itt[6]\t$itt[7]\n";
		}
	}
}
close I1;
close O;

open INN,"<$opts{o2}/visualize_cluster/inter_forR_microsyn_conserved.cluster";
open O,">$opts{o2}/visualize_cluster/final_forR_microsyn_conserved.cluster";
$chr_name_before="hahaha";my %h_bi;my %h_dute;
while ($a=<INN>){
	chomp $a;
	if ($a!~/^name/ and $a!~/gene_type$/) {
	}
	my @itt=split/\t/,$a;
	if ($itt[6] ne $chr_name_before) {
		$h_bi{$itt[6]}=$itt[1];
	}
	else{
		if ($itt[1]<$h_bi{$itt[6]}) {
			$h_dute{$itt[6]}=0;
		}
	}
	$chr_name_before=$itt[6];
}
close INN;
open INN,"<$opts{o2}/visualize_cluster/inter_forR_microsyn_conserved.cluster";
while ($a=<INN>){
	chomp $a;
	my @itt=split/\t/,$a;
	if (exists $h_dute{$itt[6]}) {
		print O "$itt[0]\t-$itt[2]\t-$itt[1]\t$itt[3]\t$itt[4]\t$itt[5]\t$itt[6]\t$itt[7]\n";
	}
	else{
		print O "$a\n";
	}
}
close O;
close INN;

system "Rscript Microsyn_visualize_cluster.R $opts{o2}/visualize_cluster $opts{o2}/visualize_cluster/final_forR_microsyn_conserved.cluster $opts{width} $opts{height}";
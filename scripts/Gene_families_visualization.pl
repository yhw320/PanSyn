#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
use List::Util qw(max);
use Cwd 'getcwd';

GetOptions(\%opts,"i1=s","o1=s","width=s","height=s","h|help");
if (!( defined $opts{i1} and defined $opts{o1})) {
		die "************************************************\n
	-i1	The parameter refers to the same parameter [-i] used in the previously executed command [Macrosyn1]
	-o1	Full path to the [outputDir1] directory containing output files (same as [-o1] in the previous command [Macrosyn1])
	Optional:
	-width	Set the width of the output image (default:18)
	-height	Set the height of the output image (default:9)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
Options:
	-i1	The parameter refers to the same parameter [-i] used in the previously executed command [Macrosyn1]
	-o1	Full path to the [outputDir1] directory containing output files (same as [-o1] in the previous commad [Macrosyn1])
	Optional:
	-width	Set the width of the output image (default:18)
	-height	Set the height of the output image (default:9)
	-h|-help	Print this help page
		*************************************************\n";
}

my $slash; 
####
if ($opts{i1} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{i1} =~ s/$slash$//;
}
####
####
if ($opts{o1} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o1} =~ s/$slash$//;
}
####


if (!(defined $opts{width})) {
	$opts{width}=18;
}
if (!(defined $opts{height})) {
	$opts{height}=9;
}

my $gene_group;

open I,"<$opts{o1}/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene_families.analysis/ancient_gene.families.\n"); 
my %h_all_spe;
while (my $a=<I>){
	chomp $a;
	my @items1=();my @items2=();my @items3=();
	@items1=split/\t/,$a;
    $gene_group=$items1[1];
	@items2=split/,/,$gene_group;
	foreach my $gene_id (@items2) {
		@items3=split/_/,$gene_id;
		$h_all_spe{$items3[0]}=0;
	}
}
close I;
my $all_spe_num = scalar keys %h_all_spe;


#barchart
open I,"<$opts{o1}/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene_families.analysis/ancient_gene.families.\n"); 

my %h_cluster_spe;my %h_cluster_bizhi;

while (my $a=<I>){
	chomp $a;
	my @items1=();my @items2=();my @items3=();
	@items1=split/\t/,$a;
    $gene_group=$items1[1];
	@items2=split/,/,$gene_group;
	foreach my $gene_id (@items2) {
		@items3=split/_/,$gene_id;
		$h_cluster_spe{$items3[0]}=0;
	}
	my $cluster_spe_num = scalar keys %h_cluster_spe;
	$h_cluster_bizhi{$items1[0]}=$cluster_spe_num/$all_spe_num;
	undef %h_cluster_spe;
}
close I;

open O,"> $opts{o1}/barchart_data";

open I,"<$opts{o1}/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene_families.analysis/ancient_gene.families.\n"); 

my %h_spe_cluster;

while (my $a=<I>){
	chomp $a;
	my @items1=();my @items2=();my @items3=();
	@items1=split/\t/,$a;
    $gene_group=$items1[1];
	@items2=split/,/,$gene_group;
	foreach my $gene_id (@items2) {
		@items3=split/_/,$gene_id;
		$h_spe_cluster{$items3[0]}{$items1[0]}=$h_all_spe{$items3[0]};
	}
}
close I;

my %h_spe_jishu1;my %h_spe_jishu2;my %h_spe_jishu3;my %h_spe_jishu4;

foreach my $spe (keys %h_spe_cluster) {
	foreach my $cluster (keys %{$h_spe_cluster{$spe}}) {
		if (exists $h_cluster_bizhi{$cluster}) {
			if (0<=$h_cluster_bizhi{$cluster} and $h_cluster_bizhi{$cluster}<0.25) {
				if (exists $h_spe_jishu1{$spe}) {
					$h_spe_jishu1{$spe}=$h_spe_jishu1{$spe}+1;
				}
				else{
					$h_spe_jishu1{$spe}=1;
				}
			}
			if (0.25<=$h_cluster_bizhi{$cluster} and $h_cluster_bizhi{$cluster}<0.5) {
				if (exists $h_spe_jishu2{$spe}) {
					$h_spe_jishu2{$spe}=$h_spe_jishu2{$spe}+1;
				}
				else{
					$h_spe_jishu2{$spe}=1;
				}
			}
			if (0.5<=$h_cluster_bizhi{$cluster} and $h_cluster_bizhi{$cluster}<0.75) {
				if (exists $h_spe_jishu3{$spe}) {
					$h_spe_jishu3{$spe}=$h_spe_jishu3{$spe}+1;
				}
				else{
					$h_spe_jishu3{$spe}=1;
				}
			}
			if (0.75<=$h_cluster_bizhi{$cluster} and $h_cluster_bizhi{$cluster}<=1) {
				if (exists $h_spe_jishu4{$spe}) {
					$h_spe_jishu4{$spe}=$h_spe_jishu4{$spe}+1;
				}
				else{
					$h_spe_jishu4{$spe}=1;
				}
			}

		}
	}
}




foreach my $spe (keys %h_all_spe) {
	print O "$spe\t";
	print O "$h_spe_jishu4{$spe}\t$h_spe_jishu3{$spe}\t$h_spe_jishu2{$spe}\t$h_spe_jishu1{$spe}\n";
}

close O;


open O3,"> $opts{o1}/inputR_barchart_data";
open I2,"<$opts{o1}/barchart_data";
print O3 "Species\tNumber\tRange\tPosition\n";
while ($a=<I2>){
	chomp $a;
	my @it=split/\t/,$a;
	my $s1=$it[1]/2;
	my $s2=$it[1]+$it[2]/2;
	my $s3=$it[1]+$it[2]+$it[3]/2;
	my $s4=$it[1]+$it[2]+$it[3]+$it[4]/2;
	print O3 "$it[0]\t$it[1]\t75%-100%\t$s1\n";
	print O3 "$it[0]\t$it[2]\t50%-75%\t$s2\n";
	print O3 "$it[0]\t$it[3]\t25%-50%\t$s3\n";
	print O3 "$it[0]\t$it[4]\t0-25%\t$s4\n";
}
close O3;



my $old_dir = getcwd();

system("Gene_families $opts{o1}/ancient_gene_families.analysis/family_situation.table  $opts{o1}/inputR_barchart_data $opts{i1}/tree.nwk $opts{o1} $opts{i1}/species_classification.table $opts{width} $opts{height}");
#system("rm $opts{o1}/ancient_gene_families.analysis/family_situation.table");
system("rm $old_dir/Rplots.pdf");
system("mv $opts{o1}/inputR_barchart_data $opts{o1}/ancient_gene_families.analysis/inputR_barchart_data");
system("mv $opts{o1}/Ancient_gene_families.pdf $opts{o1}/ancient_gene_families.analysis/Ancient_gene_families.pdf");

#!/usr/bin/perl
use Getopt::Long;
my %opts;
use List::Util qw(max);
GetOptions(\%opts,"i=s","s=s","o=s","width=s","height=s","h|help");
if (!( defined $opts{i} and defined $opts{s} and defined $opts{o})) {
		die "************************************************\n
	-i	Full path to the folder containing required data files (same as [-i] in the previous script [Macrosyn1_1_blast.pl])
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the folder storing the output results (same as [-o] in the previous script [Macrosyn1_1_blast.pl])
	Optional:
	-width	Set the width of the output picture (default:18)
	-height	Set the height of the output picture (default:9)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
Options:
	-i	Full path to the folder containing required data files (same as [-i] in the previous script [Macrosyn1_1_blast.pl])
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the folder storing the output results (same as [-o] in the previous script [Macrosyn1_1_blast.pl])
	Optional:
	-width	Set the width of the output picture (default:18)
	-height	Set the height of the output picture (default:9)
	-h|-help	Print this help page
		*************************************************\n";
}
if (!(defined $opts{width})) {
	$opts{width}=18;
}
if (!(defined $opts{height})) {
	$opts{height}=9;
}

#barchart
open I,"<$opts{o}/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $opts{o}/ancient_gene_families.analysis/ancient_gene.families.\n"); 
open O,"> $opts{o}/barchart_data";
@items1=();@items2=();@items3=();

while ($a=<I>){
	chomp $a;
	@items1=split/\t/,$a;
    $gene_group=$items1[1];
	@items2=split/,/,$gene_group;
	foreach $gene_id (@items2) {
		@items3=split/_/,$gene_id;
		$h22{$items3[0]}=1;
		$h_spe_name{$items3[0]}=0;
		$h12{$items3[0]}{$items1[0]}=0;
	}
	@keys=keys %h22;
	$size=@keys;
	$h33{$items1[0]}=$size;
	undef %h22;
	undef @items2;
	undef @items1;
}
$spe_num=0;
foreach my $i (keys %h_spe_name) {
	$spe_num=$spe_num+1;
	$h6{$i}=0;
	$h7{$i}=0;
	$h8{$i}=0;
	$h9{$i}=0;
}
$num_n1=($spe_num/4)*1;
$num_n2=($spe_num/4)*2;
$num_n3=($spe_num/4)*3;
$num_n4=($spe_num/4)*4;

foreach $gene (keys %h12) {
	foreach $fam_num (keys %{$h12{$gene}}) {
		if (exists $h33{$fam_num}) {
			$nana=$gene."_".$fam_num;
			$h44{$nana}=$h33{$fam_num};
		}
	}
}
foreach $gene2 (keys %h44) {
	$num=$h44{$gene2};
	#print "$gene2\t$num\n";
	@items4=split/_/,$gene2;
	$spe_name=$items4[0];
	if (0<=$num and $num<$num_n1) {
		$h5{$spe_name}{1}=$h6{$spe_name}+1;
		$h6{$spe_name}=$h5{$spe_name}{1};
	}
	if ($num_n1<=$num and $num<$num_n2) {
		$h5{$spe_name}{2}=$h7{$spe_name}+1;
		$h7{$spe_name}=$h5{$spe_name}{2};
	}
	if ($num_n2<=$num and $num<=$num_n3) {
		$h5{$spe_name}{3}=$h8{$spe_name}+1;
		$h8{$spe_name}=$h5{$spe_name}{3};
	}
	if ($num_n3<=$num and $num<=$num_n4) {
		$h5{$spe_name}{4}=$h9{$spe_name}+1;
		$h9{$spe_name}=$h5{$spe_name}{4};
	}
}
foreach my $i (keys %h_spe_name) {
	if (exists $h5{$i} and !exists $h5{$i}{1}) {
		$h5{$i}{1}=0;
	}
	if (exists $h5{$i} and !exists $h5{$i}{2}) {
		$h5{$i}{2}=0;
	}
	if (exists $h5{$i} and !exists $h5{$i}{3}) {
		$h5{$i}{3}=0;
	}
	if (exists $h5{$i} and !exists $h5{$i}{4}) {
		$h5{$i}{4}=0;
	}
}

foreach $key1 (keys %h5) {
	print O "$key1\t";
	foreach $key2 (keys %{$h5{$key1}}) {
		print O "$h5{$key1}{$key2}\t";
	}
	print O "\n";
}
close I;
close O;


open O3,"> $opts{o}/inputR_barchart_data";
open I2,"<$opts{o}/barchart_data";
print O3 "Species\tNumber\tRange\tPosition\n";
while ($a=<I2>){
	chomp $a;
	@it=split/\t/,$a;
	$s1=$it[1]/2;
	$s2=$it[1]+$it[2]/2;
	$s3=$it[1]+$it[2]+$it[3]/2;
	$s4=$it[1]+$it[2]+$it[3]+$it[4]/2;
	print O3 "$it[0]\t$it[1]\t75%-100%\t$s1\n";
	print O3 "$it[0]\t$it[2]\t50%-75%\t$s2\n";
	print O3 "$it[0]\t$it[3]\t25%-50%\t$s3\n";
	print O3 "$it[0]\t$it[4]\t0-25%\t$s4\n";
}
close O3;




system("Rscript $opts{s}/Gene_families.R $opts{o}/ancient_gene_families.analysis/family_situation.table  $opts{o}/inputR_barchart_data $opts{i}/tree.nwk $opts{o} $opts{i}/species_classification.table $opts{width} $opts{height}");
#system("rm $opts{o}/ancient_gene_families.analysis/family_situation.table");
system("rm $opts{o}/Rplots.pdf");
system("mv $opts{o}/inputR_barchart_data $opts{o}/ancient_gene_families.analysis/inputR_barchart_data");
system("mv $opts{o}/Ancient_gene_families.pdf $opts{o}/ancient_gene_families.analysis/Ancient_gene_families.pdf");

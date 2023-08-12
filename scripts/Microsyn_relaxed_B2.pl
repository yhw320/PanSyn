#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"pansyn=s","o2=s","node=s","n=s","tre=s","pl=s","ps=s","pk=s","pr=s","spenum2=s","h|help");
if (!( defined $opts{o2} and defined $opts{pansyn} and defined $opts{node} and defined $opts{tre} and defined $opts{n})) {
		die "************************************************\n
	-o2	Full path to the [outputDir2_S30] directory
	-node	The name of the evolution nodes of interest (e.g. Ancestor1)
	-tre	Full path to the [*.tre] file
	-n	The distance threshold above which syntenic blocks will be split (the value in the x-left column of the [nmax_test.inflection.csv])
	-pansyn Full path the [scripts] provided by PanSyn
	optional:
	-pl	Minimum number of OG co-occurrences between two extant species blocks (default: 3)
	-ps	Minimum overlap coefficient between the OG occurrences of two extant species blocks (scaffolds, default: 3)
	-pk	Minimum number of edges per community/multispecies block (if k = 3 then a multispecies block is only retained if it is conserved across at least three species, default: 3)
	-pr	Percentage of OGs of the ancestral syntenic block that an extant species block should possess to be retained (default: 0.3)
	-spenum2	Minimum number of species of a phylogenetic clade that must possess a syntenic OG pair (default: 3)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-o2	Full path to the [outputDir2_S30] directory
	-node	The name of the evolution nodes of interest (e.g. Ancestor1)
	-tre	Full path to the [*.tre] file
	-n	The distance threshold above which syntenic blocks will be split (the value in the x-left column of the [nmax_test.inflection.csv])
	-pansyn Full path the [scripts] provided by PanSyn
	optional:
	-pl	Minimum number of OG co-occurrences between two extant species blocks (default: 3)
	-ps	Minimum overlap coefficient between the OG occurrences of two extant species blocks (scaffolds, default: 3)
	-pk	Minimum number of edges per community/multispecies block (if k = 3 then a multispecies block is only retained if it is conserved across at least three species, default: 3)
	-pr	Percentage of OGs of the ancestral syntenic block that an extant species block should possess to be retained (default: 0.3)
	-spenum2	Minimum number of species of a phylogenetic clade that must possess a syntenic OG pair (default: 3)
	-h|-help Print this help page
		*************************************************\n";
}

my $slash; 
####
if ($opts{o2} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o2} =~ s/$slash$//;
}
####

####
if ($opts{pansyn} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{pansyn} =~ s/$slash$//;
}
####


if (!(defined $opts{pl})) {
	$opts{pl}=3;
}
if (!(defined $opts{ps})) {
	$opts{ps}=0.5;
}
if (!(defined $opts{pk})) {
	$opts{pk}=3;
}
if (!(defined $opts{pr})) {
	$opts{pr}=0.3;
}
if (!(defined $opts{spenum2})) {
	$opts{spenum2}=3;
}

#chdir "$opts{o2}";
system "python $opts{pansyn}/step3_find_og_commus.py $opts{o2}/$opts{node}-myresults/bynode/$opts{node}/m_$opts{spenum2}/$opts{node}-myresults.m_$opts{spenum2}.$opts{node}.dist  -n $opts{n} -o $opts{o2}/$opts{node}-myresults/bynode/$opts{node}.og_commus";

system "python $opts{pansyn}/step4_OG_communities_to_blocks_graph_check.py $opts{o2}/$opts{node}-myresults/bynode/$opts{node}.og_commus.csv -g $opts{o2}/$opts{node}-myresults/bynode/$opts{node}.og_commus.gpickle -c $opts{o2}/$opts{node}-myresults/chromdata.pickle -o $opts{o2}/$opts{node}-myresults/$opts{node}_blocks -l $opts{pl} -s $opts{ps} -k $opts{pk} -r $opts{pr}";

system "python $opts{pansyn}/BlocksByNode.py -c $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.clusters -b $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.synt -s $opts{tre} -n $opts{node} -m $opts{spenum2} -r blocks_list -t total > $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_total.synt";
system "python $opts{pansyn}/BlocksByNode.py -c $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.clusters -b $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.synt -s $opts{tre} -n $opts{node} -m $opts{spenum2} -r blocks_list -t novel > $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_novel.synt";
system "python $opts{pansyn}/BlocksByNode.py -c $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.clusters -b $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.synt -s $opts{tre} -n $opts{node} -m $opts{spenum2} -r blocks_list -t ancestral > $opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_ancestral.synt";


##
#system "python $opts{s}/analysis_intervening_genes.py -c $opts{node}-myresults/chromdata.pickle -sy $opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_total.synt -ms $opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.clusters -o $opts{node}-myresults/$opts{node}_blocks.intervening";
##

open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_total.synt" or die;
open O,">","$opts{o2}/$opts{node}-myresults/$opts{node}.summary" or die;

my %h1=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h1{$it[0]}=1;
}
close I;

open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_novel.synt" or die;
my %h2=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h2{$it[0]}=1;
}
close I;

open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len$opts{pl}.ol$opts{ps}.taxonomy_filtered_ancestral.synt" or die;
my %h3=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h3{$it[0]}=1;
}
close I;

my $tatol=0;
foreach my $i (keys %h1) {
	$tatol=$tatol+1;
}

my $novel=0;
foreach my $i (keys %h2) {
	$novel=$novel+1;
}

my $ancestral=0;
foreach my $i (keys %h3) {
	$ancestral=$ancestral+1;
}

print O "tatol_number=$tatol\nnovel_number=$novel\nancestral_number=$ancestral\n";
close O;






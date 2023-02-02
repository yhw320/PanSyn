#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"o2=s","s=s","node=s","n=s","tre=s","h|help");
if (!( defined $opts{o2} and defined $opts{s} and defined $opts{node} and defined $opts{tre} and defined $opts{n})) {
		die "************************************************\n
	-o2	Full path to the new folder containing results (same as [-o2] in the previous script [Ancient_microsyntenic_blocks 1.pl])
	-s	Full path to the [syn_scripts] file PanSyn provided
	-node	Node_names: space-separated names of the nodes for which you want to build ancestral OGs networks.
	-tre	Full path to the tree file
	-n	distance threshold above which syntenic blocks will be considered split. (Can be estimated using  script).
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-o2	Full path to the new folder containing results (same as [-o2] in the previous script [Ancient_microsyntenic_blocks 1.pl])
	-s	Full path to the [syn_scripts] file PanSyn provided
	-node	Node_names: space-separated names of the nodes for which you want to build ancestral OGs networks.
	-tre	Full path to the tree file
	-n	distance threshold above which syntenic blocks will be considered split. (Can be estimated using  script).
	-h|-help Print this help page
		*************************************************\n";
}

chdir "$opts{o2}";
system "python $opts{s}/step3_find_og_commus.py $opts{node}-myresults/bynode/$opts{node}/m_3/$opts{node}-myresults.m_3.$opts{node}.dist  -n $opts{n} -o $opts{node}-myresults/bynode/$opts{node}.og_commus";
system "python $opts{s}/step4_OG_communities_to_blocks_graph_check.py $opts{node}-myresults/bynode/$opts{node}.og_commus.csv -g $opts{node}-myresults/bynode/$opts{node}.og_commus.gpickle -c $opts{node}-myresults/chromdata.pickle -o $opts{node}-myresults/$opts{node}_blocks";
system "python $opts{s}/BlocksByNode.py -c $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.clusters -b $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.synt -s $opts{tre} -n $opts{node} -m 3 -r blocks_list -t total > $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered_total.synt";
system "python $opts{s}/BlocksByNode.py -c $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.clusters -b $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.synt -s $opts{tre} -n $opts{node} -m 3 -r blocks_list -t novel > $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered_novel.synt";
system "python $opts{s}/BlocksByNode.py -c $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.clusters -b $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.synt -s $opts{tre} -n $opts{node} -m 3 -r blocks_list -t ancestral > $opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered.synt_ancestral.synt";


open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered_total.synt" or die;
open O,">","$opts{o2}/$opts{node}-myresults/$opts{node}.summary" or die;

my %h1=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h1{$it[0]}=1;
}
close I;

open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered_novel.synt" or die;
my %h2=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h2{$it[0]}=1;
}
close I;

open I,"<","$opts{o2}/$opts{node}-myresults/$opts{node}_blocks.len3.ol0.5.taxonomy_filtered.synt_ancestral.synt" or die;
my %h3=();
while(my $a=<I>){
	chomp $a;
	my @it=split/\t/,$a;
	$h3{$it[0]}=1;
}
close I;

$tatol=0;
foreach my $i (keys %h1) {
	$tatol=$tatol+1;
}

$novel=0;
foreach my $i (keys %h2) {
	$novel=$novel+1;
}

$ancestral=0;
foreach my $i (keys %h3) {
	$ancestral=$ancestral+1;
}

print O "tatol_number=$tatol\nnovel_number=$novel\nancestral_number=$ancestral\n";
close O;

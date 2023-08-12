#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"g1=s","g2=s","a=s","o=s","pansyn=s","h|help");
if (!( defined $opts{g1} and defined $opts{g2} and defined $opts{a} and defined $opts{o} and defined $opts{pansyn})) {
		die "************************************************\n
	-g1	Full path to the [node*-Gain.OG.genes] file
	-g2	Full path to the [node*-Loss.OG.genes] file
	-a	Full path to the [Gain_Loss.emapper.annotations] file
	-o	Full path to the [outputDir_S33] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-g1	Full path to the [node*-Gain.OG.genes] file
	-g2	Full path to the [node*-Loss.OG.genes] file
	-a	Full path to the [Gain_Loss.emapper.annotations] file
	-o	Full path to the [outputDir_S33] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
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
####
if ($opts{pansyn} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{pansyn} =~ s/$slash$//;
}
####

my $name; my $chu; my $num; my $cha;

open I,"<$opts{g1}" or die("Could not open $opts{g1}\n"); 
my %h1=();my %h1_name=();

while (my $a=<I>) {
	chomp $a;
	my @it=split/\t/,$a;
	if (!exists $h1{$it[0]}) {
		$h1{$it[0]}=1;
	}
	else{
		$h1{$it[0]}=$h1{$it[0]}+1;
	}
	$h1_name{$it[1]}=$it[0];
}
close I;


open I,"<$opts{g2}" or die("Could not open $opts{g2}\n"); 
my %h2=();my %h2_name=();

while (my $a=<I>) {
	chomp $a;
	my @it=split/\t/,$a;
	if (!exists $h2{$it[0]}) {
		$h2{$it[0]}=1;
	}
	else{
		$h2{$it[0]}=$h2{$it[0]}+1;
	}
	$h2_name{$it[1]}=$it[0];
}
close I;

############count######
my %h1_count=();
foreach my $i (keys %h1_name) {
	$name=$h1_name{$i};
	$h1_count{$i}=$h1{$name};
}

my %h2_count=();
foreach my $i (keys %h2_name) {
	$name=$h2_name{$i};
	$h2_count{$i}=$h2{$name};
}

###choose annotation#####
open I,"<$opts{a}" or die("Could not open $opts{a}\n"); 
open O1,">$opts{o}/node_young.emapper.annotations" or die("Could not open $opts{o}/node_young.emapper.annotations\n"); 
open O2,">$opts{o}/node_old.emapper.annotations" or die("Could not open $opts{o}/node_old.emapper.annotations\n"); 

while (my $a=<I>) {
	chomp $a;
	if ($a!~/^#/) {
		my @it=split/\t/,$a;
		if (exists $h1_name{$it[0]}) {
			print O1 "$a\n";
		}
		if (exists $h2_name{$it[0]}) {
			print O2 "$a\n";
		}
	}
}
close I;close O1;close O2;
############annotation###
my @all = ("A", "B", "C","D", "E", "F", "G", "H","I", "J", "K", "L", "M","N", "O", "P", "Q", "S","T", "U", "V", "W","Y", "Z");

open I,"<$opts{o}/node_young.emapper.annotations" or die("Could not open $opts{o}/node_young.emapper.annotations\n"); 
my %count1=();

while (my $a=<I>) {
	chomp $a;
	if ($a!~/^#/) {
		my @it=split/\t/,$a;
		if ($it[6] !~/-/) {
			my @wo=split//,$it[6];
			my $num=@wo;
			foreach my $i (@wo) {
				if (!exists $count1{$i}) {
					my $chu=1/$h1_count{$it[0]};
					$count1{$i}=$chu/$num;
				}
				else{
					$chu=1/$h1_count{$it[0]};
					$count1{$i}=$count1{$i}+$chu/$num;
				}
			}
		}
	}
}
close I;

my %ha=();
open O,">$opts{o}/node_young_categories.counts" or die("Could not open $opts{o}/node_young_categories.counts\n"); 
foreach my $i (@all) {
	if (exists $count1{$i}) {
		print O "$i\t$count1{$i}\n";
		$ha{$i}=$count1{$i};
	}
	else{
		print O "$i\t0\n";
		$ha{$i}=0;
	}
}
close O;


##node2
open I,"<$opts{o}/node_old.emapper.annotations" or die("Could not open $opts{o}/node_old.emapper.annotations\n"); 
my %count2=();
#$na=0;
while (my $a=<I>) {
	chomp $a;
	if ($a!~/^#/) {
		my @it=split/\t/,$a;
		if ($it[6] !~/-/) {
			my @wo=split//,$it[6];
			$num=@wo;
			foreach my $i (@wo) {
				#if ($i eq "A") {
				#	print "$it[0]\n";
				#	$na=$na+1;
				#}
				if (!exists $count2{$i}) {
					$chu=1/$h2_count{$it[0]};
					$count2{$i}=$chu/$num;
				}
				else{
					$chu=1/$h2_count{$it[0]};
					$count2{$i}=$count2{$i}+$chu/$num;
				}
			}
		}
	}
}
close I;
my %hb=();
#print "$na\n";
open O,">$opts{o}/node_old_categories.counts" or die("Could not open $opts{o}/node_old_categories.counts\n"); 
foreach my $i (@all) {
	if (exists $count2{$i}) {
		print O "$i\t$count2{$i}\n";
		$hb{$i}=$count2{$i};
	}
	else{
		print O "$i\t0\n";
		$hb{$i}=0;
	}
}
close O;

open O,">$opts{o}/node_categories_variation.counts" or die("Could not open $opts{o}/node_categories_variation.counts\n"); 
print O "Categories\tValues\n";
foreach my $i (@all) {
	$cha=$ha{$i}-$hb{$i};
	print O "$i\t$cha\n";
}

system("Rscript $opts{pansyn}/Categories.r $opts{o} $opts{o}/node_categories_variation.counts");

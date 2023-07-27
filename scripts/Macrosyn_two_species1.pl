#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"pairs=s","m=s","n=s","gffm=s","gffn=s","o5=s","len=s","h|help");
if (!(defined $opts{pairs} and defined $opts{m} and defined $opts{n} and defined $opts{gffm} and defined $opts{gffn} and defined $opts{len} and defined $opts{o5})) {
		die "************************************************\n
	-pairs	Full path to the [*.pairs] file 
	-m	Enter the abbreviation for the name of the reference genome A (example: HSap)
	-n	Enter the abbreviation for the name of the genome B (example: MMus)
	-gffm	Full path to the simplified GFF file of the reference genome A (example: HSap_simplified.gff)
	-gffn	Full path to the simplified GFF file of the genome B (example: MMus_simplified.gff)
	-len	Full path to chromosomal length file of the reference genome A (example: HSap.len)
	-o5	Full path to the new [outputDir5] directory containing the output files
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-pairs	Full path to the [*.pairs] file 
	-m	Enter the abbreviation for the name of the reference genome A (example: HSap)
	-n	Enter the abbreviation for the name of the genome B (example: MMus)
	-gffm	Full path to the simplified GFF file of the reference genome A (example: HSap_simplified.gff)
	-gffn	Full path to the simplified GFF file of the genome B (example: MMus_simplified.gff)
	-len	Full path to chromosomal length file of the reference genome A (example: HSap.len)
	-o5	Full path to the new [outputDir5] directory containing the output files
	-h|-help Print this help page
		*************************************************\n";
}
my $slash;
####
if ($opts{o5} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o5} =~ s/$slash$//;
}
########

my $spename1=$opts{m};
my $spename2=$opts{n};
my $pathway_gffm="$opts{gffm}";
my $pathway_gffn="$opts{gffn}";
my $nv_cg=$spename1."_".$spename2;


if (-d "$opts{o5}/$nv_cg") {
	print "The $opts{o5}/$nv_cg directory already exists!\nThe old directory has been deleted.\n";
	system "rm -r $opts{o5}/$nv_cg";
}


system"mkdir $opts{o5}/$nv_cg";
print "The [$opts{o5}/$nv_cg] directory has been created!\n";

open IN,"<$opts{len}" or die("Could not open $opts{len}.\n"); 
open O,">$opts{o5}/$nv_cg/$opts{m}.block" or die("$opts{o5}/$nv_cg/$opts{m}.block.\n"); 
$n=0;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$n=$n+1;
	print O "$n\t$a[0]\t1\t$a[1]\n";
}
close IN;
close O;
$pathway_b="$opts{o5}/$nv_cg/$opts{m}.block";
get_matched_scaff($spename1,$spename2,$nv_cg,$pathway_gffm,$pathway_gffn);

sort_get_table($nv_cg);

#sub get_matched_scaff_2
sub get_matched_scaff{
    my $name=$_[2];
    my $spe1=$_[0];
    my $spe2=$_[1];
	my $gff_pathm=$_[3];
	my $gff_pathn=$_[4];
	my %blocks;
	open IN,"<$pathway_b" or die("Could not open $pathway_b.\n"); 
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    push @{$blocks{$a[1]}},[$a[0],$a[2],$a[3]];
    }
    close IN;
    my %gffPos=();
	open IN,"<$gff_pathm" or die("Could not open $gff_pathm.\n"); 
    while(my $bi=<IN>){
	    chomp $bi;
		if ($bi!~/^#/) {
			my @a=split /\t/,$bi;
			$gffPos{$a[1]}="$a[0]\t$a[2]\t$a[3]";
		}
    }
    close IN;
	open IN,"<$gff_pathn" or die("Could not open $gff_pathn.\n"); 
    while(my $bi2 =<IN>){
	    chomp $bi2;
		if ($bi2!~/^#/) {
			my @a=split /\t/, $bi2;
			$gffPos{$a[1]}="$a[0]\t$a[2]\t$a[3]";
		}
    }
    close IN;
    my %store=();
    my %nums=();
    open IN,"<","$opts{pairs}" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    next unless (exists $gffPos{$a[0]} && exists $gffPos{$a[1]});
	    my @cut1=split /\t/,$gffPos{$a[0]};
	    my @cut2=split /\t/,$gffPos{$a[1]};
	    if(exists $blocks{$cut1[0]}){
		    #print "match block:$_\n";
		    for(my $i=0;$i<@{$blocks{$cut1[0]}};$i++){
			    my $bid=${$blocks{$cut1[0]}}[$i][0];
			    my $sta=${$blocks{$cut1[0]}}[$i][1];
			    my $end=${$blocks{$cut1[0]}}[$i][2];
			    if($cut1[1]>=$sta && $cut1[2]<=$end){
				    #print "in block:$_\n";
				    $store{$cut2[0]}{$bid}++;
			    }
		    }
	    }
    }
    close IN;
    open O,"> $opts{o5}/$nv_cg/$name.scaff.match";
	foreach my $k1(sort keys%store){
	    foreach my $k2(keys%{$store{$k1}}){
		    print O "$k1\t$k2\t$store{$k1}{$k2}\n";
	    }
	}
	close O;
}
#sort_get_table.pl
sub sort_get_table{
	my $name=$_[0];my $num_cut=0;
    my %hash=();
    my %store2=();
    open IN,"<$opts{o5}/$nv_cg/$name.scaff.match" or die;
    while(my $a=<IN>){
	    chomp $a;
	    my @a=split/\t/,$a;
	    push @{$hash{$a[0]}},[$a[2],$a[1]];
	    $store2{$a[1]}{$a[0]}=$a[2];
    }
    open O,"> $opts{o5}/$nv_cg/$name.scaff.match.table";
    my %final=();
    foreach my $k1(sort keys%hash){
	    @{$hash{$k1}}=sort{$b->[0] <=> $a->[0]}@{$hash{$k1}};
	    my $blockid=${$hash{$k1}}[0][1];
	    my $score=${$hash{$k1}}[0][0];
	    #push @{$final{$blockid}},[$k1,$score];	
	    $final{$blockid}{$k1}=$score if ($score>=$num_cut);
    }
    my %list;
    print O "blockid";
    my @arr_sort;
    foreach my $k1(sort{$a<=>$b}keys%final){
	    foreach my $k2(sort{$final{$k1}{$b} <=> $final{$k1}{$a}} keys%{$final{$k1}}){
		    my $sca=$k2;
		    print O "\t$sca";
		    push @arr_sort,$sca;
	    }
    }
    print O "\n";
    foreach my $key1(sort{$a<=>$b}keys%store2){
	    print O "$key1";
	    foreach my $f(@arr_sort){
		    my $num=(exists $store2{$key1}{$f})?$store2{$key1}{$f}:0;
		    print O "\t$num";
	    }
	    print O "\n";
    }
    close O;close IN;
}


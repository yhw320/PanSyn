#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i1=s","m=s","n=s","gff=s","o1=s","h|help");
if (!(defined $opts{i1} and defined $opts{m} and defined $opts{n} and defined $opts{gff} and defined $opts{o1})) {
		die "************************************************\n
	-i1	Full path to the [inputDir1_S17] directory
	-m	The abbreviation for the name of the species that represents the ancestral genome (e.g. NVec)
	-n	The abbreviation for the name of the interested species (e.g. HSap)
	-gff	Full path to the gene coordinates file of the interested species (e.g. HSap_simplified.gff)
	-o1	Full path to the [outputDir1_S17] directory
	Optional:
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i1	Full path to the [inputDir1_S17] directory
	-m	The abbreviation for the name of the species that represents the ancestral genome (e.g. NVec)
	-n	The abbreviation for the name of the interested species (e.g. HSap)
	-gff	Full path to the gene coordinates file of the interested species (e.g. HSap_simplified.gff)
	-o1	Full path to the [outputDir1_S17] directory
	Optional:
	-h|-help Print this help page
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

my $spename1=$opts{m};
my $spename2=$opts{n};
my $gff1="Ancestors_database/".$spename1."/$spename1\_simplified.gff";
my $pathway_gffm="$opts{i1}/$gff1";
my $pathway_gffn="$opts{gff}";
my $pathway=$opts{o1};
my $block_name="Ancestors_database/".$spename1."/$spename1.block";
my $pathway_b="$opts{i1}/$block_name";
my $nv_cg=$spename1."_".$spename2;
if (-d "$pathway/$nv_cg") {
	print "The $pathway/$nv_cg directory already exists!\nThe old directory has been deleted.\n";
	system "rm -r $pathway/$nv_cg";
}
system"mkdir $pathway/$nv_cg";
print "The [$pathway/$nv_cg] directory has been created!\n";

system"cat $pathway/blastp/$spename1.blastp $pathway/blastp/$spename2.blastp > $pathway/$nv_cg/$nv_cg.blastp";
get_pairs_from_cluster($nv_cg,$spename1,$spename2);
#system"perl scripts/get_pairs_from_cluster.2.pl single_cluster/ancient.cluster2 $nv_cg/$nv_cg.blastp $spename1 $spename2 > $nv_cg/$nv_cg.pairs";
synteny_match_id($spename1,$spename2,$nv_cg);
#system"perl scripts/AJ-synteny-match-id2.pl -i ids_match/$spename1.ids.match -t ids_match/$spename2.ids.match -q $nv_cg/$nv_cg.pairs -o $nv_cg/$nv_cg.pairs2";
get_matched_scaff($spename1,$spename2,$nv_cg,$pathway_gffm,$pathway_gffn);
#system"perl scripts/get_matched_scaff_2.pl block/$spename1.block $nv_cg/$nv_cg.pairs2 formatted_gff/$spename1.gff gff/$spename2.gff> $nv_cg/$nv_cg.scaff.match";
sort_get_table($nv_cg);
#system"perl scripts/sort_get_table.pl $nv_cg/$nv_cg.scaff.match 0 > $nv_cg/$nv_cg.scaff.match.table";
#sub get_pairs_from_cluster.2.pl
sub get_pairs_from_cluster {
    my $one=$_[0];
    my $p1=$_[1];
    my $p2=$_[2];
    my %matchs=();
    open IN,"<","$pathway/$one/$one.blastp" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    next unless (($a[0]=~/^$p1\_/ && $a[1]=~/^$p2\_/) || ($a[0]=~/^$p2\_/ && $a[1]=~/^$p1\_/));
	    $matchs{$a[0]}{$a[1]}=$a[11];
    }
    close IN;
    my %list=();
	open IN,"<$pathway/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $pathway/ancient_gene_families.analysis/ancient_gene.families.\n"); 
    open O,"> $pathway/$one/$one.pairs";#NV_GG.pairs
    while(<IN>){
	    chomp;
	    s/,$//;
	    my @a=split /\t/;
	    my @b=split /,/,$a[1];
	    my @t1=();
	    my @t2=();
	    foreach my $f1(@b){
		    push @t1,$f1 if ($f1=~/^$p1\_/);
		    push @t2,$f1 if ($f1=~/^$p2\_/);
	    }
	    next unless (scalar(@t1)>0 && scalar(@t2)>0);
	    my %tmp;
	    foreach my $r1(@t1){
		    foreach my $r2(@t2){
			    if (exists $matchs{$r1}{$r2}){
				    push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r1}{$r2}];
			    }else{
				    push @{$tmp{'hh'}},[$r1,$r2,0]
			    }
			    if(exists $matchs{$r2}{$r1}){
				    push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r2}{$r1}];
			    }else{
				    push @{$tmp{'hh'}},[$r1,$r2,0]
			    }
		    }
	    } 
	    foreach my $k2(keys%tmp){
		    @{$tmp{$k2}}=sort{$b->[2]<=>$a->[2]}@{$tmp{$k2}};
		    for(my $i=0;$i<@{$tmp{$k2}};$i++){
			    my $g1=${$tmp{$k2}}[$i][0];
			    my $g2=${$tmp{$k2}}[$i][1];
			    if(!exists $list{$g1} && !exists $list{$g2}){
				    print O "$g1\t$g2\n";
				    $list{$g1}=1;
				    $list{$g2}=1;
			    }
		    }
	    }	
    }
    close IN;
	close O;
}
#sub synteny-match-id2.pl
sub synteny_match_id{
    my %h1;my %h2;my $spe_id11=();my $spe_id22=();
	my $spe1=$_[0];
    my $spe2=$_[1];
    my $name=$_[2];
    open I,"< $pathway/ids_match/$spe1.ids.match";#NV.id
    while (my $a=<I>) {
	    my $spe_id1=$a;
	    $spe_id11.=$spe_id1;
    }
    open T,"< $pathway/ids_match/$spe2.ids.match";#CG.id
    while (my $b=<T>) {
        my  $spe_id2=$b;
	    $spe_id22.=$spe_id2;
    }
    my $all_id=$spe_id11.$spe_id22;
    my @kk=split/\n/,$all_id;
    foreach $iii(@kk){
	    my @itms=split/\t/,$iii;
		$h1{$itms[1]}=$itms[0];
    }

    open Q,"< $pathway/$name/$name.pairs";#NV_GG.pairs
    open O,"> $pathway/$name/$name.pairs2";#NV_GG.pairs2
    while (my $a=<Q>) {
	    chomp $a;
	    my @itms=split/\t/,$a;
	    if (exists $h1{$itms[0]} and exists $h1{$itms[1]}) {
		    print O "$h1{$itms[0]}\t$h1{$itms[1]}\n";
	    }
	    else {print O "error\t$a\n";}
    }
    close I;close T;close Q;close O;
}
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
    open IN,"<","$pathway/$name/$name.pairs2" or die;
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
    open O,"> $pathway/$name/$name.scaff.match";
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
    open IN,"<$pathway/$name/$name.scaff.match" or die;
    while(my $a=<IN>){
	    chomp $a;
	    my @a=split/\t/,$a;
	    push @{$hash{$a[0]}},[$a[2],$a[1]];
	    $store2{$a[1]}{$a[0]}=$a[2];
    }
    open O,"> $pathway/$name/$name.scaff.match.table";
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

#system "rm -r $pathway/score";
#system "rm -r $pathway/cluster1";
#system "rm -r $pathway/table";
#system "rm -r $pathway/single_cluster";

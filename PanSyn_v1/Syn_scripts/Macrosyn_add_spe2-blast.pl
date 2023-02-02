#!/usr/bin/perl
use Getopt::Long;
my %opts;
use lib "/mnt/inspurfs/home/liyl/perl5/lib/perl5/SVG-2.86/lib";
use SVG;
GetOptions(\%opts,"m=s","n=s","i=s","c=s","color=s","g=s","gff=s","pep=s","l=s","s=s","f=s","o=s","t=s","evalue=s","e=s","p=s","name1=s","name2=s","h|help");
if (!(defined $opts{m} and defined $opts{n}  and defined $opts{color} and defined $opts{i} and defined $opts{g} and defined $opts{gff} and defined $opts{pep} and defined $opts{s} and defined $opts{f} and defined $opts{o})) {
		die "************************************************\n
	-m	Selected Ancestor'name by abbreviations (example: NVec)
	-n	Interested species'name byabbreviations (example: PAbe)
	-i	Full path to the [Ancestor_database] file
	-g	Full path to the interested species genome sequence file
	-gff	Full path to the interested species gff file
	-pep	Full path to the interested species protein sequence file
	-s	Full path to the [syn_scripts] file PanSyn provided
	-color	Full path to the Ancestor'chr color
	-f	Full path to the previous folder (same as [-o] in the previous script [Macrosyn1_1_blast.pl])
	-o  Full path to the new folder storing the output results
	Optional:
	-c	Specify whether to cluster the scaffolds (1 or 2)
		1: For species without chromosome-level assemblies
		2: For species with chromosome-level assemblies
	-t	Blast alignment threads (default:12)
	-evalue	Blast alignment evalue (default:1e-5)
	-e	Value indicating the tree-cutting threshold of 0.2-0.5 (default:0.5)
	-p	Pvalue for significance test (default: 0.05)
	-l	Full path where the cluster3.0 software command resides (example:/mnt/inspurfs/home/liyl/yuhw/software/cluster/bin/cluster)
	-name1	Ancestor 'name shown on the resulting graph (default:same as -m)
	-name2	Interested species'name shown on the resulting graph (default: same as -n)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
Options:
	-m	Selected Ancestor'name by abbreviations (example: NVec)
	-n	Interested species'name byabbreviations (example: PAbe)
	-i	Full path to the [Ancestor_database] file
	-g	Full path to the interested species genome sequence file
	-gff	Full path to the interested species gff file
	-pep	Full path to the interested species protein sequence file
	-s	Full path to the [syn_scripts] file PanSyn provided
	-f	Full path to the previous folder (same as [-o] in the previous script [Macrosyn1_1_blast.pl])
	-color	Full path to the Ancestor'chr color
	-o  Full path to the new folder storing the output results
	Optional:
	-c	Specify whether to cluster the scaffolds (1 or 2)
		1: For species without chromosome-level assemblies
		2: For species with chromosome-level assemblies
	-t	Blast alignment threads (default:12)
	-evalue	Blast alignment evalue (default:1e-5)
	-e	Value indicating the tree-cutting threshold of 0.2-0.5 (default:0.5)
	-p	Pvalue for significance test (default: 0.05)
	-l	Full path where the cluster3.0 software command resides (example:/mnt/inspurfs/home/liyl/yuhw/software/cluster/bin/cluster)
	-name1	Ancestor'name shown on the resulting graph (default: -m)
	-name2	Interested species'name shown on the resulting graph (default: -n)
	-h|-help Print this help page
		*************************************************\n";
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{evalue})) {
	$opts{evalue}=1e-5;
}
if (!(defined $opts{e})) {
	$opts{e}=0.5;
}
if (!(defined $opts{c})) {
	$opts{c}=2;
}
if (!(defined $opts{p})) {
	$opts{p}=0.05;
}
if (!(defined $opts{name1})) {
	$opts{name1}=$opts{m};
}
if (!(defined $opts{name2})) {
	$opts{name2}=$opts{n};
}
if ((defined $opts{l})) {
	my $cluster_pathway=$opts{l};
}
###############color file############################
%anchr_col=();
open I,"<$opts{color}";
while ($b=<I>){
	chomp $b;
	my @ittt=split/\t/,$b;
	$anchr_col{$ittt[0]}=$ittt[1];
}		
close I;
###########################################################33

system "mkdir $opts{o}/$opts{m}\_$opts{n}";
$fileExist = -e "$opts{pep}.psq";
if ( $fileExist ) {
}
else {
system "makeblastdb -in $opts{pep} -dbtype prot";
}

$fileExist = -e "$opts{f}/add_species/$opts{m}/all.peps.psq";
if ( $fileExist ) {
}
else {
system "makeblastdb -in $opts{f}/add_species/$opts{m}/all.peps -dbtype prot";
}


system "blastp -query $opts{pep} -db $opts{f}/add_species/$opts{m}/all.peps -outfmt 6 -out $opts{o}/$opts{m}\_$opts{n}/$opts{n}_all.blast -evalue $opts{e} -num_threads $opts{t}";
system "blastp -query $opts{f}/add_species/$opts{m}/all.peps -db $opts{pep} -outfmt 6 -out $opts{o}/$opts{m}\_$opts{n}/all_$opts{n}.blast -evalue $opts{e} -num_threads $opts{t}";


$opts{p}=sprintf("%.3f",$opts{p});

my %h1=();my $score;
open I,"< $opts{o}/$opts{m}\_$opts{n}/$opts{n}_all.blast";
while (my $a=<I>) {
	chomp $a;
	my @itms=split/\t/,$a;
	if (exists $h1{$itms[0]}) {
		if ($itms[11]>$score) {
			$h1{$itms[0]}=$itms[1];
			$score=$itms[11];
		}
	}
	else{
		$h1{$itms[0]}=$itms[1];
		$score=$itms[11];
	}
}
close I;
my %h2=();my $score2;
open I,"< $opts{o}/$opts{m}\_$opts{n}/all_$opts{n}.blast";
while (my $a=<I>) {
	chomp $a;
	my @itms=split/\t/,$a;
	if (exists $h2{$itms[0]}) {
		if ($itms[11]>$score2) {
			$h2{$itms[0]}=$itms[1];
			$score2=$itms[11];
		}
	}
	else{
		$h2{$itms[0]}=$itms[1];
		$score2=$itms[11];
	}
}
close I;

my %h_gff=();
open I,"<$opts{i}/$opts{m}/$opts{m}.gff" or die("Could not open $opts{i}/$opts{m}/$opts{m}.gff.\n"); 
while (my $a=<I>) {
	chomp $a;
	my @itms=split/\t/,$a;
	$h_gff{$itms[1]}=0;
}
close I;

my %h_f2=();
open I,"<$opts{f}/add_species/$opts{m}/all.pairs2" or die("Could not open $opts{f}/add_species/$opts{m}/all.pairs2.\n"); 
while (my $a=<I>) {
	chomp $a;
	my @itms=split/\t/,$a;
	$h_f2{$itms[1]}=$itms[0];
}
close I;

my %h_f1=();
foreach my $i1 (keys %h1) {
	my $zhi=$h1{$i1};
	if (exists $h2{$zhi} and $h2{$zhi} eq $i1) {
		$h_f1{$i1}=$zhi;
	}
}
my %h_f3=();
foreach my $i (keys %h_f1) {
	my $v= $h_f1{$i};
	if (exists $h_f2{$v}) {
		my $k=$h_f2{$v};
		$h_f3{$k}=$i;
	}
	if (exists $h_gff{$v}) {
		$h_f3{$v}=$i;
	}
}
open O,">$opts{o}/$opts{m}\_$opts{n}/$opts{m}\_$opts{n}.pairs2";
foreach my $i1 (keys %h_f3) {
	print O "$i1\t$h_f3{$i1}\n";
}
close O;

#上面是到生成pairs2的

my $spename1=$opts{m};
my $spename2=$opts{n};
my $pathway_gffm="$opts{i}/$opts{m}/$opts{m}.gff";
my $pathway_gffn="$opts{gff}";
my $pathway_b="$opts{i}/$opts{m}/$opts{m}.block";
get_matched_scaff($spename1,$spename2,$pathway_gffm,$pathway_gffn);
#system"perl scripts/get_matched_scaff_2.pl block/$spename1.block $nv_cg/$nv_cg.pairs2 formatted_gff/$spename1.gff gff/$spename2.gff> $nv_cg/$nv_cg.scaff.match";
sort_get_table();
#system"perl scripts/sort_get_table.pl $nv_cg/$nv_cg.scaff.match 0 > $nv_cg/$nv_cg.scaff.match.table";

sub get_matched_scaff{
    my $spe1=$_[0];
    my $spe2=$_[1];
	my $gff_pathm=$_[2];
	my $gff_pathn=$_[3];
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
    open IN,"<$opts{o}/$opts{m}\_$opts{n}/$opts{m}\_$opts{n}.pairs2" or die;
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
    open O,"> $opts{o}/$opts{m}\_$opts{n}/$opts{m}\_$opts{n}.scaff.match";
	foreach my $k1(sort keys%store){
	    foreach my $k2(keys%{$store{$k1}}){
		    print O "$k1\t$k2\t$store{$k1}{$k2}\n";
	    }
	}
	close O;
}
#sort_get_table.pl
sub sort_get_table{
	my $num_cut=0;
    my %hash=();
    my %store2=();
    open IN,"<$opts{o}/$opts{m}\_$opts{n}/$opts{m}\_$opts{n}.scaff.match" or die;
    while(my $a=<IN>){
	    chomp $a;
	    my @a=split/\t/,$a;
	    push @{$hash{$a[0]}},[$a[2],$a[1]];
	    $store2{$a[1]}{$a[0]}=$a[2];
    }
    open O,"> $opts{o}/$opts{m}\_$opts{n}/$opts{m}\_$opts{n}.scaff.match.table";
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

###下面开始画图，进行Macrosyn3.pl的步骤

my $pathway_block=$pathway_b;
my $pathway_pepm="$opts{i}/$spename1/$spename1.pep";
my $pathway_pepn="$opts{pep}";

my $cutoff=$opts{e};
my $chromosome=$opts{c};
my $pathway_genome=$opts{g};
my $pathway=$opts{o};
my $nv_cg=$spename1."_".$spename2;
my $c_value;
my $nv_cgn=$nv_cg."_scaff_c";
my $nv_cg3=$nv_cg."_scaff_c1";
my $info=$spename2."_chr2scaf";
my $scr=$opts{s};	$ci_num2="P_value<=".$opts{p};
my $circle_num=0;my $n2;my $i2;my $final;


if ($chromosome==1) {
    system"$cluster_pathway -f $pathway/$nv_cg/$nv_cg.scaff.match.table -m a -e 1"; 
    testDB_scriptN_scaffold_clustering($nv_cg,$cutoff,$nv_cg3,$info);
    #system"perl scripts/testDB_scriptN_scaffold_clustering.pl -c $nv_cg/$nv_cg.scaff.match.cdt -a $nv_cg/$nv_cg.scaff.match.atr -p $cutoff -t $nv_cg/$nv_cg.scaff.match.table -n c1 -o $nv_cg/$nv_cg3 -e $nv_cg/$info.info";
    system"$cluster_pathway -f $pathway/$nv_cg/$nv_cg3 -m a -e 1"; 
    system"awk 'NR==1{print \$4}' $pathway/$nv_cg/$nv_cg3.atr > $pathway/$nv_cg/$nv_cg.scaff.match.atr_num";
    open NUM,"<$pathway/$nv_cg/$nv_cg.scaff.match.atr_num";
    my $num=<NUM>;
    $n2=1;
    while($num > $cutoff){
	    $i2=$nv_cg."_scaff_c".$n2;
	    $n2=$n2+1;
	    $c_value=$nv_cgn.$n2;
        testDB_scriptN_scaffold_clusteringb($nv_cg,$i2,$cutoff,$n2,$c_value,$info);
  	    #system"perl scripts/testDB_scriptN_scaffold_clustering.pl -c $nv_cg/$i2.cdt -a $nv_cg/$i2.atr -p $cutoff -t $nv_cg/$i2 -n c$n2 -o $nv_cg/$c_value -e $nv_cg/$info.info";
	    system"$cluster_pathway -f $pathway/$nv_cg/$c_value -m a -e 1"; 
        system"awk 'NR==1{print \$4}' $pathway/$nv_cg/$c_value.atr > $pathway/$nv_cg/$nv_cg.scaff.match.atr_num";
	    open NUM,"<$pathway/$nv_cg/$nv_cg.scaff.match.atr_num";
	    $num=<NUM>;
        $circle_num=$circle_num+1;
		if ($circle_num>2) {
			system"rm $pathway/$nv_cg/$nv_cg.scaff.match.cdt";
			system"rm $pathway/$nv_cg/$nv_cg.scaff.match.atr";
			system"rm $pathway/$nv_cg/$nv_cg\_scaff\_c*";
			system"rm $pathway/$nv_cg/$info.info";
			system"rm $pathway/$nv_cg/$nv_cg.scaff.match.atr_num";
			die "Please enlarge the cutoff value.\n";
		}
	}
    $final=$nv_cg."_scaff_c".$n2;
    system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/$final > $pathway/$nv_cg/$final-reverse";
    scaff_c_formatted($nv_cg,$final);
	#system"perl scripts/scaff_c_formatted.pl -i $pathway/$nv_cg/$final-reverse -o $pathway/$nv_cg/$final-and";
    #system"sed -i '1d' $pathway/$nv_cg/$final-and";
    system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/$final-and > $pathway/$nv_cg/$final";
	system "rm $pathway/$nv_cg/$final-and";
	system "rm $pathway/$nv_cg/$final-reverse";
	len ($pathway_genome);
	#system"perl scripts/len.pl $pathway_genome/$spename2.fa";
    testDB_scriptN_scaffold_clusteringn($pathway_genome,$spename2,$nv_cg,$info);
	#system"perl scripts/testDB_scriptN_scaffold_clusteringn.pl -i $pathway_genome/$spename2.fa.len -c $pathway/$nv_cg/$info.info -o $pathway/$nv_cg/$spename2.chr.len";
    plot_dot_hd1($pathway_block,$spename1,$nv_cg,$spename2,$final);
	plot_dot_hd12($pathway_block,$spename1,$nv_cg,$spename2,$final);
	#system"perl scripts/plot_dot_hd.pl $pathway_block $pathway/$nv_cg/$spename2.chr.len $pathway/$nv_cg/$final $pathway/$nv_cg/$nv_cg";
}
if ($chromosome==2) {
	len ($pathway_genome);
	#system"perl scripts/len.pl $pathway_genome/$spename2.fa";
	system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/$nv_cg.scaff.match.table > $pathway/$nv_cg/$nv_cg.scaff.match.table-reverse";
	system "awk '{print \$1\"\t\"\$1}'  $pathway/$nv_cg/$nv_cg.scaff.match.table-reverse > $pathway/$nv_cg/$spename2\_chr2scaf.info";
	system"sed -i '1d' $pathway/$nv_cg/$spename2\_chr2scaf.info";
	$nnnn=$nv_cg.".scaff.match.table";
    scaff_c_formatted($nv_cg,$nnnn);
	#system"perl scripts/scaff_c_formatted.pl -i $pathway/$nv_cg/$nv_cg.scaff.match.table-reverse -o $pathway/$nv_cg/$nv_cg.scaff.match.table-and";
    #system"sed -i '1d' $pathway/$nv_cg/$nv_cg.scaff.match.table-and";
    system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/$nv_cg.scaff.match.table-and > $pathway/$nv_cg/$nv_cg.scaff.match.table";
	system "rm $pathway/$nv_cg/$nv_cg.scaff.match.table-and";
	system "rm $pathway/$nv_cg/$nv_cg.scaff.match.table-reverse";
    plot_dot_hd2($pathway_block,$spename1,$nv_cg,$spename2,$pathway_genome);
    plot_dot_hd22($pathway_block,$spename1,$nv_cg,$spename2,$pathway_genome,$scr);
	#system"perl scripts/plot_dot_hd.pl $pathway_block/$spename1.block $pathway_genome/$spename2.fa.len $pathway/$nv_cg/$nv_cg.scaff.match.table $pathway/$nv_cg/$nv_cg";
}

#testDB_scriptN_scaffold_clustering.pl
sub testDB_scriptN_scaffold_clustering{
	my $name=$_[0];my $cut=$_[1];my $name2=$_[2];my $inf=$_[3];
    my (@sca,@aid,%h1,%h2,%h3,%h4);
    my $scaN;
    open C,"< $pathway/$nv_cg/$name.scaff.match.cdt"; #cdt file
    while (my $a=<C>) {
	    chomp $a;
	    if ($a=~/^blockid/) {
		    my @itms=split/\t/,$a;
		    $scaN=@itms-3;
		    foreach my $i (3..$#itms) {
			    push @sca,$itms[$i];
		    }
	    }
	    elsif ($a=~/^AID/) {
		    my @itms=split/\t/,$a;
		    foreach my $j (3..$#itms) {
			    push @aid,$itms[$j];
		    }
	    }
    }
    foreach my $k (0..$scaN-1) {
	    $h1{$aid[$k]}=$sca[$k];
    }
    open A,"< $pathway/$nv_cg/$name.scaff.match.atr";
        while (my $a=<A>) {
	    chomp $a;
	    my @itms=split/\t/,$a;
	    if ($itms[3]>=$cutoff) {
		    my $mem1=exists $h2{$itms[1]}?$h2{$itms[1]}:$itms[1];
		    my $mem2=exists $h2{$itms[2]}?$h2{$itms[2]}:$itms[2];
		    $h2{$itms[0]}="$mem1\t$mem2";
		    if (exists $h2{$itms[1]}) {
			    delete $h2{$itms[1]};
		    }
		    if (exists $h2{$itms[2]}) {
			    delete $h2{$itms[2]};
		    }
	    }
	    else {
		    if ($itms[1]=~/ARRY/) {
			    $h2{$itms[1]}=$itms[1];
		    }
		    if ($itms[2]=~/ARRY/) {
			    $h2{$itms[2]}=$itms[2];
		    }
	    }
    }
    my %sin;
    foreach my $k (sort keys %h2) {
	    my @tmp=split/\t/,$h2{$k};
	    foreach my $i (@tmp) {
		    $h3{$h1{$i}}=$k;
		    push @{$sin{$k}},$h1{$i};
	    }
    }
    open E,">> $pathway/$nv_cg/$inf.info";
    my $num=1;
    foreach my $a (sort keys %sin) {
	    my $chr="c1"."-chr"."$num";
	    my $LINE=join "\t",@{$sin{$a}};
	    print E "$chr\t$LINE\n";
	    $num++;
    }
    open T,"< $pathway/$nv_cg/$name.scaff.match.table";
    my @SCA;
    while (my $a=<T>) {
	    chomp $a;
	    if ($a=~/blockid/) {
		    my @itms=split/\t/,$a;
		    shift @itms;
		    @SCA=@itms;
	     }
	    else {
		    my @itms2=split/\t/,$a;
		    my $alg=shift @itms2;
		    foreach my $i (0..$#itms2) {
			    my $node=$h3{$SCA[$i]};
			    $h4{$alg}{$node}+=$itms2[$i];
		    }
	    }
    }
    open O,"> $pathway/$nv_cg/$name2";
    my @CHR;
    my @NODE=sort keys %{$h4{"1"}};
    foreach my $n (0..$#NODE) {
	    my $N=$n+1;
	    my $chr="c1-chr".$N;
	    push @CHR,$chr;
    }
    my $title=join "\t",@CHR;
    print O "blockid\t$title\n";
    foreach my $x (sort {$a<=>$b} keys %h4) {
	    my $line;
	    foreach my $y (sort keys %{$h4{$x}}) {
		$line.="$h4{$x}{$y}\t";
	    }
	    chop $line;
	    print O "$x\t$line\n";
    }
    close C;close A;close E;close T;close O;
}
sub testDB_scriptN_scaffold_clusteringb{
	my $name=$_[0];my $ii=$_[1];my $in=$_[5];my $cv=$_[4];my $zhi=$_[2];my $hh=$_[3];
    my (@sca,@aid,%h1,%h2,%h3,%h4);
    my $scaN;
    open C,"< $pathway/$nv_cg/$ii.cdt"; #cdt file
    while (my $a=<C>) {
	    chomp $a;
	    if ($a=~/^blockid/) {
		    my @itms=split/\t/,$a;
		    $scaN=@itms-3;
		    foreach my $i (3..$#itms) {
			push @sca,$itms[$i];
		    }
	    }
	    elsif ($a=~/^AID/) {
		    my @itms=split/\t/,$a;
		    foreach my $j (3..$#itms) {
			    push @aid,$itms[$j];
		    }
	    }
    }

    foreach my $k (0..$scaN-1) {
	    $h1{$aid[$k]}=$sca[$k];
    }
    open A,"< $pathway/$nv_cg/$ii.atr";
    while (my $a=<A>) {
	    chomp $a;
	    my @itms=split/\t/,$a;
	    if ($itms[3]>=$zhi) {
		    my $mem1=exists $h2{$itms[1]}?$h2{$itms[1]}:$itms[1];
		    my $mem2=exists $h2{$itms[2]}?$h2{$itms[2]}:$itms[2];
		    $h2{$itms[0]}="$mem1\t$mem2";
		    if (exists $h2{$itms[1]}) {
			    delete $h2{$itms[1]};
		    }
		    if (exists $h2{$itms[2]}) {
			    delete $h2{$itms[2]};
		    }
	    }
	    else {
		    if ($itms[1]=~/ARRY/) {
			    $h2{$itms[1]}=$itms[1];
		    }
		    if ($itms[2]=~/ARRY/) {
			    $h2{$itms[2]}=$itms[2];
		    }
	    }
    }
    my %sin;
    foreach my $k (sort keys %h2) {
	    my @tmp=split/\t/,$h2{$k};
	    foreach my $i (@tmp) {
		    $h3{$h1{$i}}=$k;
		    push @{$sin{$k}},$h1{$i};
	    }
    }
    open E,">> $pathway/$nv_cg/$in.info";
    my $num=1;
    foreach my $a (sort keys %sin) {
	    my $chr="c".$hh."-chr"."$num";
	    my $LINE=join "\t",@{$sin{$a}};
	    print E "$chr\t$LINE\n";
	    $num++;
    }

    open T,"< $pathway/$nv_cg/$ii";
    my @SCA;
    while (my $a=<T>) {
	    chomp $a;
	    if ($a=~/blockid/) {
		    my @itms=split/\t/,$a;
		    shift @itms;
		    @SCA=@itms;
	    }
	    else {
		    my @itms2=split/\t/,$a;
		    my $alg=shift @itms2;
		    foreach my $i (0..$#itms2) {
			    my $node=$h3{$SCA[$i]};
			    $h4{$alg}{$node}+=$itms2[$i];
		    }
	    }
    }
    open O,"> $pathway/$nv_cg/$cv";
    my @CHR;
    my @NODE=sort keys %{$h4{"1"}};
    foreach my $n (0..$#NODE) {
	    my $N=$n+1;
	    my $chr="c".$hh."-chr".$N;
	    push @CHR,$chr;
    }
    my $title=join "\t",@CHR;
    print O "blockid\t$title\n";

    foreach my $x (sort {$a<=>$b} keys %h4) {
	    my $line;
	    foreach my $y (sort keys %{$h4{$x}}) {
		    $line.="$h4{$x}{$y}\t";
	    }
	    chop $line;
	    print O "$x\t$line\n";
    }
	close C;close A;close E;close T;close O;

}

#scaff_c_formatted.pl
sub scaff_c_formatted{
	my $na=$_[0]; my $f=$_[1];my $n_hang=0;
    my $all_sum=0;my $all_max=0;my $max=0;my $num_max=0;$nn=0.001;my @items;my $num2;my $n;my %h11;my @items2;my $name_one;my %h21;my $fin;my $hang;my $fir;
    open I,"<$pathway/$na/$f-reverse";
    while (my $a=<I>){
		$n_hang=$n_hang+1;
	    chomp $a;
	    if ($a=~/^blockid/) {
            @items=split/\t/,$a;
		    $fir=0;
		    $h11{$fir}=$items[0];
	    }
	    if ($a!~/^blockid/) {
            @items=split/\t/,$a;
		    $num2=@items;
            for($n=1;$n<$num2;$n=$n+1){
				$all_sum=$all_sum+$items[$n];
                if ($items[$n]>$max) {
				    $max=$items[$n];
				    $num_max=$n;
					$haxi{$n_hang}=$max;
                }
			    else{
				    $max=$max;
			    }
		    }
		    if (exists $h11{$num_max}) {
                $num_max=$num_max+$nn;
			    $nn=$nn+0.001;
		    }
	        $h11{$num_max}=$items[0];
	        $max=0;
		}
    }
    close I;
    open I,"$pathway/$na/$f-reverse";
    while (my $b=<I>){
	    chomp $b;
        @items2=split/\t/,$b;
	    $name_one=$items2[0];
	    $h21{$name_one}=$b;
    }
    close I;
    open O,">$pathway/$na/$f-and";
    foreach my $i(sort {$a <=> $b} keys %h11){
	    $fin=$h11{$i};
        $hang=$h21{$fin};
        print O "$hang\n";
    }
    close O;
	foreach $kkey (keys %haxi) {
		$all_max=$all_max+$haxi{$kkey};
	}
	#print "$all_max\n";
	#print "$all_sum\t$all_max\n";
	$papa=$all_max/$all_sum;
	$aaaa=sprintf("%.3f",$papa);
	$ci_num="Conservation_Index(CI)=".$aaaa;
	#$ci_num="Conservation_Index(CI)=".$aaaa." "."($all_max/$all_sum)";
	$ci_nump="Conservation_Index(CI)=".$aaaa;
}
#len.pl
sub len{
	my $pa=$_[0];
	open I,"<$pa" or die("Could not open $pa.\n"); 
    while(<I>){
        chomp;
        if (/>(\S+)/) { $id=$1; }
        else { $h{$id} .= $_; } 
    }
    close I;
    open O,">$pa.len";
    foreach $k (sort keys %h){
	    $len=length($h{$k});
	    print O "$k\t$len\n";
    }
    close O;
}
#testDB_scriptN_scaffold_clusteringn.pl
sub testDB_scriptN_scaffold_clusteringn{
    my $pa_ge=$_[0];my $spe2=$_[1];my $name=$_[2];my $in=$_[3];
    my (%H1,%H2);
    open I,"< $pa_ge.len";#scaffold length file
    while (my $a=<I>) {
	    chomp $a;
	    my @itms=split/\t/,$a;
	    $H1{$itms[0]}=$itms[1];
    }
    open C,"< $pathway/$nv_cg/$in.info";
    while (my $a=<C>) {
	    chomp $a;
	    if ($a=~/^c1/) {
		    my @itms=split/\t/,$a;
		    my $cid=shift @itms;
		    push @{$H2{$cid}},@itms;
	    }
	    else {
		    my @itms2=split/\t/,$a;
		    my $cid=shift @itms2;
		    foreach my $i (0..$#itms2) {
			    push @{$H2{$cid}},@{$H2{$itms2[$i]}};
			    delete $H2{$itms2[$i]};
		    }
	    }	
    }
    open O,"> $pathway/$nv_cg/$spe2.chr.len";
    foreach my $k (sort keys %H2) {
	    foreach my $m (@{$H2{$k}}) {
		    $H3{$k}+=$H1{$m};
	    }
	    print O "$k\t$H3{$k}\n";
    }
	close I;close C;close O;
}
#plot
sub plot_dot_hd1{
    my $path_block=$_[0];my $spe1=$_[1];my $name=$_[2];my $spe2=$_[3];my $fi=$_[4];
    my $svg;
    my ($left,$right)=(100,100);
    my ($top,$bottom)=(100,100);
    my ($width,$height)=(800,800);
    $svg=SVG->new('width',$left+$width+$right,'height',$top+$height+$bottom);
    #my $font = FontSize->new();
    my %sca_arr=();
    my @sca_sorts=();
    open IN,"<","$pathway/$nv_cg/$fi" or die;
    while(<IN>){
	    chomp;
	    if(/^blockid/){
		    my @a=split /\t/;		
		    foreach my $fa(@a){
			    $sca_arr{$fa}=1;
			    push @sca_sorts,$fa if ($fa ne 'blockid');
		    }
	    }
    }
    close IN;
    my %hsa_len=();
    open IN,"<","$path_block" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    $hsa_len{$a[0]}+=($a[3]-$a[2]+1);
    }
    close IN;

    my %lg_len=();
    open IN,"<","$pathway/$nv_cg/$spe2.chr.len" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    $lg_len{$a[0]}=$a[1];	
    }
    close IN;
    my $hsa_sum_lens=0;
    foreach my $k1(keys%hsa_len){
	    $hsa_sum_lens+=$hsa_len{$k1};
    }
    my $lg_sum_lens=0;
    foreach my $k2(keys%lg_len){
	    $lg_sum_lens+=$lg_len{$k2} if (exists $sca_arr{$k2});
    }

    my $ybin=$height/$hsa_sum_lens;
    my $xbin=$width/$lg_sum_lens;
    #print "xbin:$xbin\n";

    my $plot_y_lens=$top;
    my $font_family = "Arial";
    my $font_size = 12;

    $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
    foreach my $q1(sort {$a<=>$b}keys %hsa_len){
	    #print $q1."\n";
	    my $plot_start=$plot_y_lens;
        my $binn=$hsa_len{$q1}*$ybin;
	    $plot_y_lens+=$hsa_len{$q1}*$ybin;
	    my $name_height = 20; #$font->stringHeight($font_size);
            my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q1);
	    $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
		$q1=~s/c[0-9]-chr//;
		$svg->text('x', $left - 10 - $bac_name_width,'y',$plot_start+($binn-$name_height)/2, '-cdata',$q1,'font-family',$font_family,'font-size',$font_size);
    }

    my $plot_x_lens=$left;
    $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
    foreach my $q2(@sca_sorts){
	    my $plot_start=$plot_x_lens;
	    my $binn=$lg_len{$q2}*$xbin;
	    $plot_x_lens+=$lg_len{$q2}*$xbin;
	    my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q2);
	    $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
	    $svg->text('x', $plot_start+($binn-$bac_name_width)/2,'y',$top+$height+20, '-cdata',$q2,'font-family',$font_family,'font-size',$font_size);
    }
    my $sum_y_start=$top;
    open IN,"<","$pathway/$nv_cg/$fi" or die;
    while(<IN>){
	    chomp;
	    next if /^blockid/;
	    my @a=split /\t/;
	    my $sum_x_start=$left;
	    my $sum_y_len=$ybin*$hsa_len{$a[0]};
	    foreach my $n1(1..$#a){
		    my $sum_x_len=$xbin*$lg_len{$sca_sorts[$n1-1]};
		    #print "***********\n$sum_x_len\n";
		    foreach(1..$a[$n1]){
			    my $randx=$sum_x_start+$sum_x_len*rand(1);
			    my $randy=$sum_y_start+$sum_y_len*rand(1);
			    $svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#3399CC');
		    }
		    $sum_x_start+=$sum_x_len;
	    }
	    $sum_y_start+=$sum_y_len;
    }
    close IN;
	$svg->text('x', 350,'y',85, '-cdata',$ci_num,'font-family',$font_family,'font-size',25);
	$svg->text('x', 450,'y',950, '-cdata',$opts{name2},'font-family',$font_family,'font-size',25);
	$svg->text('x',43,'y',500, '-cdata',$opts{name1},'font-family',$font_family,'font-size',25,'transform','rotate(270,43 500)');
    open OUT,">","$pathway/$nv_cg/$name.dot.svg" or die;
    print OUT $svg->xmlify;
	close OUT;
}
sub plot_dot_hd2{
    my $path_block=$_[0];my $spe1=$_[1];my $name=$_[2];my $spe2=$_[3];my $path_geno=$_[4];
    my $svg;
    my ($left,$right)=(100,100);
    my ($top,$bottom)=(100,100);
    my ($width,$height)=(800,800);
    $svg=SVG->new('width',$left+$width+$right,'height',$top+$height+$bottom);
    #my $font = FontSize->new();
    my %sca_arr=();
    my @sca_sorts=();
    open IN,"<","$pathway/$nv_cg/$name.scaff.match.table" or die;
        while(<IN>){
	        chomp;
	        if(/^blockid/){
		        my @a=split /\t/;		
		        foreach my $fa(@a){
			        $sca_arr{$fa}=1;
			        push @sca_sorts,$fa if ($fa ne 'blockid');
		        }
	         }
        }
        close IN;
        my %hsa_len=();
        open IN,"<","$path_block" or die;
        while(<IN>){
	        chomp;
	        my @a=split /\t/;
	        $hsa_len{$a[0]}+=($a[3]-$a[2]+1);
        }
        close IN;
        my %lg_len=();
        open IN,"<","$path_geno.len" or die;
        while(<IN>){
	        chomp;
	        my @a=split /\t/;
	        $lg_len{$a[0]}=$a[1];	
        }
        close IN;
        my $hsa_sum_lens=0;
        foreach my $k1(keys%hsa_len){
	        $hsa_sum_lens+=$hsa_len{$k1};
        }
        my $lg_sum_lens=0;
        foreach my $k2(keys%lg_len){
	        $lg_sum_lens+=$lg_len{$k2} if (exists $sca_arr{$k2});
        }

        my $ybin=$height/$hsa_sum_lens;
        my $xbin=$width/$lg_sum_lens;
        #print "xbin:$xbin\n";

        my $plot_y_lens=$top;
        my $font_family = "Arial";
        my $font_size = 12;

        $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
        foreach my $q1(sort {$a<=>$b}keys %hsa_len){
	        #print $q1."\n";
	        my $plot_start=$plot_y_lens;
            my $binn=$hsa_len{$q1}*$ybin;
	        $plot_y_lens+=$hsa_len{$q1}*$ybin;
	        my $name_height = 20; #$font->stringHeight($font_size);
                my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q1);
	        $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
	        $svg->text('x', $left - 10 - $bac_name_width,'y',$plot_start+($binn-$name_height)/2, '-cdata',$q1,'font-family',$font_family,'font-size',$font_size);
        }

        my $plot_x_lens=$left;
        $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
        foreach my $q2(@sca_sorts){
	        my $plot_start=$plot_x_lens;
	        my $binn=$lg_len{$q2}*$xbin;
	        $plot_x_lens+=$lg_len{$q2}*$xbin;
	        my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q2);
	        $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
	        $svg->text('x', $plot_start+($binn-$bac_name_width)/2,'y',$top+$height+20, '-cdata',$q2,'font-family',$font_family,'font-size',$font_size);
        }
        my $sum_y_start=$top;
        open IN,"<","$pathway/$nv_cg/$name.scaff.match.table" or die;
        while(<IN>){
	        chomp;
	        next if /^blockid/;
	        my @a=split /\t/;
	        my $sum_x_start=$left;
	        my $sum_y_len=$ybin*$hsa_len{$a[0]};
	        foreach my $n1(1..$#a){
		        my $sum_x_len=$xbin*$lg_len{$sca_sorts[$n1-1]};
		        #print "***********\n$sum_x_len\n";
		        foreach(1..$a[$n1]){
			        my $randx=$sum_x_start+$sum_x_len*rand(1);
			        my $randy=$sum_y_start+$sum_y_len*rand(1);
			        $svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#3399CC');
		        }
		        $sum_x_start+=$sum_x_len;
	        }
	        $sum_y_start+=$sum_y_len;
        }
        close IN;
		$svg->text('x', 350,'y',85, '-cdata',$ci_num,'font-family',$font_family,'font-size',25);
        $svg->text('x', 450,'y',950, '-cdata',$opts{name2},'font-family',$font_family,'font-size',25);
		$svg->text('x',43,'y',500, '-cdata',$opts{name1},'font-family',$font_family,'font-size',25,'transform','rotate(270,43 500)');
		open OUT,">","$pathway/$nv_cg/$name.dot.svg" or die;
        print OUT $svg->xmlify;
		close OUT;
}

##plotPvalue
sub plot_dot_hd12{
	my $path_block=$_[0];my $spe1=$_[1];my $name=$_[2];my $spe2=$_[3];my $fi=$_[4];
	$sum_before=0;
	open I,"<$pathway/$nv_cg/$fi";
	open O,">$pathway/$nv_cg/result";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			print O "$a\tsumrow\n";
		}
		else{
			@items=split/\t/,$a;
			foreach $i (1..$#items) {
				$sum=$items[$i]+$sum_before;
				$sum_before=$sum;
			}
			print O "$a\t$sum_before\n";
		}
		$sum_before=0;
	}
	close I;
	close O;
	system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/result > $pathway/$nv_cg/result-reverse";
	$sum_before=0;
	open I,"< $pathway/$nv_cg/result-reverse";
	open O,">$pathway/$nv_cg/sum-result";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			print O "$a\tsumcol\n";
		}
		else{
			@items=split/\t/,$a;
			foreach $i (1..$#items) {
				$sum=$items[$i]+$sum_before;
				$sum_before=$sum;
			}
			print O "$a\t$sum_before\n";
		}
		$sum_before=0;
	}
	close I;
	close O;
	system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/sum-result > $pathway/$nv_cg/sum-table";
	system "rm $pathway/$nv_cg/result";
	system "rm $pathway/$nv_cg/sum-result";
	system "rm $pathway/$nv_cg/result-reverse";
	system "Rscript $opts{s}/Fish.R $pathway/$name $pathway/$nv_cg/sum-table";
	open I,"<$pathway/$nv_cg/$fi";
	open O,">$pathway/$nv_cg/first_line";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			$a=~s/blockid//g;
			$a=~s/^\t//g;
			print O "$a";
		}
	}
	print O "\n";
	close I;
	close O;

	open M,"<$pathway/$nv_cg/Pvalue-result.txt";
	open O1,">$pathway/$nv_cg/Pvalue-result1.txt";
	while ($a=<M>){
		chomp $a;
		$a=~s/ //g;
		print O1 "$a\n";
	}
	close M;
	close O1;

	system "cat $pathway/$nv_cg/first_line $pathway/$nv_cg/Pvalue-result1.txt > $pathway/$nv_cg/Pvalue-result2.txt";
	system "awk '{print \$1}' $pathway/$nv_cg/$fi >$pathway/$nv_cg/line1";
	system "paste $pathway/$nv_cg/line1 $pathway/$nv_cg/Pvalue-result2.txt > $pathway/$nv_cg/Pvalue.table";
	system "rm $pathway/$nv_cg/Pvalue-result1.txt";
	system "rm $pathway/$nv_cg/first_line";
	system "rm $pathway/$nv_cg/Pvalue-result2.txt";
	system "rm $pathway/$nv_cg/line1";

    my $svg;
    my ($left,$right)=(100,100);
    my ($top,$bottom)=(100,100);
    my ($width,$height)=(800,800);
    $svg=SVG->new('width',$left+$width+$right,'height',$top+$height+$bottom);
    #my $font = FontSize->new();
    my %sca_arr=();
    my @sca_sorts=();
    open IN,"<","$pathway/$nv_cg/$fi" or die;
    while(<IN>){
	    chomp;
	    if(/^blockid/){
		    my @a=split /\t/;		
		    foreach my $fa(@a){
			    $sca_arr{$fa}=1;
			    push @sca_sorts,$fa if ($fa ne 'blockid');
		    }
	    }
    }
    close IN;
    my %hsa_len=();
    open IN,"<","$path_block" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    $hsa_len{$a[0]}+=($a[3]-$a[2]+1);
    }
    close IN;

    my %lg_len=();
    open IN,"<","$pathway/$nv_cg/$spe2.chr.len" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    $lg_len{$a[0]}=$a[1];	
    }
    close IN;
    my $hsa_sum_lens=0;
    foreach my $k1(keys%hsa_len){
	    $hsa_sum_lens+=$hsa_len{$k1};
    }
    my $lg_sum_lens=0;
    foreach my $k2(keys%lg_len){
	    $lg_sum_lens+=$lg_len{$k2} if (exists $sca_arr{$k2});
    }

    my $ybin=$height/$hsa_sum_lens;
    my $xbin=$width/$lg_sum_lens;
    #print "xbin:$xbin\n";

    my $plot_y_lens=$top;
    my $font_family = "Arial";
    my $font_size = 12;
	my $plot_y_lens_before=0;
	my $num2=1;
	my %h2=();my %h1=();
    $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
    foreach my $q1(sort {$a<=>$b}keys %hsa_len){
	    #print $q1."\n";
	    my $plot_start=$plot_y_lens;
        my $binn=$hsa_len{$q1}*$ybin;
	    $plot_y_lens+=$hsa_len{$q1}*$ybin;
	    my $name_height = 20; #$font->stringHeight($font_size);
            my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q1);
	    $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
		$h2{$num2}=$plot_y_lens_before."\t".$plot_y_lens;
		$plot_y_lens_before=$plot_y_lens;
		$num2=$num2+1;
		$q1=~s/c[0-9]-chr//;
		$svg->text('x', $left - 10 - $bac_name_width,'y',$plot_start+($binn-$name_height)/2, '-cdata',$q1,'font-family',$font_family,'font-size',$font_size);
    }
	my $num1=1;my $plot_x_lens_before=0;
    my $plot_x_lens=$left;
    $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
    foreach my $q2(@sca_sorts){
	    my $plot_start=$plot_x_lens;
	    my $binn=$lg_len{$q2}*$xbin;
	    $plot_x_lens+=$lg_len{$q2}*$xbin;
	    my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q2);
	    $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
		$h1{$num1}=$plot_x_lens_before."\t".$plot_x_lens;
		$plot_x_lens_before=$plot_x_lens;
		$num1=$num1+1;
		$svg->text('x', $plot_start+($binn-$bac_name_width)/2,'y',$top+$height+20, '-cdata',$q2,'font-family',$font_family,'font-size',$font_size);
    }
    my $sum_y_start=$top;
    open IN,"<","$pathway/$nv_cg/Pvalue.table" or die;
	my $n1_before=0;my $sum_block=0;my %h3=();
    while(my $aa=<IN>){
		my $n2_before=1;
	    chomp $aa;
		if ($aa!~/^blockid/) {
			my @ite=split/\t/,$aa;
			foreach my $ii (1..$#ite) {
				if ($ite[$ii]<$opts{p}) {
					my $pos_id=$n1_before."\t".$n2_before;
					$n2_before=$n2_before+1;
					$sum_block=$sum_block+1;
					$h3{$pos_id}=0;
				}
				else{
					$n2_before=$n2_before+1;
				}
			}
		}
		$n1_before=$n1_before+1;
	}
	open O4,">","$pathway/$nv_cg/haha";
	print O4 "$sum_block\n";
	foreach my $yuu (sort keys %h3) {
		print O4 "$yuu\n";
	}
	close IN;
	open IN,"<","$pathway/$nv_cg/$fi" or die;
	my $yu=1;
	while(<IN>){
	chomp;
	next if /^blockid/;
	my @a=split /\t/;
	my $sum_x_start=$left;
	my $sum_y_len=$ybin*$hsa_len{$a[0]};
	foreach my $n1(1..$#a){
		my $sum_x_len=$xbin*$lg_len{$sca_sorts[$n1-1]};
		#print "***********\n$sum_x_len\n";
		foreach(1..$a[$n1]){
			my $randx=$sum_x_start+$sum_x_len*rand(1);
			my $randy=$sum_y_start+$sum_y_len*rand(1);
			my $sum_h=1;my $yh=1;
			foreach my $key4 (keys %h3) {
				my @it=split/\t/,$key4;
				my $nn1=$it[0];my $nn2=$it[1];
				my $num_heng=$h2{$nn1};
				my $num_lie=$h1{$nn2};
				my @itemm1=split/\t/,$num_heng;
				my @itemm2=split/\t/,$num_lie;
				my $h1_1=$itemm1[0];my $h1_2=$itemm1[1];
				my $h2_1=$itemm2[0];my $h2_2=$itemm2[1];
				$sum_h=$sum_h+1;
				if ($h2_1 < $randx and $randx < $h2_2 and $h1_1 < $randy and $randy< $h1_2) {
					$svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#e56161');
					$yh=2;
				}elsif ($sum_h>$sum_block) {
					if ($yh==1) {
						$svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#3399CC');
					}
				}else{$yu=$yu+1;}
			}
		}
		$sum_x_start+=$sum_x_len;
		}
		$sum_y_start+=$sum_y_len;
	}
	close IN;
	$svg->text('x', 200,'y',85, '-cdata',$ci_nump,'font-family',$font_family,'font-size',25);
	$svg->text('x', 600,'y',85, '-cdata',$ci_num2,'font-family',$font_family,'font-size',25);
	$svg->text('x', 450,'y',950, '-cdata',$opts{name2},'font-family',$font_family,'font-size',25);
	$svg->text('x',43,'y',500, '-cdata',$opts{name1},'font-family',$font_family,'font-size',25,'transform','rotate(270,43 500)');
	open OUT,">","$pathway/$nv_cg/$name.significance_test.dot.svg" or die;
	print OUT $svg->xmlify;
	system "rm $pathway/$nv_cg/haha";
	system "rm $pathway/$nv_cg/sum-table";
	system "rm $pathway/$nv_cg/Pvalue-result.txt";
}


sub plot_dot_hd22{
	my $path_block=$_[0];my $spe1=$_[1];my $name=$_[2];my $spe2=$_[3];my $path_geno=$_[4];my $scr_p=$_[5];
	$sum_before=0;
	open I,"<$pathway/$nv_cg/$name.scaff.match.table";
	open O,">$pathway/$nv_cg/result";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			print O "$a\tsumrow\n";
		}
		else{
			@items=split/\t/,$a;
			foreach $i (1..$#items) {
				$sum=$items[$i]+$sum_before;
				$sum_before=$sum;
			}
			print O "$a\t$sum_before\n";
		}
		$sum_before=0;
	}
	close I;
	close O;
	system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/result >  $pathway/$nv_cg/result-reverse";
	$sum_before=0;
	open I,"< $pathway/$nv_cg/result-reverse";
	open O,">$pathway/$nv_cg/sum-result";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			print O "$a\tsumcol\n";
		}
		else{
			@items=split/\t/,$a;
			foreach $i (1..$#items) {
				$sum=$items[$i]+$sum_before;
				$sum_before=$sum;
			}
			print O "$a\t$sum_before\n";
		}
		$sum_before=0;
	}
	close I;
	close O;
	system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/sum-result > $pathway/$nv_cg/sum-table";
	system "rm $pathway/$nv_cg/result";
	system "rm $pathway/$nv_cg/sum-result";
	system "rm $pathway/$nv_cg/result-reverse";
	system "Rscript $scr_p/Fish.R $pathway/$name $pathway/$nv_cg/sum-table";
	open I,"<$pathway/$nv_cg/$name.scaff.match.table";
	open O,">$pathway/$nv_cg/first_line";
	while ($a=<I>){
		chomp $a;
		if ($a=~/blockid/) {
			$a=~s/blockid//g;
			$a=~s/^\t//g;
			print O "$a";
		}
	}
	print O "\n";
	close I;
	close O;
	open M,"<$pathway/$nv_cg/Pvalue-result.txt";
	open O1,">$pathway/$nv_cg/Pvalue-result1.txt";
	while ($a=<M>){
		chomp $a;
		$a=~s/ //g;
		print O1 "$a\n";
	}
	close M;
	close O1;
	system "cat $pathway/$nv_cg/first_line $pathway/$nv_cg/Pvalue-result1.txt > $pathway/$nv_cg/Pvalue-result2.txt";
	system "awk '{print \$1}' $pathway/$nv_cg/$name.scaff.match.table >$pathway/$nv_cg/line1";
	system "paste $pathway/$nv_cg/line1 $pathway/$nv_cg/Pvalue-result2.txt > $pathway/$nv_cg/Pvalue.table";
	system "rm $pathway/$nv_cg/Pvalue-result1.txt";
	system "rm $pathway/$nv_cg/first_line";
	system "rm $pathway/$nv_cg/Pvalue-result2.txt";
	system "rm $pathway/$nv_cg/line1";


    my $svg;
    my ($left,$right)=(100,100);
    my ($top,$bottom)=(100,100);
    my ($width,$height)=(800,800);
    $svg=SVG->new('width',$left+$width+$right,'height',$top+$height+$bottom);
    #my $font = FontSize->new();
    my %sca_arr=();
    my @sca_sorts=();
    open IN,"<","$pathway/$nv_cg/$name.scaff.match.table" or die;
        while(<IN>){
	        chomp;
	        if(/^blockid/){
		        my @a=split /\t/;		
		        foreach my $fa(@a){
			        $sca_arr{$fa}=1;
			        push @sca_sorts,$fa if ($fa ne 'blockid');
		        }
	         }
        }
        close IN;
        my %hsa_len=();
        open IN,"<","$path_block" or die;
        while(<IN>){
	        chomp;
	        my @a=split /\t/;
	        $hsa_len{$a[0]}+=($a[3]-$a[2]+1);
        }
        close IN;
        my %lg_len=();
        open IN,"<","$pathway_genome.len" or die;
        while(<IN>){
	        chomp;
	        my @a=split /\t/;
	        $lg_len{$a[0]}=$a[1];	
        }
        close IN;
        my $hsa_sum_lens=0;
        foreach my $k1(keys%hsa_len){
	        $hsa_sum_lens+=$hsa_len{$k1};
        }
        my $lg_sum_lens=0;
        foreach my $k2(keys%lg_len){
	        $lg_sum_lens+=$lg_len{$k2} if (exists $sca_arr{$k2});
        }

        my $ybin=$height/$hsa_sum_lens;
        my $xbin=$width/$lg_sum_lens;
        #print "xbin:$xbin\n";

        my $plot_y_lens=$top;
        my $font_family = "Arial";
        my $font_size = 12;
		my %h2=();my %h1=();
	my $plot_y_lens_before=0;
	my $num2=1;

        $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
        foreach my $q1(sort {$a<=>$b}keys %hsa_len){
	        #print $q1."\n";
	        my $plot_start=$plot_y_lens;
            my $binn=$hsa_len{$q1}*$ybin;
	        $plot_y_lens+=$hsa_len{$q1}*$ybin;
	        my $name_height = 20; #$font->stringHeight($font_size);
                my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q1);
	        $svg->line('x1',$left,'y1',$plot_y_lens,'x2',$left+$width,'y2',$plot_y_lens,'stroke','grey','stroke-width',0.5);
	        	$h2{$num2}=$plot_y_lens_before."\t".$plot_y_lens;
			$plot_y_lens_before=$plot_y_lens;
			$num2=$num2+1;
			$svg->text('x', $left - 10 - $bac_name_width,'y',$plot_start+($binn-$name_height)/2, '-cdata',$q1,'font-family',$font_family,'font-size',$font_size);
        }
		my $num1=1;my $plot_x_lens_before=0;
        my $plot_x_lens=$left;
        $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
        foreach my $q2(@sca_sorts){
	        my $plot_start=$plot_x_lens;
	        my $binn=$lg_len{$q2}*$xbin;
	        $plot_x_lens+=$lg_len{$q2}*$xbin;
	        my $bac_name_width = 5; #$font->stringWidth($font_family,$font_size,$q2);
	        $svg->line('x1',$plot_x_lens,'y1',$top,'x2',$plot_x_lens,'y2',$top+$height,'stroke','grey','stroke-width',0.5);
	    	$h1{$num1}=$plot_x_lens_before."\t".$plot_x_lens;
		$plot_x_lens_before=$plot_x_lens;
		$num1=$num1+1;
			$svg->text('x', $plot_start+($binn-$bac_name_width)/2,'y',$top+$height+20, '-cdata',$q2,'font-family',$font_family,'font-size',$font_size);
        }
        my $sum_y_start=$top;
        open IN,"<","$pathway/$nv_cg/Pvalue.table" or die;
		my $n1_before=0;my $sum_block=0;my %h3=();
		while(my $aa=<IN>){
		my $n2_before=1;
	    chomp $aa;
		if ($aa!~/^blockid/) {
			my @ite=split/\t/,$aa;
			foreach my $ii (1..$#ite) {
				if ($ite[$ii]<$opts{p}) {
					my $pos_id=$n1_before."\t".$n2_before;
					$n2_before=$n2_before+1;
					$sum_block=$sum_block+1;
					$h3{$pos_id}=0;
				}
				else{
					$n2_before=$n2_before+1;
				}
			}
		}
		$n1_before=$n1_before+1;
	}


	open O4,">","$pathway/$nv_cg/haha";
	print O4 "$sum_block\n";
	foreach my $yuu (sort keys %h3) {
		print O4 "$yuu\n";
	}	


	close IN;
	open IN,"<","$pathway/$nv_cg/$name.scaff.match.table" or die;
	my $yu=1;
	while(<IN>){
		chomp;
		next if /^blockid/;
		my @a=split /\t/;
		my $sum_x_start=$left;
		my $sum_y_len=$ybin*$hsa_len{$a[0]};
		foreach my $n1(1..$#a){
			my $sum_x_len=$xbin*$lg_len{$sca_sorts[$n1-1]};
			#print "***********\n$sum_x_len\n";
			foreach(1..$a[$n1]){
				my $randx=$sum_x_start+$sum_x_len*rand(1);
				my $randy=$sum_y_start+$sum_y_len*rand(1);
				my $sum_h=1;my $yh=1;
				foreach my $key4 (keys %h3) {
					my @it=split/\t/,$key4;
					my $nn1=$it[0];my $nn2=$it[1];
					my $num_heng=$h2{$nn1};
					my $num_lie=$h1{$nn2};
					my @itemm1=split/\t/,$num_heng;
					my @itemm2=split/\t/,$num_lie;
					my $h1_1=$itemm1[0];my $h1_2=$itemm1[1];
					my $h2_1=$itemm2[0];my $h2_2=$itemm2[1];
					$sum_h=$sum_h+1;
					if ($h2_1 < $randx and $randx < $h2_2 and $h1_1 < $randy and $randy< $h1_2) {
						$svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#e56161');
						$yh=2;
					}elsif ($sum_h>$sum_block) {
						if ($yh==1) {
							$svg->circle('cx',$randx,'cy',$randy,'r',2,'fill','#3399CC');
						}
					}else{$yu=$yu+1;}
				}
			}
			$sum_x_start+=$sum_x_len;
		}
		$sum_y_start+=$sum_y_len;
	}
	close IN;
	$svg->text('x', 200,'y',85, '-cdata',$ci_nump,'font-family',$font_family,'font-size',25);
	$svg->text('x', 600,'y',85, '-cdata',$ci_num2,'font-family',$font_family,'font-size',25);
	$svg->text('x', 450,'y',950, '-cdata',$opts{name2},'font-family',$font_family,'font-size',25);
	$svg->text('x',43,'y',500, '-cdata',$opts{name1},'font-family',$font_family,'font-size',25,'transform','rotate(270,43 500)');
	open OUT,">","$pathway/$nv_cg/$name.significance_test.dot.svg" or die;
	print OUT $svg->xmlify;
	system "rm $pathway/$nv_cg/haha";
	system "rm $pathway/$nv_cg/sum-table";
	system "rm $pathway/$nv_cg/Pvalue-result.txt";
}


getseqpair($opts{c});
#system"perl get_all_gene_and_dia_pair_seq.pl";
sub getseqpair{
	system "mkdir $opts{o}/$nv_cg/Macrosyn_genes";
	$c_chmun=$_[0];
	if ($c_chmun==1) {
		id_chr_match($pathway_gffm,$pathway_block);
		#system "perl  id_chr_match.pl -g /public/home/ws/ws_data/yuhw/synteny/connect86/gff/PY.chr.gff -b /public/home/ws/ws_data/yuhw/synteny/connect86/PY.block -o try";
		sub id_chr_match{
		my $gff1=$_[0];my $bl=$_[1];my %h1;my %h2;
		open G,"<$gff1";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		while ($a=<G>){
			chomp $a;
			my @items=split/\t/,$a;
			open B,"<$bl";
			while ($b=<B>){
				chomp $b;
				my @items2=split/\t/,$b;
				if ($items[0] eq $items2[1]) {
					if ($items2[2]<=$items[2] and $items2[3]>=$items[3]) {
						print O "$items2[0]\t$items[1]\n";
					}
				}
			}
			close B;
		}
		close G;close O;
		}
		c_scaffold2();
		#system "perl /public/home/wangrj/synteny/testDB/tiqv/scripts/c123-scaffold2.pl -i /public/home/wangrj/synteny/testDB/plot/BF/BF-AC/AC_chr2scaf.info";
		sub c_scaffold2{
		open I1,"<$opts{o}/$nv_cg/$spename2\_chr2scaf.info"; 
		$lv=1;
		while(<I1>){
			chomp;
			my @ar=split/\t/,$_;
			if ($ar[0]=~/^c1-chr/ ) { 
				for ($i=1;$i<=$#ar;$i++) {
					$h{$ar[0]}.="\t"."$ar[$i]";
				}
			} elsif ($ar[0]=~/^c2-chr/ ) { 
				$lv=2;
				for ($j=1;$j<=$#ar;$j++) {
					$h2{$ar[0]}.="$ar[$j]"."\t";
				}
			} elsif ($ar[0]=~/^c3-chr/ ) {
				$lv=3;
			}
		}
		close I1;
		open I2,"<$opts{o}/$nv_cg/$spename2\_chr2scaf.info"; 
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/scaf_out"; 
		while(<I2>){
			chomp;
			my @ar2=split/\t/,$_;
			if ( ($ar2[0]=~/^c1-chr/) && ($lv==1)) {
				print O "$ar2[0]$h{$ar2[0]}\n";
			} elsif ( ($ar2[0]=~/^c2-chr/) && $lv==2) { 
				print O "$ar2[0]";
				for ($k=1;$k<=$#ar2;$k++) {
					print O "$h{$ar2[$k]}";
				}
				print O "\n";
			} elsif (($ar2[0]=~/^c3-chr/) && ($lv==3) ) { 
				print O "$ar2[0]";
				for ($l=1;$l<=$#ar2;$l++) {
					@ar3=split/\t/,$h2{$ar2[$l]};
					for ($m=0;$m<=$#ar2;$m++) {
						print O "$h{$ar3[$m]}";
					} 
				}
				print O "\n";
			}
		}
		close I2;
		close O;
		}
		bf($pathway_pepn);
		#system "perl BF.pl -i /public/home/wangrj/synteny/testDB/plot/BF/BF-AC/BF_AC.pairirs2 -m ../BFpair.txt -o BF_AC.pair3";
		sub bf{
		my $a;my @items;my $b;my $b2;my $b3;my %h1;
		open M,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		while ($a=<M>){
			chomp $a;
			my @items=split/\t/,$a;
			$h1{$items[1]}=$items[0];
		}
		open I,"<$opts{o}/$nv_cg/$spename1\_$spename2.pairs2";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/pair3";
		while ($b=<I>){
			chomp $b;
			my @items2=split/\t/,$b;
			$b2=$items2[0];
			$b3=$items2[1];
			if (exists $h1{$b2}){
				$b4=$h1{$b2};
				print O "$b2\t$b3\t$b4\n";
			}
		}
		close M;
		close I;
		close O;
		}
		id_chr_match2();
		#system "perl  id_chr_match2.pl -g /public/home/ws/ws_data/yuhw/synteny/connect86/gff/CG.gff -o try";
		sub id_chr_match2{
		open G,"<$pathway_gffn";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		while ($a=<G>){
			chomp $a;
			my @items=split/\t/,$a;
			print O "$items[0]\t$items[1]\n";
		}
		close G;
		close O;
		}
		cg();
		#system "perl /public/home/wangrj/synteny/testDB/tiqv/scripts/CG2.pl -i /public/home/wangrj/synteny/testDB/tiqv/pairir/PCpair -m BF_AC_chr2scaf_info_out -n BF_AC.pairir3 -o BF-ACdot";
		sub cg{
		my $a;my $b;my $b2;my $b3;my %h1;my %h2;my %h3;my $c;my $value;my $result;my $result2;my $n;my $i2;
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
		print O "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
		while ($a=<I>){
			chomp $a;
			my @items=split/\t/,$a;
			$h1{$items[1]}=$items[0];
		}
		open M,"<$opts{o}/$nv_cg/Macrosyn_genes/scaf_out";
		while ($c=<M>){
			chomp $c;
			my @items2=split/\t/,$c;
			foreach my $ii(@items2) {
				$h2{$ii}=$items2[0]; 
			}
		}
		foreach my $b3 (sort keys %h1) {
			$value=$h1{$b3};
			if (exists $h2{$value}){
				$h3{$b3}=$h2{$value};
			}
		}
		open N,"<$opts{o}/$nv_cg/Macrosyn_genes/pair3";
		while ($n=<N>){
			chomp $n;
			my @items4=split/\t/,$n;
			$i2=$items4[1];
			if (exists $h3{$i2}){
				$result=$items4[0]."\t".$items4[1]."\t".$items4[2];
				$result2=$h3{$i2};
				if ($result2=~/c(\d)-chr/) {
					$hahaa=$_;
					$result2=~s/$hahaa//;
				}
				print O "$result\t$result2\n";
			}
		}
		close M;
		close I;
		close N;
		close O;
		}
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/scaf_out";
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/pair3";
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		match_seqs($opts{m},$opts{n});
		#perl /public/home/wangrj/synteny/testDB/tiqv/scripts/match_seqs.pl -i BF-ACdot
		sub match_seqs{
		$sp1=$_[0];
		$sp2=$_[1];
		my %h1;my %h2;my @items;my $id;
		open P1,"<$pathway_pepm";
		while ($a=<P1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$value=$a;
				$h1{$id}.=$value;
			}
		}
		open P2,"<$pathway_pepn";
		while ($a=<P2>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$value=$a;
				$h2{$id}.=$value;
			}
		}
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
		open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_all_gene_pairs.pep";
		open O2,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{n}_all_gene_pairs.pep";
		while ($b=<I>){
			chomp $b;
			@items=split/\t/,$b;
			if (exists $h1{$items[0]}) {
				print O1 ">$items[0]\n$h1{$items[0]}\n";
			}
			if (exists $h2{$items[1]}) {
				print O2 ">$items[1]\n$h2{$items[1]}\n";
			}
		}
		close I;close P1;close P2;close O1;close O2;
		}
		get_dia_gene();
		#system "perl get_dia_gene.pl -matrix /public/home/ws/ws_data/yuhw/synteny/connect86/result2/PY_CG/PY_CG_scaff_c3 -o /public/home/ws/ws_data/yuhw/synteny/connect86/result/tiqv";
		sub get_dia_gene{
			open Y,"<$pathway/$nv_cg/Pvalue.table";
			my %name_chr=();
			while ($y=<Y>){
				chomp $y;
				if ($y=~/^blockid/){
					my @iee=split/\t/,$y;
					foreach my $yuew (1..$#iee) {
						$name_chr{$yuew}=$iee[$yuew];
					}
				}
				else {
					my @iee=split/\t/,$y;
					foreach my $yuew (1..$#iee) {
						if ($iee[$yuew]<=$opts{p}) {
							$ny=$iee[0];
							$uak=$name_chr{$yuew};
							$uak=~s/c1-chr//;
							$uak=~s/c2-chr//;
							$uak=~s/c3-chr//;
							$h_y1{$ny}{$uak}=0;
						}
					}
				}
			}
			close Y;
			system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $opts{o}/$nv_cg/$final > $opts{o}/$nv_cg/Macrosyn_genes/try";
			my $max=0;my $num_max=0;my @items;my $num;my $n;my %h1;my %h2;my $fir;my @items2;my %h3;
			open IN,"<$opts{o}/$nv_cg/Macrosyn_genes/try";
			open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_diagonal_gene_pairs.dot";
			open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_significant_gene_pairs.dot";
			print O "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
			print O1 "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
			while (my $a=<IN>){
				chomp $a;
				if ($a=~/^blockid/) {
					@items=split/\t/,$a;
					$fir=0;
					$h1{$fir}=$items[0];
				}
				else {
					@items=split/\t/,$a;
					$num=@items;
					for($n=1;$n<$num;$n=$n+1){
						if ($items[$n]>$max) {
							$max=$items[$n];
							$num_max=$n;
						}
						else{
							$max=$max;
						}
					}
					$items[0]=~s/c1-chr//;
					$items[0]=~s/c2-chr//;
					$items[0]=~s/c3-chr//;
					$h1{$num_max}=$items[0];
					#print  "$num_max\t$items[0]\n";
					$h2{$num_max."\t".$items[0]}=0;
					$max=0;
					#print "$num_max\t$items[0]\n";
				}
			}
			close IN;
			system "rm $opts{o}/$nv_cg/Macrosyn_genes/try";
			open IN2,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
			while (my $a=<IN2>){
				chomp $a;
				@items2=split/\t/,$a;
				if (exists $h2{$items2[2]."\t".$items2[3]}) {
					print O "$a\n";
				}
				if (exists $h_y1{$items2[2]} and exists $h_y1{$items2[2]}{$items2[3]}) {
					print O1 "$a\n";
				}
			}
			close IN2;
			close O;
			close O1;
		}
		match_seqs2($opts{m},$opts{n});
		#perl /public/home/wangrj/synteny/testDB/tiqv/scripts/match_seqs.pl -i BF-ACdot
		sub match_seqs2{
		$sp1=$_[0];
		$sp2=$_[1];
		my %h1;my %h2;my @items;my $id;
		open P1,"<$pathway_pepm";
		while ($a=<P1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$v=$a;
				$h1{$id}.=$v;
			}
		}
		open P2,"<$pathway_pepn";
		while ($a=<P2>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$v=$a;
				$h2{$id}.=$v;
			}
		}
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_diagonal_gene_pairs.dot";
		open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_diagonal_gene_pairs.pep";
		open O2,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{n}_diagonal_gene_pairs.pep";
		while ($b=<I>){
			chomp $b;
			@items=split/\t/,$b;
			if (exists $h1{$items[0]}) {
				print O1 ">$items[0]\n$h1{$items[0]}\n";
			}
			if (exists $h2{$items[1]}) {
				print O2 ">$items[1]\n$h2{$items[1]}\n";
			}
		}
		close I;close P1;close P2;close O1;close O2;
		}
	}


	#开始新的
	if ($c_chmun==2) {
		open Y,"<$pathway/$nv_cg/Pvalue.table";
		my %name_chr=();
		while ($y=<Y>){
			chomp $y;
			if ($y=~/^blockid/){
				my @iee=split/\t/,$y;
				foreach my $yuew (1..$#iee) {
					$name_chr{$yuew}=$iee[$yuew];
				}
			}
			else {
				my @iee=split/\t/,$y;
				foreach my $yuew (1..$#iee) {
					if ($iee[$yuew]<=$opts{p}) {
						$ny=$iee[0];
						$uak=$name_chr{$yuew};
						$h_y2{$ny}{$uak}=0;
					}
				}
			}
		}
		close Y;
		id_chr_match1($pathway_gffm,$pathway_block);
		#system "perl  id_chr_match.pl -g /public/home/ws/ws_data/yuhw/synteny/connect86/gff/PY.chr.gff -b /public/home/ws/ws_data/yuhw/synteny/connect86/PY.block -o try";
		sub id_chr_match1{
		my $gff1=$_[0];my $bl=$_[1];my %h1;my %h2;
		open G,"<$gff1";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		while ($a=<G>){
			chomp $a;
			my @items=split/\t/,$a;
			open B,"<$bl";
			while ($b=<B>){
				chomp $b;
				my @items2=split/\t/,$b;
				if ($items[0] eq $items2[1]) {
					if ($items2[2]<=$items[2] and $items2[3]>=$items[3]) {
						print O "$items2[0]\t$items[1]\n";
					}
				}
			}
			close B;
		}
		close G;close O;
		}
		bf1($pathway_pepn);
		#system "perl BF.pl -i /public/home/wangrj/synteny/testDB/plot/BF/BF-AC/BF_AC.pairirs2 -m ../BFpair.txt -o BF_AC.pair3";
		sub bf1{
		my $a;my @items;my $b;my $b2;my $b3;my %h1;my @items2;my $b4;
		open M,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		while ($a=<M>){
			chomp $a;
			@items=split/\t/,$a;
			$h1{$items[1]}=$items[0];
		}
		open I,"<$opts{o}/$nv_cg/$nv_cg.pairs2";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/pair3";
		while ($b=<I>){
			chomp $b;
			@items2=split/\t/,$b;
			$b2=$items2[0];
			$b3=$items2[1];
			if (exists $h1{$b2}){
				$b4=$h1{$b2};
				print O "$b2\t$b3\t$b4\n";
			}
		}
		close M;
		close I;
		close O;
		}
		id_chr_match21();
		#system "perl  id_chr_match2.pl -g /public/home/ws/ws_data/yuhw/synteny/connect86/gff/CG.gff -o try";
		sub id_chr_match21{
		my @items;
		open G,"<$pathway_gffn";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		while ($a=<G>){
			chomp $a;
			@items=split/\t/,$a;
			print O "$items[0]\t$items[1]\n";
		}
		close G;
		close O;
		}
		cg1();
		#system "perl /public/home/wangrj/synteny/testDB/tiqv/scripts/CG2.pl -i /public/home/wangrj/synteny/testDB/tiqv/pairir/PCpair -m BF_AC_chr2scaf_info_out -n BF_AC.pairir3 -o BF-ACdot";
		sub cg1{
		my $a;my @items;my $b;my $b2;my $b3;my %h1;my @items2;my %h2;my %h3;my $c;my $value;my @items4;my $result;my $result2;my $n;my $i2;
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
		print O "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
		while ($a=<I>){
			chomp $a;
			@items=split/\t/,$a;
			$h1{$items[1]}=$items[0];
		}
		open M,"<$opts{o}/$nv_cg/$spename2\_chr2scaf.info";
		while ($c=<M>){
			chomp $c;
			@items2=split/\t/,$c;
			foreach my $ii(@items2) {
				$h2{$ii}=$items2[0]; 
			}
		}
		foreach my $b3 (sort keys %h1) {
			$value=$h1{$b3};
			if (exists $h2{$value}){
				$h3{$b3}=$h2{$value};
			}
		}
		open N,"<$opts{o}/$nv_cg/Macrosyn_genes/pair3";
		while ($n=<N>){
			chomp $n;
			@items4=split/\t/,$n;
			$i2=$items4[1];
			if (exists $h3{$i2}){
				$result=$items4[0]."\t".$items4[1]."\t".$items4[2];
				$result2=$h3{$i2};
				if ($result2=~/c(\d)-chr/) {
					$result2=~s/$_//;
				}
				print O "$result\t$result2\n";
			}
		}
		close M;
		close I;
		close N;
		close O;
		}
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/$opts{m}.pair";
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/pair3";
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/spe2.pair";
		match_seqs1($opts{m},$opts{n});
		#perl /public/home/wangrj/synteny/testDB/tiqv/scripts/match_seqs.pl -i BF-ACdot
		sub match_seqs1{
		$sp1=$_[0];
		$sp2=$_[1];
		my %h1;my %h2;my @items;my $id;
		open P1,"<$pathway_pepm";
		while ($a=<P1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$value=$a;
				$h1{$id}.=$value;
			}
		}
		open P2,"<$pathway_pepn";
		while ($a=<P2>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$value=$a;
				$h2{$id}.=$value;
			}
		}
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
		open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_all_gene_pairs.pep";
		open O2,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{n}_all_gene_pairs.pep";
		while ($b=<I>){
			chomp $b;
			@items=split/\t/,$b;
			if (exists $h1{$items[0]}) {
				print O1 ">$items[0]\n$h1{$items[0]}\n";
			}
			if (exists $h2{$items[1]}) {
				print O2 ">$items[1]\n$h2{$items[1]}\n";
			}
		}
		close I;close P1;close P2;close O1;close O2;
		}
		get_dia_gene1();
		#system "perl get_dia_gene.pl -matrix /public/home/ws/ws_data/yuhw/synteny/connect86/result2/PY_CG/PY_CG_scaff_c3 -o /public/home/ws/ws_data/yuhw/synteny/connect86/result/tiqv";
		sub get_dia_gene1{
		system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $opts{o}/$nv_cg/$nv_cg.scaff.match.table > $opts{o}/$nv_cg/Macrosyn_genes/try";
		my $max=0;my $num_max=0;my @items;my $num;my $n;my %h1;my %h2;my $fir;my @items2;my %h3;
		open IN,"<$opts{o}/$nv_cg/Macrosyn_genes/try";
		open O,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_diagonal_gene_pairs.dot";
		print O "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
		open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_significant_gene_pairs.dot";
		print O1 "Ancestor-geneID\tspecies-geneID\tAncestor-chr\tspecies-chr\n";
		while (my $a=<IN>){
			chomp $a;
			if ($a=~/^blockid/) {
				@items=split/\t/,$a;
				$fir=0;
				$h1{$fir}=$items[0];
			}
			else {
				@items=split/\t/,$a;
				$num=@items;
				for($n=1;$n<$num;$n=$n+1){
					if ($items[$n]>$max) {
						$max=$items[$n];
						$num_max=$n;
					}
					else{
						$max=$max;
					}
				}
				$items[0]=~s/c1-chr//;
				$items[0]=~s/c2-chr//;
				$items[0]=~s/c3-chr//;
				$h1{$num_max}=$items[0];
				#print  "$num_max\t$items[0]\n";
				$h2{$num_max."\t".$items[0]}=0;
				$max=0;
			}
		}
		close IN;
		system "rm $opts{o}/$nv_cg/Macrosyn_genes/try";
		open IN2,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
		while (my $a=<IN2>){
			chomp $a;
			@items2=split/\t/,$a;
			if (exists $h2{$items2[2]."\t".$items2[3]}) {
				print O "$a\n";
			}
			if (exists $h_y2{$items2[2]} and exists $h_y2{$items2[2]}{$items2[3]}) {
				print O1 "$a\n";
			}
		}
		close IN2;
		close O;
		close O1;
		}
		match_seqs21($opts{m},$opts{n});
		#perl /public/home/wangrj/synteny/testDB/tiqv/scripts/match_seqs.pl -i BF-ACdot
		sub match_seqs21{
		$sp1=$_[0];
		$sp2=$_[1];
		my %h1;my %h2;my @items;my $id;
		open P1,"<$pathway_pepm";
		while ($a=<P1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$v=$a;
				$h1{$id}.=$v;
			}
		}
		open P2,"<$pathway_pepn";
		while ($a=<P2>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$id=$1;
			}
			else{
				$v=$a;
				$h2{$id}.=$v;
			}
		}
		open I,"<$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_$opts{n}_diagonal_gene_pairs.dot";
		open O1,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{m}_diagonal_gene_pairs.pep";
		open O2,">$opts{o}/$nv_cg/Macrosyn_genes/$opts{n}_diagonal_gene_pairs.pep";
		while ($b=<I>){
			chomp $b;
			@items=split/\t/,$b;
			if (exists $h1{$items[0]}) {
				print O1 ">$items[0]\n$h1{$items[0]}\n";
			}
			if (exists $h2{$items[1]}) {
				print O2 ">$items[1]\n$h2{$items[1]}\n";
			}
		}
		close I;close P1;close P2;close O1;close O2;
		}
	}

}

open I,"<","$pathway/$nv_cg/Pvalue.table" or die;
while ($a=<I>){
	chomp $a;
	if ($a!~/^blockid/) {
		$n=0;
		my @items=split/\t/,$a;
		foreach my $j (1..$#items) {
			if ($items[$j]<=$opts{p}) {
				$n=$n+1;
				$h1a{$items[0]}=$n;
			}
		}
	}
}
close I;
system "awk -F '\t'  '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS \$i:\$i}END{for(i=0;i++<NF;)print a[i]}' $pathway/$nv_cg/Pvalue.table > $pathway/$nv_cg/Pvalue.table-reverse";

open I,"<$pathway/$nv_cg/Pvalue.table-reverse";
while ($b=<I>){
	chomp $b;
	if ($b!~/^blockid/) {
		my @it=split/\t/,$b;
		foreach my $j2 (1..$#it) {
			if ($it[$j2]<=$opts{p}) {
				if (exists $h2a{$it[0]}) {
					$h2a{$it[0]}=$h2a{$it[0]}."\t".$j2;
				}
				else{$h2a{$it[0]}=$j2;}
			}
		}
	}
}
close I;
system "mkdir $pathway/$nv_cg/chr_breakage_fusion_result";
open O,">$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R1";
#print O "Chr\tNumber\tAncestor\n";

foreach $i (keys %h2a) {
	if ($h2a{$i}=~/\t/) {
		my @num=split/\t/,$h2a{$i};
		foreach $nu (@num) {
			if (exists $h1a{$nu}) {
				if ($h1a{$nu}==1) {
					$chu=1/$h1a{$nu};
					print O "$i\t$chu\t$nu\_all\n";
				}
				else{
					$chu=1/$h1a{$nu};
					$chu=sprintf "%.5f",$chu;
					print O "$i\t$chu\t$nu\_part\n";			
				}
			}
		}
	}
	else{
		$nu2=$h2a{$i};
		if (exists $h1a{$nu2}) {
			$chu=1/$h1a{$nu2};
			if ($h1a{$nu2}==1) {
				print O "$i\t$chu\t$h2a{$i}\_all\n";
			}
			else{
				print O "$i\t$chu\t$h2a{$i}\_part\n";			
			}
		}
	}
}
close O;

system "sort -Vk1,1 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R1 >$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R2";

open O,">$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R3";
print O "Chr\tNumber\tAncestor\n";
close O;
system "sort -Vk3,3 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R1 >$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R4";
system "cat $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R3 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R4 >$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R5";

system "cat $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R3 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R2 >$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R";

open I,"<$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R";
open O,">$pathway/$nv_cg/chr_breakage_fusion_result/$opts{n}_chr_breakage_fusion.result";

while ($c=<I>){
	chomp $c;
	@ai=split/\t/,$c;
	if (exists $h3a{$ai[0]}) {
		$h3a{$ai[0]}=$h3a{$ai[0]}."\t".$ai[2];
	}
	else{$h3a{$ai[0]}=$ai[2];}
}
close I;
foreach $ia (sort keys %h3a) {
	print O "$ia\t$h3a{$ia}\n";
}
close O;

$h_nv=0;
open I,"<$pathway_block";
while (my $c=<I>){
	chomp $c;
	@ai2=split/\t/,$c;
	if (!exists $h_nv_num{$ai2[0]}) {
		$h_nv_num{$ai2[0]}=0;
		$h_nv=$h_nv+1;
	}
}
close I;
open I,"<$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R";
while (my $c=<I>){
	chomp $c;
	if ($c!~/Ancestor$/) {
		my @ai2=split/\t/,$c;
		if (exists $h_chr_num2{$ai2[0]}) {
			$h_chr_num2{$ai2[0]}=$h_chr_num2{$ai2[0]}+1;
		}
		else{
			$h_chr_num2{$ai2[0]}=1;
		}
	}
}
close I;
open I,"<$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R";
open O,">$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R6";
while (my $c=<I>){
	chomp $c;
	if ($c!~/Ancestor$/) {
		my @ai2=split/\t/,$c;
		if (exists $h_chr_num2{$ai2[0]}) {
			$bizhi=1/$h_chr_num2{$ai2[0]};
			print O "$ai2[0]\t$bizhi\t$ai2[2]\n";
		}
	}
}
close I;
close O;
system "cat $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R3 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R6 >$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R7";

system "rm $pathway/$nv_cg/Pvalue.table-reverse";

open I,"<$pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R7";
open O,">$opts{o}/formatted_ancestor_chr_color.txt";
my %h_un=();my %h_numpai=();

while (my $b=<I>){
	chomp $b;
	my @it=split/\t/,$b;
	if ($it[2] ne "Ancestor") {
		my @it2=split/_/,$it[2];
		if (exists $anchr_col{$it2[0]}) {
			$h_numpai{$it2[0]}=$it[2];
			if (!exists $h_un{$it[2]}) {
				$h_un{$it[2]}=$anchr_col{$it2[0]};
			}
		}
	}
}		
close I;

my $nnn=1;
foreach my $i (sort keys %h_numpai) {
	$ss=$h_numpai{$i};
	$ssa=qq{"$ss"};
	if ($nnn==1) {
		$sum_color1=$ssa."="."$h_un{$ss}";
	}
	else{
		$sum_color2=$ssa."="."$h_un{$ss}";
		$sum_color=$sum_color1.",".$sum_color2;
		$sum_color1=$sum_color;
	}
	$nnn=$nnn+1;
	print O "$ss\t$h_un{$ss}\n";
}
close O;

system "Rscript $opts{s}/Chr_breakage_fusion.R $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R $pathway/$nv_cg/chr_breakage_fusion_result $opts{m}_$opts{n} $h_nv $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R5 $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R7 $opts{o}/formatted_ancestor_chr_color.txt";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R1";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R2";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R3";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R4";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R5";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R6";
system "rm $pathway/$nv_cg/chr_breakage_fusion_result/chr_breakage_fusion_result_for_R7";

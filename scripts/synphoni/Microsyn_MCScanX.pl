#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'getcwd';

my %opts;
GetOptions(\%opts,"i=s","o=s","a=s","t=s","e=s","s=s","h|help","mk=s","mg=s","ms=s","me=s","mm=s","mw=s","ma=s","mb=s");
if(!(defined $opts{i} and defined $opts{o})){
	die"**********************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-o	Full path to the [outputDir] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5)
	-s	The number of best non-self alignment hits (default:5)
	-mk	MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_PENALTY (default:50)
	-mg	GAP_PENALTY, gap penalty (default:-1)
	-ms	MATCH_SIZE, number of genes required to call a collinear block (default:5)
	-me	E_VALUE, alignment significance (default:1e-5/0.001)
	-mm	MAX_GAPS, maximum gaps allowed (default:25)
	-mw	OVERLAP_WINDOW, maximum distance (of genes) to collapse BLAST matches (default:5)
	-ma	1:Not only builds the pairwise blocks; 2:Only builds the pairwise blocks(default:1)
	-mb	Patterns of collinear blocks. 0:intra- and inter-species (default); 1:intra-species; 2:inter-species
		*********************************************\n";
}

if (defined $opts{h} or defined $opts{help}) {
	die"**********************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-o	Full path to the [outputDir] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5)
	-s	The number of best non-self alignment hits (default:5)
	-mk	MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_PENALTY (default:50)
	-mg	GAP_PENALTY, gap penalty (default:-1)
	-ms	MATCH_SIZE, number of genes required to call a collinear block (default:5)
	-me	E_VALUE, alignment significance (default:1e-5/0.001)
	-mm	MAX_GAPS, maximum gaps allowed (default:25)
	-mw	OVERLAP_WINDOW, maximum distance (of genes) to collapse BLAST matches (default:5)
	-ma	1:Not only builds the pairwise blocks; 2:Only builds the pairwise blocks (default:1)
	-mb	Patterns of collinear blocks. 0:intra- and inter-species; 1:intra-species; 2:inter-species (default:0)
		*********************************************\n";
}

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
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####

################################
if (!(defined $opts{mk})) {
	$opts{mk}=50;
}
if (!(defined $opts{mg})) {
	$opts{mg}=-1;
}
if (!(defined $opts{ms})) {
	$opts{ms}=5;
}

if (!(defined $opts{me})) {
	$opts{me}=1e-5;
}
if (!(defined $opts{mm})) {
	$opts{mm}=25;
}
if (!(defined $opts{mw})) {
	$opts{mw}=5;
}
if (!(defined $opts{mb})) {
	$opts{mb}=0;
}
if (!(defined $opts{ma})) {
	$opts{ma}=1;
}

#################################
if (!(defined $opts{e})) {
	if ($opts{a} eq "blast") {
		$opts{e}=1e-5;
	}
	if ($opts{a} eq "diamond") {
		$opts{e}=0.001;
	}
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{s})) {
	$opts{s}=5;
}

my $old_dir = getcwd();

open P2,">$opts{o}/Microsynteny_statistics.txt";
print P2 "Species\tAll_genes_numbers\tCollinear_genes_numbers\tMicrosyn_percentage\tBlock_numbers\tMax_score\tMin_score\tAverage_score\tMedian_score\n";

my $percentage_zhi; my $num_col; my $num_all;

if ($opts{a} eq "blast") {
	my $n=1; my %h=();
	my $dir1=$opts{i};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $file=$dir1."/$c";
			my $dirphr = "$opts{i}/$1.pep.phr";
			my $dirpin = "$opts{i}/$1.pep.pin";
			my $dirpsq = "$opts{i}/$1.pep.psq";
			if (-e $dirphr and -e $dirpin and -e $dirpsq) {
			}
			else{
				system "makeblastdb -in $file -dbtype prot";
			}
			$h{$n}=$1;
			$n=$n+1;
		}
	}
	close P;
	my $num=keys %h;
	system "mkdir $opts{o}/microsyn_genes_links";
	print "The [$opts{o}/microsyn_genes_links] directory has been created!\n";

	my $dirli1="$opts{o}/MCScanX_result";
	system "mkdir $opts{o}/MCScanX_result";
	print "The [$opts{o}/MCScanX_result] directory has been created!\n";

	my $n1;
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "mkdir $opts{o}/$h{$i}_$h{$n1}";
				print "The [$opts{o}/$h{$i}_$h{$n1}] directory has been created!\n";

				system "blastp -query $opts{i}/$h{$i}.pep -db $opts{i}/$h{$n1}.pep -outfmt 6 -out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast -evalue $opts{e} -num_threads $opts{t} -max_target_seqs $opts{s}";
				system "cat $opts{i}/$h{$i}_simplified.gff $opts{i}/$h{$n1}_simplified.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
				chdir "$opts{o}/$h{$i}_$h{$n1}";
				###
				if ($opts{ma}==1) {
					system "MCScanX $h{$i}_$h{$n1} -k $opts{mk} -g $opts{mg} -s $opts{ms} -e $opts{me} -m $opts{mm} -w $opts{mw} -b $opts{mb}";
				}
				if ($opts{ma}==2) {
					system "MCScanX $h{$i}_$h{$n1} -k $opts{mk} -g $opts{mg} -s $opts{ms} -e $opts{me} -m $opts{mm} -w $opts{mw} -b $opts{mb} -a";
				}
				chdir($old_dir);

				###Percentage start#
				open Y1,"<$opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity";
				my @score = (); my @sorted_score=();
				while (my $y1=<Y1>){
					chomp $y1;
					if ($y1=~/collinear(\s)genes:(\s)(\S+),/) {
						$num_col=$3;
					}
					if ($y1=~/all(\s)genes:(\s)(\S+)/) {
						$num_all=$3;
					}
					if ($y1=~/Percentage:(\s)(\S+)$/) {
						$percentage_zhi=$2;
					}
					if ($y1=~/score=(\S+)(\s)/) {
						push @score, $1;
					}
				}
				close Y1;
				###
				# 数组的数量
				my $count = scalar(@score);
				#print "数量：$count\n";

				# 数组的最大值
				my $max = $score[0];
				foreach my $number (@score) {
					if ($number > $max) {
						$max = $number;
					}
				}
				#print "最大值：$max\n";

				# 数组的最小值
				my $min = $score[0];
				foreach my $number (@score) {
					if ($number < $min) {
						$min = $number;
					}
				}
				#print "最小值：$min\n";

				# 数组的总和和平均值
				my $sum = 0;
				foreach my $number (@score) {
					$sum += $number;
				}
				my $average = $sum / $count;
				#print "平均值：$average\n";

				# 数组的中位数
				my $median = 0;
				@sorted_score = sort {$a <=> $b} @score;
				if ($count % 2 == 0) {
					$median = ($sorted_score[$count/2 - 1] + $sorted_score[$count/2]) / 2;
				} else { 
					$median = $sorted_score[($count-1)/2];
				}
				#print "中位数：$median\n";

				# 数组的方差和标准差
				my $variance = 0;
				foreach my $number (@score) {
					$variance += ($number - $average) ** 2;
				}
				$variance = $variance / $count;
				my $standard_deviation = sqrt($variance);
				#print "方差：$variance\n";
				#print "标准差：$standard_deviation\n";

				print P2 "$h{$i}_$h{$n1}\t$num_all\t$num_col\t$percentage_zhi\t$count\t$max\t$min\t$average\t$median\n";
				###Percentage end#

				system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/MCScanX_result";
				posi2($h{$i},$h{$n1}); 
			}
		}
	}
	sub posi2{
		my %h1=();my %h2=();
		my $spe1=$_[0];	my $spe2=$_[1];
		open B1,"<$opts{i}/$spe1\_simplified.gff";
			while (my $a=<B1>){
			chomp $a;
			my @items1=split/\t/,$a;
			$h1{$items1[1]}=$items1[0]."\t".$items1[2]."\t".$items1[3];
		}
		close B1;
		open B2,"<$opts{i}/$spe2\_simplified.gff";
		while (my $a2=<B2>){
			chomp $a2;
			my @items2=split/\t/,$a2;
			$h2{$items2[1]}=$items2[0]."\t".$items2[2]."\t".$items2[3];
		}
		close B2;

		open O,">$opts{o}/microsyn_genes_links/$spe1\_$spe2\_microsyn_genes.links";
		open I,"<$opts{o}/MCScanX_result/$spe1\_$spe2.collinearity";
		while (my $b=<I>){
			chomp $b;
			if ($b!~/^#/) {
				my @items3=split/\t/,$b;
				if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
					print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=chr6\n";
				}
				if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
					print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=chr6\n";
				}
			}
		}
		close I;
		close O;
	}
	system "cat $opts{o}/microsyn_genes_links/*_microsyn_genes.links > $opts{o}/microsyn_genes_links/All_species_microsyn_genes.links";
}


if ($opts{a} eq "diamond") {
	my $n=1;my %h=();
	my $dir1=$opts{i};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $file=$dir1."/$c";
			my $dirdb = "$opts{i}/$1.pep.database.dmnd";
			if (-e $dirdb) {
			}
			else{
				system "diamond makedb --in $file -d $file.database";
			}
			$h{$n}=$1;
			$n=$n+1;
		}
	}
	close P;
	my $num=keys %h;
	system "mkdir $opts{o}/microsyn_genes_links";
	print "The [$opts{o}/microsyn_genes_links] directory has been created!\n";

	system "mkdir $opts{o}/MCScanX_result";
	print "The [$opts{o}/MCScanX_result] directory has been created!\n";

	my $n1;
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "mkdir $opts{o}/$h{$i}_$h{$n1}";
				print "The [$opts{o}/$h{$i}_$h{$n1}] directory has been created!\n";

				system "diamond blastp --query $opts{i}/$h{$i}.pep --db $opts{i}/$h{$n1}.pep.database --outfmt 6 --out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs $opts{s}";
				system "cat $opts{i}/$h{$i}_simplified.gff $opts{i}/$h{$n1}_simplified.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
				chdir "$opts{o}/$h{$i}_$h{$n1}";
				
				###
				if ($opts{ma}==1) {
					system "MCScanX $h{$i}_$h{$n1} -k $opts{mk} -g $opts{mg} -s $opts{ms} -e $opts{me} -m $opts{mm} -w $opts{mw} -b $opts{mb}";
				}
				if ($opts{ma}==2) {
					system "MCScanX $h{$i}_$h{$n1} -k $opts{mk} -g $opts{mg} -s $opts{ms} -e $opts{me} -m $opts{mm} -w $opts{mw} -b $opts{mb} -a";
				}
				###
				chdir($old_dir);
				###Percentage start#
				open Y1,"<$opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity";
				my @score = (); my @sorted_score=();
				while (my $y1=<Y1>){
					chomp $y1;
					if ($y1=~/collinear(\s)genes:(\s)(\S+),/) {
						$num_col=$3;
					}
					if ($y1=~/all(\s)genes:(\s)(\S+)/) {
						$num_all=$3;
					}
					if ($y1=~/Percentage:(\s)(\S+)$/) {
						$percentage_zhi=$2;
					}
					if ($y1=~/score=(\S+)(\s)/) {
						push @score, $1;
					}
				}
				close Y1;
				###
				# 数组的数量
				my $count = scalar(@score);
				#print "数量：$count\n";

				# 数组的最大值
				my $max = $score[0];
				foreach my $number (@score) {
					if ($number > $max) {
						$max = $number;
					}
				}
				#print "最大值：$max\n";

				# 数组的最小值
				my $min = $score[0];
				foreach my $number (@score) {
					if ($number < $min) {
						$min = $number;
					}
				}
				#print "最小值：$min\n";

				# 数组的总和和平均值
				my $sum = 0;
				foreach my $number (@score) {
					$sum += $number;
				}
				my $average = $sum / $count;
				#print "平均值：$average\n";

				# 数组的中位数
				my $median = 0;
				@sorted_score = sort {$a <=> $b} @score;
				if ($count % 2 == 0) {
					$median = ($sorted_score[$count/2 - 1] + $sorted_score[$count/2]) / 2;
				} else { 
					$median = $sorted_score[($count-1)/2];
				}
				#print "中位数：$median\n";

				# 数组的方差和标准差
				my $variance = 0;
				foreach my $number (@score) {
					$variance += ($number - $average) ** 2;
				}
				$variance = $variance / $count;
				my $standard_deviation = sqrt($variance);
				#print "方差：$variance\n";
				#print "标准差：$standard_deviation\n";

				print P2 "$h{$i}_$h{$n1}\t$num_all\t$num_col\t$percentage_zhi\t$count\t$max\t$min\t$average\t$median\n";
				###Percentage end#

				system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/MCScanX_result";
				posi22($h{$i},$h{$n1}); 
			}
		}
	}
	sub posi22{
		my %h1=();my %h2=();
		my $spe1=$_[0];	my $spe2=$_[1];
		open B1,"<$opts{i}/$spe1\_simplified.gff";
			while ($a=<B1>){
			chomp $a;
			my @items1=split/\t/,$a;
			$h1{$items1[1]}=$items1[0]."\t".$items1[2]."\t".$items1[3];
		}
		close B1;
		open B2,"<$opts{i}/$spe2\_simplified.gff";
		while (my $a2=<B2>){
			chomp $a2;
			my @items2=split/\t/,$a2;
			$h2{$items2[1]}=$items2[0]."\t".$items2[2]."\t".$items2[3];
		}
		close B2;

		open O,">$opts{o}/microsyn_genes_links/$spe1\_$spe2\_microsyn_genes.links";
		open I,"<$opts{o}/MCScanX_result/$spe1\_$spe2.collinearity";
		while (my $b=<I>){
			chomp $b;
			if ($b!~/^#/) {
				my @items3=split/\t/,$b;
				if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
					print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=chr6\n";
				}
				if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
					print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=chr6\n";
				}
			}
		}
		close I;
		close O;
	}
	system "cat $opts{o}/microsyn_genes_links/*_microsyn_genes.links > $opts{o}/microsyn_genes_links/All_species_microsyn_genes.links";
}
close P2;
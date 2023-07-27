#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'getcwd';

my %opts;
GetOptions(\%opts,"i1=s","a=s","t=s","ca=s","ci=s","e=s","o=s","h|help","mk=s","mg=s","ms=s","me=s","mm=s","mw=s","ma=s","mb=s");
if(!(defined $opts{i1} and defined $opts{o})){

	die"**********************************************\n
	-i1	Full path to the [inputDir1] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
	-ca	The color of the lines representing conserved both mi-crosynteny and macrosynteny genes in the resulting graph (default: chr5)
	-ci	The color of the lines representing only conserved microsynteny genes in the resulting graph (default: chr24)
	-mk	MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_PENALTY (default:50)
	-mg	GAP_PENALTY, gap penalty (default:-1)
	-ms	MATCH_SIZE, number of genes required to call a microsynteny block (default:5)
	-me	E_VALUE, alignment significance (default:1e-5)
	-mm	MAX_GAPS, maximum gaps allowed (default:25)
	-mw	OVERLAP_WINDOW, maximum distance (of genes) to collapse BLAST matches (default:5)
	-ma	1:Not builds the pairwise blocks; 2:Only builds the pairwise blocks(default:1)
	-mb	Patterns of microsynteny blocks. 0:intra- and inter-species (default); 1:intra-species; 2:inter-species
	-h|-help Print this help page
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i1	Full path to the [inputDir1] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
	-ca	The color of the lines representing conserved both mi-crosynteny and macrosynteny genes in the resulting graph (default: chr5)
	-ci	The color of the lines representing only conserved microsynteny genes in the resulting graph (default: chr24)
	-mk	MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_PENALTY (default:50)
	-mg	GAP_PENALTY, gap penalty (default:-1)
	-ms	MATCH_SIZE, number of genes required to call a microsynteny block (default:5)
	-me	E_VALUE, alignment significance (default:1e-5)
	-mm	MAX_GAPS, maximum gaps allowed (default:25)
	-mw	OVERLAP_WINDOW, maximum distance (of genes) to collapse BLAST matches (default:5)
	-ma	1:Not builds the pairwise blocks; 2:Only builds the pairwise blocks(default:1)
	-mb	Patterns of microsynteny blocks. 0:intra- and inter-species (default); 1:intra-species; 2:inter-species
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

if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{e})) {
	if ($opts{a} eq "blast") {
		$opts{e}=1e-5;
	}
	if ($opts{a} eq "diamond") {
		$opts{e}=0.001;
	}
}
if (!(defined $opts{ca})) {
	$opts{ca}="chr5";
}
if (!(defined $opts{ci})) {
	$opts{ci}="chr24";
}

######################
my $old_dir = getcwd();

my $dir = "$opts{o}"; 


#######################
my $count1;my $count2;
my $count3;
my $count4;
my $overlap_n1; my $overlap_n2;  my %overlap1=(); my %overlap2=();
my %ha_num=(); my %hb_num=();

my %h=(); my $num; my $n1; my %h_gene1=(); my %h_gene2=();my %h1=(); my %h2=(); my $a2; my %gene_chr_p=();my $dir1; my $name;

my $n=1;


if ($opts{a} eq "blast") {
	$dir1=$opts{i1};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ($c!~/_/) {
			if ( $c=~/^(\S+).pep$/ ) {
				my $file=$dir1."/$c";

				my $dirphr = "$opts{i1}/$1.pep.phr";
				my $dirpin = "$opts{i1}/$1.pep.pin";
				my $dirpsq = "$opts{i1}/$1.pep.psq";
				if (-e $dirphr and -e $dirpin and -e $dirpsq) {
				}
				else{
					system "makeblastdb -in $file -dbtype prot";
				}
				$h{$n}=$1;
				$n=$n+1;
			}
		}
	}
	close P;

	$num=keys%h;
	system "mkdir $opts{o}/microsyn_genes.links";
	print "The [$opts{o}/microsyn_genes_links] directory has been created!\n";
	system "mkdir $opts{o}/collinearity";
	print "The [$opts{o}/collinearity] directory has been created!\n";
	
	open S,">$opts{o}/microsyn_genes.links/summary.txt";
	print S "SpeciesA\tSpeciesB\tA_Macrosyn_genes_numbers\tB_Macrosyn_genes_numbers\tA_Microsyn_genes_numbers\tB_Microsyn_genes_numbers\tA_Overlap_genes_numbers\tB_Overlap_genes_numbers\n";
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "mkdir $opts{o}/$h{$i}_$h{$n1}";
				print "The [$opts{o}/$h{$i}_$h{$n1}] directory has been created!\n";
				system "blastp -query $opts{i1}/$h{$i}.pep -db $opts{i1}/$h{$n1}.pep -outfmt 6 -out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast -evalue $opts{e} -num_threads $opts{t} -max_target_seqs 5";
				system "cat $opts{i1}/$h{$i}_simplified.gff $opts{i1}/$h{$n1}_simplified.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
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
				system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/collinearity";
				color1($h{$i},$h{$n1}); 
				posi21($h{$i},$h{$n1}); 
			}
		}
	}
	system "cat $opts{o}/microsyn_genes.links/*_microsyn_genes.links > $opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
	sub color1{
		my $spe1=$_[0];	my $spe2=$_[1];%h_gene1=();%h_gene2=();
		open GENE1,"<$opts{i1}/$spe1\_diagonal_gene_pairs.pep";
		while (my $a=<GENE1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$h_gene1{$1}=1;
			}
		}
		close GENE1;
		open GENE2,"<$opts{i1}/$spe2\_diagonal_gene_pairs.pep";
		while (my $b=<GENE2>){
			chomp $b;
			if ($b=~/>(\S+)/) {
				$h_gene2{$1}=1;
			}
		}
		close GENE2;
		$count1 = scalar keys %h_gene1;
		$count2 = scalar keys %h_gene2;
		print S "$spe1\t$spe2\t$count1\t$count2\t";
	}

	sub posi21{
		%h1=();%h2=(); %ha_num=(); %hb_num=(); %overlap1=();%overlap2=();
		my $spe1=$_[0];	my $spe2=$_[1];
		open B1,"<$opts{i1}/$spe1\_simplified.gff";
			while ($a=<B1>){
			chomp $a;
			my @items1=split/\t/,$a;
			$h1{$items1[1]}=$items1[0]."\t".$items1[2]."\t".$items1[3];
		}
		close B1;
		open B2,"<$opts{i1}/$spe2\_simplified.gff";
		while ($a2=<B2>){
			chomp $a2;
			my @items2=split/\t/,$a2;
			$h2{$items2[1]}=$items2[0]."\t".$items2[2]."\t".$items2[3];
		}
		close B2;
		open O,">$opts{o}/microsyn_genes.links/$spe1\_$spe2\_microsyn_genes.links";
		open I,"<$opts{o}/collinearity/$spe1\_$spe2.collinearity";
		while ($b=<I>){
			chomp $b;
			if ($b!~/^#/) {
				my @items3=split/\t/,$b;
				if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
					$ha_num{$items3[1]}=0; $hb_num{$items3[2]}=0;
					if (exists $h_gene1{$items3[1]} and exists $h_gene2{$items3[2]}) {
						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ca}\n";
					}
					else{

						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ci}\n";
					}
					if (exists $h_gene1{$items3[1]}) {
						$overlap1{$items3[1]}=0;
					}
					if (exists $h_gene2{$items3[2]}) {
						$overlap2{$items3[2]}=0;
					}
				}
				if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
					$ha_num{$items3[2]}=0; $hb_num{$items3[1]}=0;
					if (exists $h_gene2{$items3[1]} and exists $h_gene1{$items3[2]}) {
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ca}\n";
					}
					else{
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ci}\n";
					}
					if (exists $h_gene1{$items3[2]}) {
						$overlap1{$items3[2]}=0;
					}
					if (exists $h_gene2{$items3[1]}) {
						$overlap2{$items3[1]}=0;
					}
				}
			}
		}
		close I;
		close O;
		$count3 = scalar keys %ha_num;
		$count4 = scalar keys %hb_num;
		$overlap_n1 = scalar keys %overlap1;
		$overlap_n2 = scalar keys %overlap2;
		print S "$count3\t$count4\t$overlap_n1\t$overlap_n2\n";
	}
	#micro_highlight
	open I,"<$opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
	open O,">$opts{o}/microsyn_genes.links/1";
	while (my $a=<I>){
		chomp $a;
		my @yu=split/\t/,$a;
		print O "$yu[0]\t$yu[1]\t$yu[2]\n";
		print O "$yu[3]\t$yu[4]\t$yu[5]\n";
	}
	close I;close O;
	system "sort -k 1,1 -k 2n,2 $opts{o}/microsyn_genes.links/1| uniq >$opts{o}/microsyn_genes.links/micro_highlight.txt";
	system "rm $opts{o}/microsyn_genes.links/1";

	#macro_highlight
	open O,">$opts{o}/microsyn_genes.links/macro_highlight.txt";
	my $dir1=$opts{i1};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_diagonal_gene_pairs.pep$/ ) {
			$name=$1;
			my $file=$dir1."/$c";
			open I,"<$file";
			open I2,"<$opts{i1}/$name\_simplified.gff";
			while (my $a=<I2>){
				chomp $a;
				my @ui=split/\t/,$a;
				$gene_chr_p{$ui[1]}=$ui[0]."\t".$ui[2]."\t".$ui[3];
			}
			while (my $a=<I>){
				chomp $a;
				if ($a=~/>(\S+)/) {
					if (exists $gene_chr_p{$1}) {
						print O "$gene_chr_p{$1}\n";
					}
				}
			}
			close I;
			close I2;
		}
	}
	close P;
	close O;
	close S;
}
if ($opts{a} eq "diamond") {
	$dir1=$opts{i1};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ($c!~/_/) {
			if ( $c=~/^(\S+).pep$/ ) {
				my $file=$dir1."/$c";

				my $dirdb = "$opts{i1}/$1.pep.database.dmnd";
				if (-e $dirdb) {
				}
				else{
					system "diamond makedb --in $file -d $file.database";
				}
				$h{$n}=$1;
				$n=$n+1;
			}
		}
	}
	close P;

	$num=keys%h;
	system "mkdir $opts{o}/microsyn_genes.links";
	print "The [$opts{o}/microsyn_genes.links] directory has been created!\n";
	system "mkdir $opts{o}/collinearity";
	print "The [$opts{o}/collinearity] directory has been created!\n";
	open S,">$opts{o}/microsyn_genes.links/summary.txt";
	print S "SpeciesA\tSpeciesB\tA_Macrosyn_genes_numbers\tB_Macrosyn_genes_numbers\tA_Microsyn_genes_numbers\tB_Microsyn_genes_numbers\tA_Overlap_genes_numbers\tB_Overlap_genes_numbers\n";
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "mkdir $opts{o}/$h{$i}_$h{$n1}";
				print "The [$opts{o}/$h{$i}_$h{$n1}] directory has been created!\n";
				system "diamond blastp --query $opts{i1}/$h{$i}.pep --db $opts{i1}/$h{$n1}.pep.database --outfmt 6 --out $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs 5";
				system "cat $opts{i1}/$h{$i}_simplified.gff $opts{i1}/$h{$n1}_simplified.gff > $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.gff";
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
				system "cp $opts{o}/$h{$i}_$h{$n1}/$h{$i}_$h{$n1}.collinearity $opts{o}/collinearity";
				color($h{$i},$h{$n1}); 
				posi2($h{$i},$h{$n1}); 
			}
		}
	}
	system "cat $opts{o}/microsyn_genes.links/*_microsyn_genes.links > $opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
	sub color{
		my $spe1=$_[0];	my $spe2=$_[1];%h_gene1=();%h_gene2=();
		open GENE1,"<$opts{i1}/$spe1\_diagonal_gene_pairs.pep";
		while (my $a=<GENE1>){
			chomp $a;
			if ($a=~/>(\S+)/) {
				$h_gene1{$1}=1;
			}
		}
		close GENE1;
		open GENE2,"<$opts{i1}/$spe2\_diagonal_gene_pairs.pep";
		while (my $b=<GENE2>){
			chomp $b;
			if ($b=~/>(\S+)/) {
				$h_gene2{$1}=1;
			}
		}
		close GENE2;
		$count1 = scalar keys %h_gene1;
		$count2 = scalar keys %h_gene2;
		print S "$spe1\t$spe2\t$count1\t$count2\t";
	}

	sub posi2{
		%h1=();%h2=(); $overlap_n1=0; $overlap_n2=0; my %overlap1=(); my %overlap2=();
		my $spe1=$_[0];	my $spe2=$_[1];
		open B1,"<$opts{i1}/$spe1\_simplified.gff";
			while ($a=<B1>){
			chomp $a;
			my @items1=split/\t/,$a;
			$h1{$items1[1]}=$items1[0]."\t".$items1[2]."\t".$items1[3];
		}
		close B1;
		open B2,"<$opts{i1}/$spe2\_simplified.gff";
		while ($a2=<B2>){
			chomp $a2;
			my @items2=split/\t/,$a2;
			$h2{$items2[1]}=$items2[0]."\t".$items2[2]."\t".$items2[3];
		}
		close B2;
		open O,">$opts{o}/microsyn_genes.links/$spe1\_$spe2\_microsyn_genes.links";
		open I,"<$opts{o}/collinearity/$spe1\_$spe2.collinearity";
		while ($b=<I>){
			chomp $b;
			if ($b!~/^#/) {
				my @items3=split/\t/,$b;
				if (exists $h1{$items3[1]} and exists $h2{$items3[2]}) {
					$ha_num{$items3[1]}=0; $hb_num{$items3[2]}=0;
					if (exists $h_gene1{$items3[1]} and exists $h_gene2{$items3[2]}) {
						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ca}\n";
					}
					else{
						print O "$h1{$items3[1]}\t$h2{$items3[2]}\tcolor=$opts{ci}\n";
					}
					if (exists $h_gene1{$items3[1]}) {
						$overlap1{$items3[1]}=0;
					}
					if (exists $h_gene2{$items3[2]}) {
						$overlap2{$items3[2]}=0;
					}
				}
				if (exists $h2{$items3[1]} and exists $h1{$items3[2]}) {
					$ha_num{$items3[2]}=0; $hb_num{$items3[1]}=0;
					if (exists $h_gene2{$items3[1]} and exists $h_gene1{$items3[2]}) {
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ca}\n";
					}
					else{
						print O "$h2{$items3[1]}\t$h1{$items3[2]}\tcolor=$opts{ci}\n";
					}
					if (exists $h_gene1{$items3[2]}) {
						$overlap1{$items3[2]}=0;
					}
					if (exists $h_gene2{$items3[1]}) {
						$overlap2{$items3[1]}=0;
					}
				}
			}
		}
		close I;
		close O;
		$count3 = scalar keys %ha_num;
		$count4 = scalar keys %hb_num;
		$overlap_n1 = scalar keys %overlap1;
		$overlap_n2 = scalar keys %overlap2;
		print S "$count3\t$count4\t$overlap_n1\t$overlap_n2\n";
	}
	#micro_highlight
	open I,"<$opts{o}/microsyn_genes.links/All_species_microsyn_genes.links";
	open O,">$opts{o}/microsyn_genes.links/1";
	while (my $a=<I>){
		chomp $a;
		my @yu=split/\t/,$a;
		print O "$yu[0]\t$yu[1]\t$yu[2]\n";
		print O "$yu[3]\t$yu[4]\t$yu[5]\n";
	}
	close I;close O;
	system "sort -k 1,1 -k 2n,2 $opts{o}/microsyn_genes.links/1| uniq >$opts{o}/microsyn_genes.links/micro_highlight.txt";
	system "rm $opts{o}/microsyn_genes.links/1";

	#macro_highlight
	open O,">$opts{o}/microsyn_genes.links/macro_highlight.txt";
	$dir1=$opts{i1};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_diagonal_gene_pairs.pep$/ ) {
			$name=$1;
			my $file=$dir1."/$c";
			open I,"<$file";
			open I2,"<$opts{i1}/$name\_simplified.gff";
			while (my $a=<I2>){
				chomp $a;
				my @ui=split/\t/,$a;
				$gene_chr_p{$ui[1]}=$ui[0]."\t".$ui[2]."\t".$ui[3];
			}
			while (my $a=<I>){
				chomp $a;
				if ($a=~/>(\S+)/) {
					if (exists $gene_chr_p{$1}) {
						print O "$gene_chr_p{$1}\n";
					}
				}
			}
			close I;
			close I2;
		}
	}
	close P;
	close O;
	close S;
}
#!/usr/bin/perl -w
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","m=s","n=s","gff=s","g=s","c=s","s=s","o=s","width=s","height=s","h|help");
if (!(defined $opts{i} and defined $opts{m} and defined $opts{n} and defined $opts{gff} and defined $opts{c} and defined $opts{g} and defined $opts{s} and defined $opts{o})) {
		die "************************************************\n
	-i	Full path to the folder containing required data files (same as [-i] in the previous script [Macrosyn1_1_blast.pl])
	-m	Selected Ancestor'name by abbreviations (example: NVec)
	-n	Interested species'name by abbreviations (example: HSap)
	-gff Full path to the interested species' s gff file
	-g	Full path to the interested species' s genome sequence file
	-c	Full path to the [chr_id.txt] file
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the folder storing the output results (same as [-o] in the previous script [Macrosyn3.pl] or [Macrosyn_add_spe2-blast.pl])
	Optional:
	-width	Sets the width of the output graph (default: 10)
	-height	Sets the height of the output graph (default: 6)
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the folder containing required data files (same as [-i] in the previous script [Macrosyn1_1_blast.pl])
	-m	Selected Ancestor'name by abbreviations (example: NVec)
	-n	Interested species'name by abbreviations (example: HSap)
	-gff Full path to the interested species' s gff file
	-g	Full path to the interested species' s genome sequence file
	-c	Full path to the [chr_id.txt] file
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the folder storing the output results (same as [-o] in the previous script [Macrosyn3.pl] or [Macrosyn_add_spe2-blast.pl])
	Optional:
	-width	Sets the width of the output graph (default: 10)
	-height	Sets the height of the output graph (default: 6)
		*************************************************\n";
}

if (!defined $opts{width}) {
	$opts{width}=10;
}
if (!defined $opts{width}) {
	$opts{width}=6;
}


open I,"<$opts{c}" or die("Could not open $opts{c}.\n"); 
my %h_chrnum;
while (my $a=<I>) {
	chomp $a;
	if ($a!~/chr_number/) {
		my @it=split/\t/,$a;
		$h_chrnum{$it[0]}=$it[1];
	}
}
close I;


my $spename1=$opts{m};
my $spename2=$opts{n};
my $block_name="Ancestors_database/".$spename1."/$spename1.block";
my $gff_an="Ancestors_database/".$spename1."/$spename1.gff";
my $pathway_block="$opts{i}/$block_name";
open I,"<$pathway_block" or die("Could not open $pathway_block.\n"); 
$chr_num_max=0;
while (my $a=<I>) {
	   chomp $a;
	    my @chr_num=split/\t/,$a;
		if ($chr_num[0]>$chr_num_max) {
			$chr_num_max=$chr_num[0];
		}
}
close I;
chromosome($opts{m},$opts{n});
sub chromosome{
	$spe1=$_[0];$spe2=$_[1];$zero=0;
	system "mkdir $opts{o}/$spe1\_$spe2/chromosome_comparison";
	open I,"<$opts{gff}" or die("Could not open $opts{gff}.\n"); 
	while (my $a=<I>){
		chomp $a;
		@items1=split/\t/,$a;
		$h1{$items1[0]}=0;
	}
	close I;

	open I,"<$opts{g}";
	while(<I>){
		chomp;
		if (/>(\S+)/) { $id=$1; }
		else { $h2{$id} .= $_; } 
	}
	close I;

	foreach $k (sort keys %h2){
		$len=length($h2{$k});
		$h3{$k}=$len;
	}

	open O,"> $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt";
	#print O "Chr\tStart\tEnd\n";
	foreach my $i (sort keys %h3) {
		if (exists $h1{$i}) {
			print O "$i\t$zero\t$h3{$i}\n";
		}
	}
	close O;
	system "sort -Vk 1,1 -k 2n,2 $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt >$opts{o}/$spe1\_$spe2/chromosome_comparison/1";
	system "rm $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt";
	open I,"<$opts{o}/$spe1\_$spe2/chromosome_comparison/1" or die("$opts{o}/$spe1\_$spe2/chromosome_comparison/1.\n"); 
	open O,">$opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt" or die("$opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt.\n"); 
	print O "Chr\tStart\tEnd\n";
	while (my $a=<I>){
		chomp $a;
		my @it=split/\t/,$a;
		if (exists $h_chrnum{$it[0]}) {
			print O "$h_chrnum{$it[0]}\t$it[1]\t$it[2]\n";
		}
	}
	close I;
	close O;
	system "rm $opts{o}/$spe1\_$spe2/chromosome_comparison/1";

	open I,"<$opts{i}/$gff_an" or die("Could not open $opts{i}/$gff_an\n"); 
	while (my $a=<I>){
		chomp $a;
		@items2=split/\t/,$a;
		$key=$items2[1];
		open M,"<$pathway_block";
		while (my $b=<M>) {
			@it=split/\t/,$b;
			if ($items2[0] eq $it[1]) {
				if ($items2[2]>=$it[2] and $items2[3]<=$it[3] ) {
					$h4_py_gene_pos{$key}=$it[0];
				}
			}
		}
		close M;
	}
	close I;
	open I,"<$opts{gff}";
	while (my $a=<I>){
		chomp $a;
		@items3=split/\t/,$a;
		$h5_bp_gene_pos{$items3[1]}=$items3[0]."\t".$items3[2]."\t".$items3[3];
	}
	close I;

	open I,"<$opts{o}/$opts{m}_$opts{n}/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
	open I2,"<$opts{o}/$opts{m}_$opts{n}/Macrosyn_genes/$opts{m}_$opts{n}_significant_gene_pairs.dot";
	open O,"> $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_density.txt";
	open O2,"> $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_density_sig.txt";
	print O "Chr\tStart\tEnd\tAncestor_chr\n";
	while (my $a=<I>){
		if ($a!~/^ALG/) {
			chomp $a;
			@items4=split/\t/,$a;
			if (exists $h4_py_gene_pos{$items4[0]}) {
				if (exists $h5_bp_gene_pos{$items4[1]}) {
					@items5=split/\t/,$h5_bp_gene_pos{$items4[1]};
					if (exists $h_chrnum{$items5[0]}) {
					print O "$h_chrnum{$items5[0]}\t$items5[1]\t$items5[2]\t$h4_py_gene_pos{$items4[0]}\n";
					}
				}
			}
		}
	}
	close I;
	close O;
	print O2 "Chr\tStart\tEnd\tAncestor_chr\n";
	while (my $a=<I2>){
		if ($a!~/^ALG/) {
			chomp $a;
			@items4=split/\t/,$a;
			if (exists $h4_py_gene_pos{$items4[0]}) {
				if (exists $h5_bp_gene_pos{$items4[1]}) {
					@items5=split/\t/,$h5_bp_gene_pos{$items4[1]};
					if (exists $h_chrnum{$items5[0]}) {
					print O2 "$h_chrnum{$items5[0]}\t$items5[1]\t$items5[2]\t$h4_py_gene_pos{$items4[0]}\n";
					}
				}
			}
		}
	}
	close I2;
	close O2;


	open I,"<$pathway_block";
	open O,"> $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_density.txt";
	print O "Chr\tStart\tEnd\tAncestor_chr\n";
	while (my $a=<I>){
			chomp $a;
			@items_b=split/\t/,$a;
			$aaa=0;
			$chr=$items_b[0];
			if (exists $h_chr{$chr}) {
				$h_chr{$chr}=$h_chr{$chr}+$items_b[3]-$items_b[2];
			}
			else{
				$h_chr{$chr}=$items_b[3]-$items_b[2];
			}
	}
	foreach my $iii (sort {$a<=>$b} keys %h_chr) {
		$aaa=0;
		$hhhh=$h_chr{$iii}+1;
		print O "$iii\t$aaa\t$hhhh\t$iii\n";
	}		

	close I;
	close O;

	open O,"> $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt";
	#print O "Chr\tStart\tEnd\n";

	foreach my $ii (sort {$a<=>$b} keys %h_chr) {
		$aa=0;
		$hhh=$h_chr{$ii}+1;
		print O "$ii\t$aa\t$hhh\n";
	}		

	close I;
	close O;
	system "sort -Vk1,1 -k 2n,2 $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt >$opts{o}/$spe1\_$spe2/chromosome_comparison/1";
	system "rm $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt";
	open I,"<$opts{o}/$spe1\_$spe2/chromosome_comparison/1" or die("$opts{o}/$spe1\_$spe2/chromosome_comparison/1.\n"); 
	open O,">$opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt" or die("$opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt.\n"); 
	print O "Chr\tStart\tEnd\n";
	while (my $a=<I>){
		print O "$a";
	}
	close I;
	close O;
	system "rm $opts{o}/$spe1\_$spe2/chromosome_comparison/1";
	system("Rscript $opts{s}/Chr_comparison.R $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_density.txt $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_karyotype.txt $opts{o}/$spe1\_$spe2/chromosome_comparison $spe1  $chr_num_max $opts{width} $opts{height}");
	system("Rscript $opts{s}/Chr_comparison.R $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_density.txt $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt $opts{o}/$spe1\_$spe2/chromosome_comparison $spe1\_$spe2  $chr_num_max $opts{width} $opts{height}");
	system("Rscript $opts{s}/Chr_comparison.R $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_density_sig.txt $opts{o}/$spe1\_$spe2/chromosome_comparison/$spe1\_$spe2\_karyotype.txt $opts{o}/$spe1\_$spe2/chromosome_comparison $spe1\_$spe2\_sig  $chr_num_max $opts{width} $opts{height}");
	
}

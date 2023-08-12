#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"pansyn=s","i1=s","m=s","color=s","n=s","gff=s","g=s","c=s","ab=s","width=s","height=s","h|help");
if (!(defined $opts{i1} and defined $opts{m} and defined $opts{n} and defined $opts{gff} and defined $opts{color} and defined $opts{c} and defined $opts{g} and defined $opts{ab} and defined $opts{pansyn})) {
		die "************************************************\n
	-i1	Full path to the [outputDir1_S17] directory
	-m	The abbreviation for the name of the species representing the ancestral genome (e.g. NVec)
	-n	The abbreviation for the name of the interested species (e.g. HSap)
	-gff	Full path to the gene coordinates file of the interested species (e.g. HSap_simplified.gff)
	-g	Full path to the genome file of the interested species (e.g. HSap.fa)
	-c	Full path to the [*_chr_id.txt] file
	-color	Full path to the [*_ancestor_chr_color.txt] file
	-ab	 Full path to the [speAname_speBname] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
	-width	The width of the output graph (default: 10)
	-height	The height of the output graph (default: 6)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i1	Full path to the [outputDir1_S17] directory
	-m	The abbreviation for the name of the species representing the ancestral genome (e.g. NVec)
	-n	The abbreviation for the name of the interested species (e.g. HSap)
	-gff	Full path to the gene coordinates file of the interested species (e.g. HSap_simplified.gff)
	-g	Full path to the genome file of the interested species (e.g. HSap.fa)
	-c	Full path to the [*_chr_id.txt] file
	-color	Full path to the [*_ancestor_chr_color.txt] file
	-ab	 Full path to the [speAname_speBname] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
	-width	The width of the output graph (default: 10)
	-height	The height of the output graph (default: 6)
	-h|-help	Print this help page
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
########
if ($opts{ab} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{ab} =~ s/$slash$//;
}
########
####
if ($opts{pansyn} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{pansyn} =~ s/$slash$//;
}
####

if (!defined $opts{width}) {
	$opts{width}=10;
}
if (!defined $opts{height}) {
	$opts{height}=6;
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
my $gff_an="Ancestors_database/".$spename1."/$spename1\_simplified.gff";
my $pathway_block="$opts{i1}/$block_name";
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

	my $dir = "$opts{ab}/chromosome_comparison";
	if (-d $dir) {
		system "rm -r $dir";
		print "The $opts{ab}/chromosome_comparison already exists!\nThe old directory has been deleted.\n";
		system "mkdir $opts{ab}/chromosome_comparison";
		print "The [$opts{ab}/chromosome_comparison] directory has been created!\n";

	} else {
		system "mkdir $opts{ab}/chromosome_comparison";
		print "The [$opts{ab}/chromosome_comparison] directory has been created!\n";
	}

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

	open O,"> $opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt";
	#print O "Chr\tStart\tEnd\n";
	foreach my $i (sort keys %h3) {
		if (exists $h1{$i}) {
			print O "$i\t$zero\t$h3{$i}\n";
		}
	}
	close O;
	system "sort -Vk 1,1 -k 2n,2 $opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt >$opts{ab}/chromosome_comparison/1";
	system "rm $opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt";
	open I,"<$opts{ab}/chromosome_comparison/1" or die("$opts{ab}/chromosome_comparison/1.\n"); 
	open O,">$opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt" or die("$opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt.\n"); 
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
	system "rm $opts{ab}/chromosome_comparison/1";

	open I,"<$opts{i1}/$gff_an" or die("Could not open $opts{i1}/$gff_an\n"); 
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

	open I,"<$opts{ab}/Macrosyn_genes/$opts{m}_$opts{n}_all_gene_pairs.dot";
	open I2,"<$opts{ab}/Macrosyn_genes/$opts{m}_$opts{n}_significant_gene_pairs.dot";
	open O,"> $opts{ab}/chromosome_comparison/$spe1\_$spe2\_density.txt";
	open O2,"> $opts{ab}/chromosome_comparison/$spe1\_$spe2\_density_sig.txt";
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
	open O,"> $opts{ab}/chromosome_comparison/$spe1\_density.txt";
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

	open O,"> $opts{ab}/chromosome_comparison/$spe1\_karyotype.txt";
	#print O "Chr\tStart\tEnd\n";

	foreach my $ii (sort {$a<=>$b} keys %h_chr) {
		$aa=0;
		$hhh=$h_chr{$ii}+1;
		print O "$ii\t$aa\t$hhh\n";
	}		

	close I;
	close O;
	system "sort -Vk1,1 -k 2n,2 $opts{ab}/chromosome_comparison/$spe1\_karyotype.txt >$opts{ab}/chromosome_comparison/1";
	system "rm $opts{ab}/chromosome_comparison/$spe1\_karyotype.txt";
	open I,"<$opts{ab}/chromosome_comparison/1" or die("$opts{ab}/chromosome_comparison/1.\n"); 
	open O,">$opts{ab}/chromosome_comparison/$spe1\_karyotype.txt" or die("$opts{ab}/chromosome_comparison/$spe1\_karyotype.txt.\n"); 
	print O "Chr\tStart\tEnd\n";
	while (my $a=<I>){
		print O "$a";
	}
	close I;
	close O;
	system "rm $opts{ab}/chromosome_comparison/1";
	system("Rscript $opts{pansyn}/Chr_comparison.R $opts{ab}/chromosome_comparison/$spe1\_density.txt $opts{ab}/chromosome_comparison/$spe1\_karyotype.txt $opts{ab}/chromosome_comparison $spe1  $chr_num_max $opts{width} $opts{height} $opts{color}");
	system("Rscript $opts{pansyn}/Chr_comparison.R $opts{ab}/chromosome_comparison/$spe1\_$spe2\_density.txt $opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt $opts{ab}/chromosome_comparison $spe1\_$spe2 $chr_num_max $opts{width} $opts{height} $opts{color}");
	system("Rscript $opts{pansyn}/Chr_comparison.R $opts{ab}/chromosome_comparison/$spe1\_$spe2\_density_sig.txt $opts{ab}/chromosome_comparison/$spe1\_$spe2\_karyotype.txt $opts{ab}/chromosome_comparison $spe1\_$spe2\_sig $chr_num_max $opts{width} $opts{height} $opts{color}");
	
}

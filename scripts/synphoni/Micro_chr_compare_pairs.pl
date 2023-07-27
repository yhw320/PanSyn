#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","o=s","color=s","r=s","w=s","sl=s","m=s","h|help");
if (!( defined $opts{i} and defined $opts{o} and defined $opts{color} and defined $opts{r})) {
		die "************************************************\n
	-i	Full path to the [inputDir] directoty containing input files
	-r	Enter the abbreviated name of the reference spe-cies (example: HSap)
	-o	Full path to the [outputDir] directoty containing output files
	-color	Full path to the [*.color] file containing the color information for chromosomes
	Optional:
	-w	Window size (default:10)
	-sl	Slide size (default:5)
	-m	Max distance (bp) between orthologs (default:50000000)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [inputDir] directoty containing input files
	-r	Enter the abbreviated name of the reference spe-cies (example: HSap)
	-o	Full path to the [outputDir] directoty containing output files
	-color	Full path to the [*.color] file containing the color information for chromosomes
	Optional:
	-w	Window size (default:10)
	-sl	Slide size (default:5)
	-m	Max distance (bp) between orthologs (default:50000000)
	-h|-help	Print this help page
		*************************************************\n";
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

if (!(defined $opts{w})) {
	$opts{w}=10;
}
if (!(defined $opts{sl})) {
	$opts{sl}=5;
}
if (!(defined $opts{m})) {
	$opts{m}=50000000;
}
######################
my $dir = "$opts{o}"; 

opendir my $dh, $dir or die "can't open $dir: $!";

my @files = readdir($dh);       
closedir($dh);   

if (@files == 2) {
} 
else {
	system "rm -r $opts{o}/*";
}

#######################

my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^$opts{r}_(\S+).pairs$/ ) {
		my $mm_name=$1;
		my %h_hs_gene_pos=();my %h_mm_gene_pos=();
		open I1,"< $opts{i}/$opts{r}_simplified_chr.gff";
		while (my $a=<I1>) {
			chomp $a;
			my @it=split/\t/,$a;
			#if ($it[4]=~/"+"/) {
				my $zhi=1;
			#	print "$zhi\n";
			#}
			#else{my $zhi=0;}
			$h_hs_gene_pos{$it[1]}=$it[0]."\t".$it[2]."\t".$it[3]."\t".$zhi;
		}
		close I1;
		open I1,"< $opts{i}/$mm_name\_simplified_chr.gff";
		while (my $a=<I1>) {
			chomp $a;
			my @it=split/\t/,$a;
			#if ($it[4]=~/"+"/) {
				my $zhi=1;
			#}
			#else{my $zhi=0;}
			$h_mm_gene_pos{$it[1]}=$it[0]."\t".$it[2]."\t".$it[3]."\t".$zhi;
		}
		close I1;
		
		open O,"> $opts{o}/$opts{r}_$mm_name.lab";
		open I,"< $opts{i}/$c";
		my $n=1;my %h_n=();
		while (my $a=<I>) {
			chomp $a;
			my @it=split/\t/,$a;
			if (exists $h_hs_gene_pos{$it[0]}) {
				my $pos=$h_hs_gene_pos{$it[0]};
				my @it2=split/\t/,$pos;
				print O "$n\t$it2[0]\t$a\n";
				$h_n{$it[0]}{$it[1]}=$n;
				$n=$n+1;
			}
		}
		close I;
		close O;

		open I,"< $opts{color}";
		while (my $a=<I>) {
			chomp $a;
			my @ha=split/\t/,$a;
			$chr_color{$ha[0]}=$ha[1];
		}
		close I;
			
		open O,"> $opts{o}/$opts{r}_$mm_name.msynt";
		open I,"< $opts{i}/$c";
		while (my $a=<I>) {
			chomp $a;
			my @it=split/\t/,$a;
			if (exists $h_hs_gene_pos{$it[0]}) {
				my @it1=split/\t/,$h_hs_gene_pos{$it[0]};
				if (exists $chr_color{$it1[0]}) {
					print O "$it1[0]\t$it[0]\t$it1[1]\t$h_n{$it[0]}{$it[1]}\t$it1[3]\t$chr_color{$it1[0]}\t$it1[0]\n";
				}
			}
		}
		close O;

		open O,"> $opts{o}/$mm_name\_$opts{r}.msynt";
		open I,"< $opts{i}/$c";
		while (my $a=<I>) {
			chomp $a;
			my @it=split/\t/,$a;
			if (exists $h_mm_gene_pos{$it[1]} and exists $h_hs_gene_pos{$it[0]}) {
				my @it1=split/\t/,$h_mm_gene_pos{$it[1]};
				my @it2=split/\t/,$h_hs_gene_pos{$it[0]};
				if (exists $chr_color{$it2[0]}) {
						print O "$it1[0]\t$it[1]\t$it1[1]\t$h_n{$it[0]}{$it[1]}\t$it1[3]\t$chr_color{$it2[0]}\t$it1[0]\n";
				}
			}
		}
		close O;
		system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}_$mm_name.msynt >$opts{o}/$opts{r}_$mm_name.sorted.msynt";
		system "sort -k 1,1 -k 3n,3 $opts{o}/$mm_name\_$opts{r}.msynt >$opts{o}/$mm_name\_$opts{r}.sorted.msynt";
		system "perl drawCLGContrib2.4_2.pl $opts{o}/$opts{r}_$mm_name.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$mm_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$mm_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} $opts{o}/$mm_name\_$opts{r}.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$mm_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$mm_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}_$mm_name.svg";
	}
}
###############################3
open I,"< $opts{i}/$opts{r}_simplified_chr.gff";
my %chr_min=();
while (my $a=<I>) {
	chomp $a;
	my @aa=split/\t/,$a;
	if (exists $chr_min{$aa[0]}) {
		if ($aa[2]<$chr_min{$aa[0]}) {
			$chr_min{$aa[0]}=$aa[2];
		}
	}
	else{$chr_min{$aa[0]}=$aa[2];}

}
close I;

open O,"> $opts{o}/$opts{r}_add.gff";
open I,"< $opts{i}/$opts{r}_simplified_chr.gff";
while (my $a=<I>) {
	chomp $a;
	my @aa=split/\t/,$a;
	print O "$aa[0]\t$aa[1]\t$aa[2]\n";
}
close I;
my $gene_number;

foreach my $i (keys %chr_min) {
	$gene_number=1;
	for ($g_start=0;$g_start<$chr_min{$i};$g_start=$g_start+20000) {
		$add_gene_name=$i."_"."gene".$gene_number;
		print O "$i\t$add_gene_name\t$g_start\n";
		$gene_number=$gene_number+1;
	}
}

close O;
close I;

system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}_add.gff >$opts{o}/$opts{r}_add.sorted.gff";

			open I,"< $opts{color}";
			while (my $a=<I>) {
				chomp $a;
				my @ha=split/\t/,$a;
				$chr_color{$ha[0]}=$ha[1];
			}
			close I;
my %hhh=();
open O,"> $opts{o}/$opts{r}_final.lab";
open O1,"> $opts{o}/$opts{r}_final.msynt";
open I,"< $opts{o}/$opts{r}_add.sorted.gff";
while (my $a=<I>) {
	chomp $a;
	my @aa=split/\t/,$a;
	if (!exists $hhh{$aa[0]}) {
		$jishu=0;
		$hhh{$aa[0]}=0;
	}
	if (exists $hhh{$aa[0]}) {
		$jishu=$jishu+1;
	}	
	print O1 "$aa[0]\t$aa[1]\t$aa[2]\t$jishu\t1\t$chr_color{$aa[0]}\t$aa[0]\n";
	print O "$jishu\t$aa[0]\t$aa[1]\n";
}
close I;
close O;close O1;
system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}_final.msynt >$opts{o}/$opts{r}_final.msynt.sorted";
system "sort -k 1n,1 $opts{o}/$opts{r}_final.lab >$opts{o}/$opts{r}_final.lab.sorted";

system "perl drawCLGContrib2.4_2.pl $opts{o}/$opts{r}_final.msynt.sorted:,algcolor=$opts{o}/$opts{r}_final.msynt.sorted:type=alg,file=$opts{o}/$opts{r}_final.lab.sorted,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}_reference.svg";





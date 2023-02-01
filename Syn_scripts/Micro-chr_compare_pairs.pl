#!/usr/bin/perl -w
#use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","o=s","color=s","a=s","s=s","w=s","sl=s","m=s","h|help");
if (!( defined $opts{i} and defined $opts{o} and defined $opts{color} and defined $opts{a} and defined $opts{s})) {
		die "************************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-a	Reference species by name abbreviations (example: HSap)
	-o	Full path to the [outputDir] folder containing output files
	-color Full path to the chromosomal color file
	-s	Full path to the [syn_scripts] file PanSyn provided
	Optional:
	-w	Window size (default:10)
	-sl	Slide size (default:5)
	-m	Max distance (bp) between orthologs (default:50000000)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-a	Reference species by name abbreviations (example: HSap)
	-o	Full path to the [outputDir] folder containing output files
	-color Full path to the chromosomal color file
	-s	Full path to the [syn_scripts] file PanSyn provided
	Optional:
	-w	Window size (default:10)
	-sl	Slide size (default:5)
	-m	Max distance (bp) between orthologs (default:50000000)
	-h|-help	Print this help page
		*************************************************\n";
}
if (!(defined $opts{w})) {
	$opts{w}=10;
}
if (!(defined $opts{sl})) {
	$opts{sl}=5;
}
if (!(defined $opts{m})) {
	$opts{m}=50000000;
}

my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^$opts{a}_(\S+).pairs$/ ) {
		my $mm_name=$1;
		my %h_hs_gene_pos=();my %h_mm_gene_pos=();
		open I1,"< $opts{i}/$opts{a}.gff";
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
		open I1,"< $opts{i}/$mm_name.gff";
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
		
		open O,"> $opts{o}/$opts{a}_$mm_name.lab";
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
			
		open O,"> $opts{o}/$opts{a}_$mm_name.msynt";
		open I,"< $opts{i}/$c";
		while (my $a=<I>) {
			chomp $a;
			my @it=split/\t/,$a;
			if (exists $h_hs_gene_pos{$it[0]}) {
				my @it1=split/\t/,$h_hs_gene_pos{$it[0]};
				print O "$it1[0]\t$it[0]\t$it1[1]\t$h_n{$it[0]}{$it[1]}\t$it1[3]\t$chr_color{$it1[0]}\t$it1[0]\n";
			}
		}
		close O;

		open O,"> $opts{o}/$mm_name\_$opts{a}.msynt";
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
		system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{a}_$mm_name.msynt >$opts{o}/$opts{a}_$mm_name.sorted.msynt";
		system "sort -k 1,1 -k 3n,3 $opts{o}/$mm_name\_$opts{a}.msynt >$opts{o}/$mm_name\_$opts{a}.sorted.msynt";
		system "perl $opts{s}/drawCLGContrib2.4_2.pl $opts{o}/$opts{a}_$mm_name.sorted.msynt:,algcolor=$opts{o}/$opts{a}_$mm_name.sorted.msynt:type=alg,file=$opts{o}/$opts{a}_$mm_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} $opts{o}/$mm_name\_$opts{a}.sorted.msynt:,algcolor=$opts{o}/$opts{a}_$mm_name.sorted.msynt:type=alg,file=$opts{o}/$opts{a}_$mm_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{a}_$mm_name.svg";
	}
}
###############################3
open I,"< $opts{i}/$opts{a}.gff";
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

open O,"> $opts{o}/$opts{a}_add.gff";
open I,"< $opts{i}/$opts{a}.gff";
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

system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{a}_add.gff >$opts{o}/$opts{a}_add.sorted.gff";

			open I,"< $opts{color}";
			while (my $a=<I>) {
				chomp $a;
				my @ha=split/\t/,$a;
				$chr_color{$ha[0]}=$ha[1];
			}
			close I;
my %hhh=();
open O,"> $opts{o}/$opts{a}_final.lab";
open O1,"> $opts{o}/$opts{a}_final.msynt";
open I,"< $opts{o}/$opts{a}_add.sorted.gff";
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
system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{a}_final.msynt >$opts{o}/$opts{a}_final.msynt.sorted";
system "sort -k 1n,1 $opts{o}/$opts{a}_final.lab >$opts{o}/$opts{a}_final.lab.sorted";

system "perl $opts{s}/drawCLGContrib2.4_2.pl $opts{o}/$opts{a}_final.msynt.sorted:,algcolor=$opts{o}/$opts{a}_final.msynt.sorted:type=alg,file=$opts{o}/$opts{a}_final.lab.sorted,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{a}_reference.svg";





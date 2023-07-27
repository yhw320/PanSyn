#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","t=s","e=s","o=s","s=s","a=s","h|help");
if(!(defined $opts{i} and defined $opts{o} and defined $opts{a})){

	die"**********************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-o	Full path to the [outputDir] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
	-s	The number of best non-self protein hits (default:5)
	-h|-help Print this help page
		*********************************************\n";
}
if((defined $opts{h} or defined $opts{help})){

	die"**********************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-o	Full path to the [outputDir] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
	-s	The number of best non-self protein hits (default:5)
	-h|-help Print this help page
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

######################



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
	my $num=keys %h; my $n1;
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "blastp -query $opts{i}/$h{$i}.pep -db $opts{i}/$h{$n1}.pep -outfmt 6 -out $opts{o}/$h{$i}_$h{$n1}_initial.blast -evalue $opts{e} -num_threads $opts{t} -max_target_seqs $opts{s}";
			}
		}
	}

	opendir P,$dir1;
	my %change=();my %chr_max=();my $nn;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_simplified_chr.gff$/ ) {
			my $file=$dir1."/$c";
			system "sort -k 1,1 -k 3n,3 $file >$opts{o}/$1.sorted.gff";
			open I,"< $opts{o}/$1.sorted.gff";
			open O,"> $opts{o}/$1.gff";
			my %h_chr=(); 
			while (my $a=<I>) {
				chomp $a;
				my @it=split/\t/,$a;
				if (!exists $h_chr{$it[0]}) {
					$nn=1;
					$h_chr{$it[0]}=0;
					print O "$it[0]\t$it[0]_$nn\t$it[2]\t$it[3]\t$it[4]\t$nn\t$it[1]\n";
					$change{$it[1]}=$it[0]."_".$nn;
					$chr_max{$it[0]}=$nn;
				}
				else{
					$nn=$nn+1;
					print O "$it[0]\t$it[0]_$nn\t$it[2]\t$it[3]\t$it[4]\t$nn\t$it[1]\n";
					$change{$it[1]}=$it[0]."_".$nn;
					if ($nn>=$chr_max{$it[0]}) {
						$chr_max{$it[0]}=$nn;
					}
				}
			}
			close I;
			system "rm $opts{o}/$1.sorted.gff";
		}
	}
	close P;
	close O;

	opendir P,$opts{o};
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_initial.blast$/ ) {
			my $file=$opts{o}."/$c";
			open I,"< $file";
			open O,"> $opts{o}/$1.blast";
			while (my $a=<I>) {
				chomp $a;
				my @it=split/\t/,$a;
				if (exists $change{$it[0]} and exists $change{$it[1]}) {
					print O "$change{$it[0]}\t$change{$it[1]}\t$it[2]\t$it[3]\t$it[4]\t$it[5]\t$it[6]\t$it[7]\t$it[8]\t$it[9]\t$it[10]\t$it[11]\n";
				}
			}
			close I;
		}
	}
	close P;close O;

	opendir P,$opts{i};
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).len$/ ) {
			my $file=$opts{i}."/$c";
			open I,"< $file";
			open O,"> $opts{o}/$1.len";
			while (my $a=<I>) {
				chomp $a;
				my @iii=split/\t/,$a;
				print O "$iii[0]\t$iii[1]\t$chr_max{$iii[0]}\n";
			}
			close I;
		}
	}
	close P;close O;
}


if ($opts{a} eq "diamond") {
	my $n=1; my %h=();
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
	my $num=keys %h;my $n1;
	foreach my $i (sort {$a<=>$b} keys %h) {
		for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
			if ($i!=$num) {
				system "diamond blastp --query $opts{i}/$h{$i}.pep --db $opts{i}/$h{$n1}.pep.database --outfmt 6 --out $opts{o}/$h{$i}_$h{$n1}_initial.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs $opts{s}";
			}
		}
	}

	opendir P,$dir1;
	my %change=();my %chr_max=();my $nn;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_simplified_chr.gff$/ ) {
			my $file=$dir1."/$c";
			system "sort -k 1,1 -k 3n,3 $file >$opts{o}/$1.sorted.gff";
			open I,"< $opts{o}/$1.sorted.gff";
			open O,"> $opts{o}/$1.gff";
			my %h_chr=();
			while (my $a=<I>) {
				chomp $a;
				my @it=split/\t/,$a;
				if (!exists $h_chr{$it[0]}) {
					$nn=1;
					$h_chr{$it[0]}=0;
					print O "$it[0]\t$it[0]_$nn\t$it[2]\t$it[3]\t$it[4]\t$nn\t$it[1]\n";
					$change{$it[1]}=$it[0]."_".$nn;
					$chr_max{$it[0]}=$nn;
				}
				else{
					$nn=$nn+1;
					print O "$it[0]\t$it[0]_$nn\t$it[2]\t$it[3]\t$it[4]\t$nn\t$it[1]\n";
					$change{$it[1]}=$it[0]."_".$nn;
					if ($nn>=$chr_max{$it[0]}) {
						$chr_max{$it[0]}=$nn;
					}
				}
			}
			close I;
			system "rm $opts{o}/$1.sorted.gff";
		}
	}
	close P;
	close O;

	opendir P,$opts{o};
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_initial.blast$/ ) {
			my $file=$opts{o}."/$c";
			open I,"< $file";
			open O,"> $opts{o}/$1.blast";
			while (my $a=<I>) {
				chomp $a;
				my @it=split/\t/,$a;
				if (exists $change{$it[0]} and exists $change{$it[1]}) {
					print O "$change{$it[0]}\t$change{$it[1]}\t$it[2]\t$it[3]\t$it[4]\t$it[5]\t$it[6]\t$it[7]\t$it[8]\t$it[9]\t$it[10]\t$it[11]\n";
				}
			}
			close I;
		}
	}
	close P;close O;

	opendir P,$opts{i};
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).len$/ ) {
			my $file=$opts{i}."/$c";
			open I,"< $file";
			open O,"> $opts{o}/$1.len";
			while (my $a=<I>) {
				chomp $a;
				my @iii=split/\t/,$a;
				if (exists $chr_max{$iii[0]}) {
					print O "$iii[0]\t$iii[1]\t$chr_max{$iii[0]}\n";
				}
			}
			close I;
		}
	}
	close P;close O;
}
#!/usr/bin/perl-w
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","t=s","e=s","o=s","h|help");
if(!(defined $opts{i} and defined $opts{o})){

	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-s	The number of best non-self BLASTP hits (default:5)
		*********************************************\n";
}
if (!(defined $opts{e})) {
	$opts{e}=1e-3;
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{s})) {
	$opts{s}=5;
}
if((defined $opts{h} or defined $opts{help})){

	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-s	The number of best non-self BLASTP hits (default:5)
		*********************************************\n";
}

$n=1;
my $dir1=$opts{i};
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).pep$/ ) {
		my $file=$dir1."/$c";
		system "diamond makedb --in $file -d $file.database";
		$h{$n}=$1;
		$n=$n+1;
	}
}
close P;
$num=keys %h;
foreach my $i (sort {$a<=>$b} keys %h) {
	for ($n1=$i+1;$n1<=$num;$n1=$n1+1) {
		if ($i!=$num) {
	        system "diamond blastp --query $opts{i}/$h{$i}.pep --db $opts{i}/$h{$n1}.pep.database --outfmt 6 --out $opts{o}/$h{$i}_$h{$n1}_initial.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs $opts{s}";
		}
	}
}

opendir P,$dir1;
my %change=();
while (my $c=readdir P) {
	if ( $c=~/^(\S+).gff$/ ) {
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
			print O "$change{$it[0]}\t$change{$it[1]}\t$it[2]\t$it[3]\t$it[4]\t$it[5]\t$it[6]\t$it[7]\t$it[8]\t$it[9]\t$it[10]\t$it[11]\n";
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

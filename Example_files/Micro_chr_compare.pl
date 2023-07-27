#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","t=s","e=s","o=s","color=s","r=s","w=s","sl=s","m=s","a=s","h|help");
if (!( defined $opts{i} and defined $opts{o} and defined $opts{color} and defined $opts{r} and defined $opts{a})) {
		die "************************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-r	Enter the abbreviated name of the reference species (example: HSap)
	-o	Full path to the [outputDir] directory containing output files
	-color	Full path to the [*.color] file containing the color information for chromosomes
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
	-w	Window size (default:10)
	-sl	Slide size (default:5)
	-m	Max distance (bp) between orthologs (default:50000000)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i	Full path to the [inputDir] directory containing input files
	-r	Enter the abbreviated name of the reference species (example: HSap)
	-o	Full path to the [outputDir] directory containing output files
	-color	Full path to the [*.color] file containing the color information for chromosomes
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-t	Protein alignment threads (default:12)
	-e	Protein alignment evalue (default:1e-5/0.001)
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


if (!(defined $opts{e})) {
	if ($opts{a} eq "blast") {
		$opts{e}=1e-5;
	}
	if ($opts{a} eq "diamond") {
		$opts{e}=0.001;
	}
}
if (!(defined $opts{e})) {
	$opts{e}=1e-5;
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

######################
my $py_pep; my $i2; my $n; my %chr_color=(); my %h_gene_id=();my $cf; my $g_start; my $jishu; my $add_gene_name;
if ($opts{a} eq "blast") {
	my $dir1=$opts{i}; 
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $pep_file=$dir1."/$c";
			my $dirphr = "$opts{i}/$1.pep.phr";
			my $dirpin = "$opts{i}/$1.pep.pin";
			my $dirpsq = "$opts{i}/$1.pep.psq";
			if (-e $dirphr and -e $dirpin and -e $dirpsq) {
			}
			else{
				system "makeblastdb -in $pep_file -dbtype prot";
			}
			if ($1 eq $opts{r}) {
				$py_pep=$dir1."/$c";
			}
		}
	}
	close P;

	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $cf_name=$1;
			if ($cf_name ne $opts{r}) {
				my $pep_file=$dir1."/$c";
				system "blastp -query $py_pep -db $pep_file -outfmt 6 -out $opts{o}/$opts{r}_$cf_name.blastp -evalue $opts{e} -num_threads $opts{t}";
				system "blastp -query $pep_file -db $py_pep -outfmt 6 -out $opts{o}/$cf_name\_$opts{r}.blastp -evalue $opts{e} -num_threads $opts{t}";
				open I,"< $opts{o}/$opts{r}_$cf_name.blastp";
				open O,"> $opts{o}/$opts{r}_$cf_name\_best.blastp";
				my @ids=();my %h1=();my %h2=();
				while (my $a=<I>) {
					chomp $a;
					my @itms=split/\t/,$a;
					push @ids,$itms[0];
					if (exists $h1{$itms[0]}) {
						if ($itms[11]>$h2{$itms[0]}) {
							$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];
						}
					}
					else {$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];}
				}
				my @uniq1;
				@uniq1=&uniq1(@ids);
				foreach my $i (@uniq1) {
					print O "$h1{$i}\n";
				}
				close I;close O;
				open I,"< $opts{o}/$cf_name\_$opts{r}.blastp";
				open O,"> $opts{o}/$cf_name\_$opts{r}_best.blastp";
				@ids=();%h1=();%h2=();
				while (my $a=<I>) {
					chomp $a;
					my @itms=split/\t/,$a;
					push @ids,$itms[0];
					if (exists $h1{$itms[0]}) {
						if ($itms[11]>$h2{$itms[0]}) {
							$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];
						}
					}
					else {$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];}
				}

				@uniq1=();
				@uniq1=&uniq1(@ids);

				foreach my $i (@uniq1) {
					print O "$h1{$i}\n";
				}
				close I;close O;
				sub uniq1 {
					my (@arr1)=@_;
					my %count;
					my @uniq11=grep {++$count{$_}<2} @arr1;
					return @uniq11;
				}
				my %h_pair1=();my %h_pair2=();my %h_pair=();
				open I1,"< $opts{o}/$opts{r}_$cf_name\_best.blastp";
				open I2,"< $opts{o}/$cf_name\_$opts{r}_best.blastp";
				while (my $a=<I1>) {
					chomp $a;
					my @it=split/\t/,$a;
					$h_pair1{$it[0]}=$it[1];
				}
				close I1;
				while (my $a=<I2>) {
					chomp $a;
					my @it=split/\t/,$a;
					$h_pair2{$it[0]}=$it[1];
				}
				close I2;
				foreach my $i (keys %h_pair1) {
					$i2=$h_pair1{$i};
					if (exists $h_pair2{$i2} and $h_pair2{$i2} eq $i) {
						$h_pair{$i}=$i2;
					}
				}
				my %py_gene_chr=();my %py_gene_pos=();
				my %cf_gene_chr=();my %cf_gene_pos=();
				open I,"< $opts{i}/$opts{r}_simplified_chr.gff";
				while (my $a=<I>) {
					chomp $a;
					my @it=split/\t/,$a;
					$py_gene_chr{$it[1]}=$it[0];
					$py_gene_pos{$it[1]}=$it[2];
				}
				close I;
				open I,"< $opts{i}/$cf_name\_simplified_chr.gff";
				while (my $a=<I>) {
					chomp $a;
					my @it=split/\t/,$a;
					$cf_gene_chr{$it[1]}=$it[0];
					$cf_gene_pos{$it[1]}=$it[2];
				}
				close I;

				$n=0;
				open O,"> $opts{o}/$opts{r}_$cf_name.lab";
				foreach my $i (keys %h_pair) {
					$n=$n+1;
					$h_gene_id{$i}=$n;
					$h_gene_id{$h_pair{$i}}=$n;
					if (exists $py_gene_chr{$i}) {
						print O "$n\t$py_gene_chr{$i}\t$i\t$h_pair{$i}\n";
					}
				}
				close O;

				open I,"< $opts{color}";
				while (my $a=<I>) {
					chomp $a;
					my @ha=split/\t/,$a;
					$chr_color{$ha[0]}=$ha[1];
				}
				close I;

				open O,"> $opts{o}/$opts{r}_$cf_name.msynt";
				foreach my $i (keys %h_pair) {
					if (exists $py_gene_chr{$i} and exists $py_gene_pos{$i} and exists $h_gene_id{$i} and exists $chr_color{$py_gene_chr{$i}}) {
						print O "$py_gene_chr{$i}\t$i\t$py_gene_pos{$i}\t$h_gene_id{$i}\t1\t$chr_color{$py_gene_chr{$i}}\t$py_gene_chr{$i}\n";
					}
				}
				close O;

				open O,"> $opts{o}/$cf_name\_$opts{r}.msynt";
				foreach my $i (keys %h_pair) {
					$cf=$h_pair{$i};
					if (exists $cf_gene_chr{$cf} and exists $py_gene_chr{$i} and exists $cf_gene_pos{$cf} and exists $h_gene_id{$cf} and exists $chr_color{$py_gene_chr{$i}}) {
						print O "$cf_gene_chr{$cf}\t$i\t$cf_gene_pos{$cf}\t$h_gene_id{$cf}\t1\t$chr_color{$py_gene_chr{$i}}\n";
					}
				}
				close O;
				system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}_$cf_name.msynt >$opts{o}/$opts{r}_$cf_name.sorted.msynt";
				system "sort -k 1,1 -k 3n,3 $opts{o}/$cf_name\_$opts{r}.msynt >$opts{o}/$cf_name\_$opts{r}.sorted.msynt";
				system "perl drawCLGContrib2.4_2.pl $opts{o}/$opts{r}_$cf_name.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$cf_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$cf_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} $opts{o}/$cf_name\_$opts{r}.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$cf_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$cf_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}_$cf_name.svg";
			}
		}
	}
	close P;

	my $dir2=$opts{o};
	my %h_msynt=();
	opendir P,$dir2;
	open O,"> $opts{o}/$opts{r}.msynt";
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_(\S+).sorted.msynt$/ ) {
			if ($1 eq $opts{r}) {
				my $py_msynt=$dir2."/$c";
				open I,"< $py_msynt";
				while (my $a=<I>) {
					chomp $a;
					my @aa=split/\t/,$a;
					if (!exists $h_msynt{$aa[1]}) {
						print O "$a\n";
						$h_msynt{$aa[1]}=0;
					}
				}
				close I;
			}
		}
	}
	close P;
	system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}.msynt >$opts{o}/$opts{r}.sorted.msynt";

	my $dir3=$opts{o};
	my %h_lab=();
	opendir P,$dir3;
	open O,"> $opts{o}/$opts{r}.lab";
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_(\S+).lab$/ ) {
			if ($1 eq $opts{r}) {
				my $py_lab=$dir3."/$c";
				open I,"< $py_lab";
				while (my $a=<I>) {
					chomp $a;
					my @aa=split/\t/,$a;
					if (!exists $h_lab{$aa[2]}) {
						print O "$a\n";
						$h_lab{$aa[2]}=0;
					}
				}
				close I;
			}
		}
	}
	close P;
	#system "perl $opts{s}/drawCLGContrib2.4_2.pl $opts{o}/$opts{r}.sorted.msynt:,algcolor=$opts{o}/$opts{r}.sorted.msynt:type=alg,file=$opts{o}/$opts{r}.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}.svg";

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
}


if ($opts{a} eq "diamond") {
	my $dir1=$opts{i};
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $pep_file=$dir1."/$c";
			my $dirdb = "$opts{i}/$1.pep.database.dmnd";
			if (-e $dirdb) {
			}
			else{
				system "diamond makedb --in $pep_file -d $pep_file.database";
			}
			if ($1 eq $opts{r}) {
				$py_pep=$dir1."/$c";
			}
		}
	}
	close P;

	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $cf_name=$1;
			if ($cf_name ne $opts{r}) {
				my $pep_file=$dir1."/$c";
				system "diamond blastp --query $py_pep --db $pep_file.database --outfmt 6 --out $opts{o}/$opts{r}_$cf_name.blastp  --evalue $opts{e} --threads $opts{t}";
				system "diamond blastp --query $pep_file --db $py_pep.database --outfmt 6 --out $opts{o}/$cf_name\_$opts{r}.blastp  --evalue $opts{e} --threads $opts{t}";
				open I,"< $opts{o}/$opts{r}_$cf_name.blastp";
				open O,"> $opts{o}/$opts{r}_$cf_name\_best.blastp";
				my @ids=();my %h1=();my %h2=();
				while (my $a=<I>) {
					chomp $a;
					my @itms=split/\t/,$a;
					push @ids,$itms[0];
					if (exists $h1{$itms[0]}) {
						if ($itms[11]>$h2{$itms[0]}) {
							$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];
						}
					}
					else {$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];}
				}
				my @uniq2;
				@uniq2=&uniq2(@ids);
				foreach my $i (@uniq2) {
					print O "$h1{$i}\n";
				}
				close I;close O;
				open I,"< $opts{o}/$cf_name\_$opts{r}.blastp";
				open O,"> $opts{o}/$cf_name\_$opts{r}_best.blastp";
				@ids=();%h1=();%h2=();
				while (my $a=<I>) {
					chomp $a;
					my @itms=split/\t/,$a;
					push @ids,$itms[0];
					if (exists $h1{$itms[0]}) {
						if ($itms[11]>$h2{$itms[0]}) {
							$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];
						}
					}
					else {$h1{$itms[0]}=$a;$h2{$itms[0]}=$itms[11];}
				}

				@uniq2=();
				@uniq2=&uniq2(@ids);

				foreach my $i (@uniq2) {
					print O "$h1{$i}\n";
				}
				close I;close O;
				sub uniq2 {
					my (@arr1)=@_;
					my %count;
					my @uniq22=grep {++$count{$_}<2} @arr1;
					return @uniq22;
				}
				my %h_pair1=();my %h_pair2=();my %h_pair=();
				open I1,"< $opts{o}/$opts{r}_$cf_name\_best.blastp";
				open I2,"< $opts{o}/$cf_name\_$opts{r}_best.blastp";
				while (my $a=<I1>) {
					chomp $a;
					my @it=split/\t/,$a;
					$h_pair1{$it[0]}=$it[1];
				}
				close I1;
				while (my $a=<I2>) {
					chomp $a;
					my @it=split/\t/,$a;
					$h_pair2{$it[0]}=$it[1];
				}
				close I2;
				foreach my $i (keys %h_pair1) {
					$i2=$h_pair1{$i};
					if (exists $h_pair2{$i2} and $h_pair2{$i2} eq $i) {
						$h_pair{$i}=$i2;
					}
				}
				my %py_gene_chr=();my %py_gene_pos=();
				my %cf_gene_chr=();my %cf_gene_pos=();
				open I,"< $opts{i}/$opts{r}_simplified_chr.gff";
				while (my $a=<I>) {
					chomp $a;
					my @it=split/\t/,$a;
					$py_gene_chr{$it[1]}=$it[0];
					$py_gene_pos{$it[1]}=$it[2];
				}
				close I;
				open I,"< $opts{i}/$cf_name\_simplified_chr.gff";
				while (my $a=<I>) {
					chomp $a;
					my @it=split/\t/,$a;
					$cf_gene_chr{$it[1]}=$it[0];
					$cf_gene_pos{$it[1]}=$it[2];
				}
				close I;

				$n=0;
				open O,"> $opts{o}/$opts{r}_$cf_name.lab";
				foreach my $i (keys %h_pair) {
					$n=$n+1;
					$h_gene_id{$i}=$n;
					$h_gene_id{$h_pair{$i}}=$n;
					if (exists $py_gene_chr{$i}) {
						print O "$n\t$py_gene_chr{$i}\t$i\t$h_pair{$i}\n";
					}
				}
				close O;

				open I,"< $opts{color}";
				while (my $a=<I>) {
					chomp $a;
					my @ha=split/\t/,$a;
					$chr_color{$ha[0]}=$ha[1];
				}
				close I;

				open O,"> $opts{o}/$opts{r}_$cf_name.msynt";
				foreach my $i (keys %h_pair) {
					if (exists $py_gene_chr{$i} and exists $py_gene_pos{$i} and exists $h_gene_id{$i} and exists $chr_color{$py_gene_chr{$i}}) {
						print O "$py_gene_chr{$i}\t$i\t$py_gene_pos{$i}\t$h_gene_id{$i}\t1\t$chr_color{$py_gene_chr{$i}}\t$py_gene_chr{$i}\n";
					}
				}
				close O;

				open O,"> $opts{o}/$cf_name\_$opts{r}.msynt";
				foreach my $i (keys %h_pair) {
					$cf=$h_pair{$i};
					if (exists $cf_gene_chr{$cf} and exists $py_gene_chr{$i} and exists $cf_gene_pos{$cf} and exists $h_gene_id{$cf} and exists $chr_color{$py_gene_chr{$i}}) {
						print O "$cf_gene_chr{$cf}\t$i\t$cf_gene_pos{$cf}\t$h_gene_id{$cf}\t1\t$chr_color{$py_gene_chr{$i}}\n";
					}
				}
				close O;
				system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}_$cf_name.msynt >$opts{o}/$opts{r}_$cf_name.sorted.msynt";
				system "sort -k 1,1 -k 3n,3 $opts{o}/$cf_name\_$opts{r}.msynt >$opts{o}/$cf_name\_$opts{r}.sorted.msynt";
				system "perl drawCLGContrib2.4_2.pl $opts{o}/$opts{r}_$cf_name.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$cf_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$cf_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} $opts{o}/$cf_name\_$opts{r}.sorted.msynt:,algcolor=$opts{o}/$opts{r}_$cf_name.sorted.msynt:type=alg,file=$opts{o}/$opts{r}_$cf_name.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}_$cf_name.svg";
			}
		}
	}
	close P;

	my $dir2=$opts{o};
	my %h_msynt=();
	opendir P,$dir2;
	open O,"> $opts{o}/$opts{r}.msynt";
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_(\S+).sorted.msynt$/ ) {
			if ($1 eq $opts{r}) {
				my $py_msynt=$dir2."/$c";
				open I,"< $py_msynt";
				while (my $a=<I>) {
					chomp $a;
					my @aa=split/\t/,$a;
					if (!exists $h_msynt{$aa[1]}) {
						print O "$a\n";
						$h_msynt{$aa[1]}=0;
					}
				}
				close I;
			}
		}
	}
	close P;
	system "sort -k 1,1 -k 3n,3 $opts{o}/$opts{r}.msynt >$opts{o}/$opts{r}.sorted.msynt";

	my $dir3=$opts{o};
	my %h_lab=();
	opendir P,$dir3;
	open O,"> $opts{o}/$opts{r}.lab";
	while (my $c=readdir P) {
		if ( $c=~/^(\S+)_(\S+).lab$/ ) {
			if ($1 eq $opts{r}) {
				my $py_lab=$dir3."/$c";
				open I,"< $py_lab";
				while (my $a=<I>) {
					chomp $a;
					my @aa=split/\t/,$a;
					if (!exists $h_lab{$aa[2]}) {
						print O "$a\n";
						$h_lab{$aa[2]}=0;
					}
				}
				close I;
			}
		}
	}
	close P;
	#system "perl $opts{s}/drawCLGContrib2.4_2.pl $opts{o}/$opts{r}.sorted.msynt:,algcolor=$opts{o}/$opts{r}.sorted.msynt:type=alg,file=$opts{o}/$opts{r}.lab,width=40,window=$opts{w},slide=$opts{sl},maxbreak=$opts{m} > $opts{o}/$opts{r}.svg";

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
}

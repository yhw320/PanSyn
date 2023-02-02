#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","o=s","a=s","m=s","t=s","e=s","max=s","h|help");
if(!(defined $opts{i} and defined $opts{o} and defined $opts{a})){

	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	-a	Reference species by name abbreviations (example: HSap)
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-max	The number of best non-self BLASTP hits (default: 5) (default:5)
	-m	Full path where the MCScanX software resides (/software/MCScanX-master/MCScanX)
		*********************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
	die"**********************************************\n
	-i	Full path to the [inputDir] folder containing input files
	-o	Full path to the [outputDir] folder containing output files
	-a	Reference species by name abbreviations (example: HSap)
	Optional:
	-t	diamond alignment threads (default:12)
	-e	diamond alignment evalue (default:1e-3)
	-max	The number of best non-self BLASTP hits (default: 5) (default:5)
	-m	Full path where the MCScanX software resides (/software/MCScanX-master/MCScanX)
		*********************************************\n";
}

if (!(defined $opts{e})) {
	$opts{e}=1e-3;
}
if (!(defined $opts{t})) {
	$opts{t}=12;
}
if (!(defined $opts{max})) {
	$opts{max}=5;
}
if (!(defined $opts{m})) {
	$opts{m}="MCScanX";
}

$dir1=$opts{i};
$num_spe=0;
opendir I,$dir1;
while ($aa=readdir I) {
	if ( $aa=~/^(\S+).pep$/ ) {
		$spe_name=$1;
		$num_spe=$num_spe+1;
		if ($spe_name ne $opts{a}) {
			$h_spe{$spe_name}=0;
			$file=$dir1."/$aa";
			system "diamond makedb --in $file -d $file.database";
			system "mkdir $opts{o}/$opts{a}_$spe_name";
	        system "diamond blastp --query $opts{i}/$opts{a}.pep --db $file.database --outfmt 6 --out $opts{o}/$opts{a}_$spe_name/$opts{a}_$spe_name.blast  --evalue $opts{e} --threads $opts{t} --max-target-seqs $opts{max}";
			system "cat $opts{i}/$opts{a}.gff $opts{i}/$spe_name.gff > $opts{o}/$opts{a}_$spe_name/$opts{a}_$spe_name.gff";
		}
	}
}
close I;
open IN,"<","$opts{i}/$opts{a}.pep" or die;
while (my $a=<IN>) {
	if ($a=~/>(\S+)/) {
		$h_gene_id{$1}=0;
	}
}
close IN;

system "mkdir $opts{o}/MCScanX_result";
foreach my $m (keys %h_spe) {
	chdir "$opts{o}/$opts{a}_$m";
	system "$opts{m} $opts{a}_$m";
	system "cp $opts{o}/$opts{a}_$m/$opts{a}_$m.collinearity $opts{o}/MCScanX_result";
}

for ($num_block=2;$num_block<100000;$num_block=$num_block+1) {
	##跑古老微观共线性第一步
	$dir1="$opts{o}/MCScanX_result";
	opendir I,$dir1;
	system "mkdir $opts{o}/$num_block\_microsyn_genes.blocks";
	while ($aa=readdir I){
		if ($aa=~/^(\S+).collinearity$/) {
			%h=();%h_bn=();
			$file=$dir1."/$aa";
			$name=$1."_all.blocks";
			open O,">$opts{o}/$num_block\_microsyn_genes.blocks/$name";
			open IN,"<","$file" or die;
			while (my $a=<IN>) {
				chomp $a;
				if ($a!~/^#/) {
					$a=~s/ //g;
					$a=~s/://g;
					@items=split/\t/,$a;
					@num_id=split/-/,$items[0];
					$key1=$num_id[0];
					$key2=$num_id[1];
					if (exists $h_gene_id{$items[1]}) {
						$h{$key1}{$key2}=$items[2]."\t".$items[1];	
					}
					if (exists $h_gene_id{$items[2]}) {
						$h{$key1}{$key2}=$items[1]."\t".$items[2];	
					}
				}
			}
			foreach my $i1 ( sort {$a<=>$b} keys %h) {
				$i_before=0;
				foreach my $i2 ( sort {$a<=>$b} keys %{$h{$i1}}) {
					if ($i2>=$i_before) {
						$i_before=$i2;
					}
					else{$i_before=$i_before;}
				}
				$h_bn{$i1}=$i_before;
			}
			foreach my $i1 ( sort {$a<=>$b} keys %h) {
				$max=$h_bn{$i1}-$num_block+1;
				for ($v_ing=0;$v_ing<=$max;$v_ing=$v_ing+1){
					$jia=$v_ing+$num_block-1;
					print O "##\n";
					for ($v=$v_ing;$v<=$jia;$v=$v+1){
						$output=$h{$i1}{$v};
						print O "$output\n";
					}
				}
			}
		}
	close I;
	close O;
	}

	##跑古老微观共线性第二步
	%h_spe2=();%h_an=();my %h1=();undef $combine2;undef $com;undef $com2;undef $combine;
	$dir1="$opts{o}/$num_block\_microsyn_genes.blocks";
	opendir I,$dir1;
	while (my $aa=readdir I) {
		if ( $aa=~/^(\S+)_all.blocks$/ ) {
			$nnae=$1;
			$file=$dir1."/$aa";
			open IN,"<","$file" or die;
			while (my $a=<IN>){
				chomp $a;
				if ($a=~/##/) {
					$combine="non";		
				}
				else{
					my @items=split/\t/,$a;
					if ($combine eq "non") {
						$combine=$items[1];
						$combine2=$items[1];
						$com=$items[0];
						$com2=$items[0];
					}
					else{
						$combine=$combine."\t".$items[1];
						$combine2=$items[1]."\t".$combine2;
						$com=$com."\t".$items[0];
						$com2=$items[0]."\t".$com2;
						my @items2=split/\t/,$combine;
						my @items3=split/\t/,$combine2;
						$items2=@items2;$items3=@items3;
						if ($items2==$num_block and $items3==$num_block) {
							$h_spe2{$nnae}{$combine}=$com;
							$h_spe2{$nnae}{$combine2}=$com2;
							if (exists $h1{$combine}) {
								$h1{$combine}=$h1{$combine}+1;
							}
							elsif(exists $h1{$combine2}){
								$h1{$combine2}=$h1{$combine2}+1;
							}
							else{
								$h1{$combine}=1;
							}
						}
					}
				}
			}
		}
	}
	system "mkdir $opts{o}/$num_block\_microsyn_genes.blocks/conserved_clusters";
	open (O,">$opts{o}/$num_block\_microsyn_genes.blocks/conserved_clusters/final_conserved.clusters");
	$nn=$num_spe-1;
	foreach $i (keys %h1) {
		if ($h1{$i}>= $nn) {
			#my @zhuan=();my @zhuan2=();
			print O "$i\n";
			$h_an{$i}=0;
			#@zhuan=split/\t/,$i;
			#@zhuan2 = reverse (@zhuan);
			#$nne=join("\t",@zhuan2);
			#$h_an{$nne}=0;
		}
	}
	close I;
	close O;
	foreach my $na (keys %h_spe2) {
		open (O,">$opts{o}/$num_block\_microsyn_genes.blocks/conserved_clusters/$na\_conserved.clusters");
			foreach my $id (keys %{$h_spe2{$na}}) {
				if (exists $h_an{$id}) {
					print O "###\n";
					my @it1=split/\t/,$id;
					my @it2=split/\t/,$h_spe2{$na}{$id};
					$haha=$num_block-1;
					for ($v=0;$v<=$haha;$v=$v+1) {
						print O "$it1[$v]\t$it2[$v]\n";
					}
				}
			}
		print O "###\n";
		close O;
	}
	$hash_size=keys%h_an;
	if ($hash_size==0) {
		system "rm -r $opts{o}/$num_block\_microsyn_genes.blocks";
		exit;
	}
}

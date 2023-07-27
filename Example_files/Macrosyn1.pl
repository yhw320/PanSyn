#!/usr/bin/perl
use Getopt::Long;
my %opts;
use List::Util qw(max);
GetOptions(\%opts,"i1=s","t=s","a=s","e=s","o1=s","h|help");
if (!( defined $opts{i1} and defined $opts{o1} and defined $opts{a})) {
		die "************************************************\n
	-i1	Full path to the [inputDir1] directory containing input files
	-o1	Full path to the [outputDir1] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-e	Protein alignment Evalue (default:0.001/1e-5)
	-t	Protein alignment threads (default:12)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i1	Full path to the [inputDir1] directory containing input files
	-o1	Full path to the [outputDir1] directory containing output files
	-a	Specify the protein alignment software (It can be set to 'diamond' or 'blast')
	Optional:
	-e	Protein alignment Evalue (default:0.001/1e-5)
	-t	Protein alignment threads (default:12)
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
if ($opts{o1} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o1} =~ s/$slash$//;
}
####


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



system "mkdir $opts{o1}/ids_match";
print "The [$opts{o1}/ids_match] directory has been created!\n";

system "mkdir $opts{o1}/score";
print "The [$opts{o1}/score] directory has been created!\n";

system "mkdir $opts{o1}/single_cluster";
print "The [$opts{o1}/single_cluster] directory has been created!\n";

system "mkdir $opts{o1}/blastp";
print "The [$opts{o1}/blastp] directory has been created!\n";

system "mkdir $opts{o1}/formatted_peps";
print "The [$opts{o1}/formatted_peps directory has been created!\n";


if ($opts{a} eq "blast") {
	my $dir1=$opts{i1}."/peps";
	my $spe_num=0;
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $file=$dir1."/$c";
			testDB_script1_get_formatted_input_pep_file1($1,$file); 
			$spe_num=$spe_num+1;
		}
	}
	close P;

	##sub testDB_script1_get_formatted_input_pep_file1.pl
	sub testDB_script1_get_formatted_input_pep_file1 {
		my $spe=$_[0];
		my $filename=$_[1];
		my $n=0;
		open IN,"<$filename" or die("Could not open $filename\n"); 
		open OUT1,"> $opts{o1}/formatted_peps/$spe.pep.fa";
		open OUT2,"> $opts{o1}/ids_match/$spe.ids.match";
		while (my $a=<IN>) {
			chomp $a;
			if ($a=~/>(\S+)/) {
				$n++;
				my $rid=$spe."_$spe"."_amo$n";
				print OUT1 ">$rid\n";
				print OUT2 "$1\t$rid\n";
			}
			else {print OUT1 "$a\n";}
		}
		close IN;
		close OUT1;
		close OUT2;
	}
	###################################################
	system "cat $opts{o1}/formatted_peps/*fa > $opts{o1}/formatted_peps/all.fa";
	$fileExist = -e "$opts{o1}/formatted_peps/all.fa.psq";
	if ( $fileExist ) {
	}
	else {
		system "makeblastdb -in $opts{o1}/formatted_peps/all.fa -dbtype prot";
	}


	my $dir2=$opts{o1}."/formatted_peps";
	opendir F,$dir2;
	while (my $d=readdir F) {
		if ( $d=~/^(\S+).pep.fa$/ ) {
			my $file2=$dir2."/$d";
			system "blastp -query $opts{o1}/formatted_peps/$1.pep.fa -db $opts{o1}/formatted_peps/all.fa -outfmt 6 -out $opts{o1}/blastp/$1.blastp -evalue $opts{e} -num_threads $opts{t}";
			system "cut -f1,2,12 $opts{o1}/blastp/$1.blastp | sort -k 1 -nk 3 > $opts{o1}/score/$1.score";
			#system"less $opts{o1}/blastp/$1.blastp |awk '{print \$1\"\t\"\$2\"\t\"\$12}'|sort -k 1 -nk 3 > $opts{o1}/score/$1.score";
			cluster_out1($1); 
		}
	}
	close F;
	##sub cluster_out1.pl
	sub cluster_out1{
		my $n=0;my %h_sub;
		my $spename=$_[0];
		open IN,"<$opts{o1}/score/$spename.score";
		while (my $a=<IN>){
			chomp $a;
			my @items=split/\t/,$a;
			$h_sub{$items[0]}=1;
			$h_sub{$items[1]}=1;
		}
		open OUT,">$opts{o1}/single_cluster/$spename.cluster2";
		foreach my $key (sort keys %h_sub) {
			if ($key=~/^$spename\_/) {
				$n=$n+1;
				print OUT "$n\t$key\n";
			}
		 }
		close IN;
		close OUT;
	}
	######################3
	my @b;my $pair;my @single_spe;my $spe1;my $spe2;my $addname;
	system"mkdir $opts{o1}/cluster1";
	print "The [$opts{o1}/cluster1] directory has been created!\n";

	system"mkdir $opts{o1}/table";
	print "The [$opts{o1}/table] directory has been created!\n";

	open T5,"<$opts{i1}/tree.nwk" or die("Could not open $opts{i1}/tree.nwk\n"); 
	while (my $a=<T5>){
		chomp $a;
		$a=~s/;//g;
		while ($a=~/,/) {
			@b=$a=~/(\w+,\w+)/g;
			foreach $pair (@b) {
				@single_spe=split/,/,$pair;
				$spe1=$single_spe[0];
				$spe2=$single_spe[1];
				$addname=$spe1."_".$spe2;
				system"cat $opts{o1}/score/$spe1.score $opts{o1}/score/$spe2.score > $opts{o1}/score/$addname.score";
				MacroSystenic_firstCluster1($spe1,$spe2,$addname);
				count_inter_outer_scores1($addname);
				MacroSystenic_step21($addname);
			}
			@b=();
			$a=~s/\((\w+),(\w+)\)/$1\_$2/g;
		}
		my $path_cate=$opts{i1}."/species_classification.table";
		get_ancietS1($addname,$path_cate);
	}
	close T5;
	##sub get_ancietS1_2.pl
	sub get_ancietS1{
		$nn_num=1;$nnn1=1;$nnn2=2;$nnn3=3;
		my $cl=$_[0];
		my $pathdy=$_[1];
		my @it;my $aa;my $bb;my %h1;
		open T,"<","$pathdy" or die;
		while(my $a=<T>){
			if ($a!~/^Species/) {
				chomp $a;
				@it=split/\t/,$a;
				$aa=$it[0];
				$bb=$it[1];
				$h1{$aa}=$bb;
				$h_name{$bb}=$aa;
			}
		}
		foreach $ia (keys %h_name) {
			$hnn{$nn_num}=$ia;
			$nn_num=$nn_num+1;
		}
		open I,"<$opts{o1}/single_cluster/$cl.cluster2" or die;
		open OUT,">$opts{o1}/ancient_gene.families";
		while(<I>){
			chomp;
			my %count=();
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f1(@b){
				my @nms=split/_/,$f1;
				my $p=$nms[0];
				my $cla=$h1{$p};
				$count{$cla}++;
			}
			my $out=$count{$hnn{$nnn1}};my $yuankou=$count{$hnn{$nnn2}};my $houkou=$count{$hnn{$nnn3}};
			if(($out>=2 && $yuankou>=2) || ($out>=2 && $houkou>=2) || ($yuankou>=2 && $houkou>=2)){
				print OUT "$_\n";
			}	
		}
		close I;close T;
		close OUT;
	}
	##sub count_inter_outer_scores1.pl
	sub count_inter_outer_scores1{
		my $aname=$_[0];
		my %list=();
		my %outMaxScore=();
		my %interMaxScore=();
		open IN,"<","$opts{o1}/cluster1/$aname.cluster1" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f(@b){
				$list{$f}=$a[0];
			}
		}
		close IN;
		open IN,"<","$opts{o1}/score/$aname.score" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			next if ($a[1] eq $a[0]);
			if(exists $list{$a[0]} && exists $list{$a[1]}){
				my $c1=$list{$a[0]};
				my $c2=$list{$a[1]};
				if($c1!=$c2){
					if(!exists $interMaxScore{$c1}{$c2}){
						$interMaxScore{$c1}{$c2}=$a[2];
					}else{
						$interMaxScore{$c1}{$c2}=$a[2] if ($a[2]>$interMaxScore{$c1}{$c2});
					}
				}
			}
			if(exists $list{$a[0]} && !exists $list{$a[1]}){
				my $c3=$list{$a[0]};
				if(!exists $outMaxScore{$c3}){
					$outMaxScore{$c3}=$a[2];
				}else{
					$outMaxScore{$c3}=$a[2] if ($a[2]>$outMaxScore{$c3});
				}
			}
			if(!exists $list{$a[0]} && exists $list{$a[1]}){
				my $c4=$list{$a[1]};
				if(!exists $outMaxScore{$c4}){
					$outMaxScore{$c4}=$a[2];
				}else{
					$outMaxScore{$c4}=$a[2] if ($a[2]>$outMaxScore{$c4});
				}
			}	 
		} 
		close IN;
		my %final;
		foreach my $k1(keys%interMaxScore){
			foreach my $k2(keys%{$interMaxScore{$k1}}){
				push @{$final{'hh'}},[$interMaxScore{$k1}{$k2},$k1,$k2];
			}   
		}
		my %pres=();
		open OUT,">$opts{o1}/table/$aname.cluster1.table";
		foreach my $key1(keys%final){
			@{$final{$key1}}=sort{$b->[0] <=> $a->[0]}@{$final{$key1}};
			for(my $i=0;$i<@{$final{$key1}};$i++){
				my ($oc1,$oc2)=(0,0);
				my $scoress=${$final{$key1}}[$i][0];
				my $id1=${$final{$key1}}[$i][1];
				my $id2=${$final{$key1}}[$i][2];
				next if(exists $pres{$id1}{$id2} || exists $pres{$id2}{$id1});
				if(exists $outMaxScore{$id1}){
					$oc1=$outMaxScore{$id1};
				}
				if(exists $outMaxScore{$id2}){
					$oc2=$outMaxScore{$id2};
				}
				my $maxoc=max($oc1,$oc2);
				if($scoress>=$maxoc){
					print OUT "$id1\t$id2\t$scoress\t$maxoc\n";
					$pres{$id1}{$id2}=1;
				}
			}
		}
		close OUT;
	}
	##sub MacroSystenic_firstCluster1.pl
	sub MacroSystenic_firstCluster1{
		my $spea=$_[0];
		my $speb=$_[1];
		my $adname=$_[2];
		my %genelist1=();
		my %genelist2=();
		my %hCluster1=();
		my %hCluster2=();
	#====================================
		open IN,"<","$opts{o1}/single_cluster/$spea.cluster2" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $n1 (0..$#b) {
				$genelist1{$b[$n1]}=$a[0];
				$hCluster1{$a[0]}{$b[$n1]}=1;
			}
		}
		close IN;
	#====================================
		open IN,"<","$opts{o1}/single_cluster/$speb.cluster2" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $n2 (0..$#b) {
				$genelist2{$b[$n2]}=$a[0];
				$hCluster2{$a[0]}{$b[$n2]}=1;
			}
		}
		close IN;
	#=======================================
		my %ScoreStore1=();
		my %ScoreStore2=();
		open IN,"<","$opts{o1}/score/$adname.score" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			if(exists $genelist1{$a[0]} && exists $genelist2{$a[1]}){
				push @{$ScoreStore1{$genelist1{$a[0]}}},[$a[2],$genelist2{$a[1]}];
			}
			if(exists $genelist2{$a[0]} && exists $genelist1{$a[1]}){
				push @{$ScoreStore2{$genelist2{$a[0]}}},[$a[2],$genelist1{$a[1]}];
			}
		}
		close IN;
	#==========================================
		my %MaxScore1=();
		my %MaxScore2=();
		foreach my $key1 (keys%ScoreStore1) {
			@{$ScoreStore1{$key1}}=sort{$b->[0]<=>$a->[0]}@{$ScoreStore1{$key1}};
			$MaxScore1{$key1}=${$ScoreStore1{$key1}}[0][1];
		}
		foreach my $key2 (keys%ScoreStore2) {
			 @{$ScoreStore2{$key2}}=sort{$b->[0]<=>$a->[0]}@{$ScoreStore2{$key2}};
			$MaxScore2{$key2}=${$ScoreStore2{$key2}}[0][1];
		}
		%ScoreStore1=();
		%ScoreStore2=();
	#=================================================
		my $count=0;
		my %final=();
		my %repeat1=();
		my %repeat2=();
		foreach my $key3 (sort{$a<=>$b}keys%MaxScore1) {
			my $matchkey=$MaxScore1{$key3};
			if(exists $MaxScore2{$matchkey}){
				if($MaxScore2{$matchkey}==$key3){
					$repeat1{$key3}=1;
					$repeat2{$matchkey}=1;
					$count++;
					foreach my $k1 (sort keys%{$hCluster1{$key3}}) {
						$final{$count}{$k1}=1;
					}
					foreach my $k2 (sort keys%{$hCluster2{$matchkey}}) {
						$final{$count}{$k2}=1;
				   }
				}
			}
		}
	#========================================================
		open OUT,">$opts{o1}/cluster1/$adname.cluster1";
		foreach my $f1 (sort{$a<=>$b}keys%final) {
			print OUT "$f1\t";
			foreach my $f2 (sort keys%{$final{$f1}}) {
				print OUT "$f2,";
			}
			print OUT "\n";
		}
		foreach my $f3 (sort{$a<=>$b}keys%hCluster1) {
			if(!exists $repeat1{$f3}){
				$count++;
				print OUT "$count\t";
				foreach my $f4 (sort keys%{$hCluster1{$f3}}) {
					print OUT "$f4,";
				}
				print OUT "\n";
			}
		}  
		foreach my $f5 (sort{$a<=>$b}keys%hCluster2) {
			if(!exists $repeat2{$f5}){
				$count++;
				print OUT "$count\t";
				foreach my $f6 (sort keys%{$hCluster2{$f5}}) {
					print OUT "$f6,";
				}
				print OUT "\n";
			}
		}
		close OUT;
	}
	##sub MacroSystenic_step2.pl
	sub MacroSystenic_step21{
		my $aname2=$_[0];
		my %genes=();
		open IN,"<","$opts{o1}/cluster1/$aname2.cluster1" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f(@b){
				$genes{$a[0]}{$f}=1;
			}
		}
		my %matchs=();
		open IN,"<","$opts{o1}/table/$aname2.cluster1.table" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			$matchs{$a[0]}{$a[1]}=1;
			$matchs{$a[1]}{$a[0]}=1;
		}
		close IN;
		my %hCluster=();
		my %list=();
		my $count=0;
		open IN,"<","$opts{o1}/table/$aname2.cluster1.table" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			if(!exists $list{$a[0]} && !exists $list{$a[1]}){
				$count++;
				$list{$a[0]}=$count;
				$list{$a[1]}=$count;
				$hCluster{$count}{$a[0]}=1;
				$hCluster{$count}{$a[1]}=1;
			}
			if(exists $list{$a[0]} && !exists $list{$a[1]}){
				my $cluid1=$list{$a[0]};
				my $flag1=0;
				foreach my $f1(keys%{$hCluster{$cluid1}}){
					$flag1++ unless(exists $matchs{$f1}{$a[1]});
				}
				if($flag1==0){
					$list{$a[1]}=$cluid1;
					$hCluster{$cluid1}{$a[1]}=1;
				}else{
					$count++;
					$list{$a[1]}=$count;
					$hCluster{$count}{$a[1]}=1;
				}
			}
			if(!exists $list{$a[0]} && exists $list{$a[1]}){
				my $cluid2=$list{$a[1]};
				my $flag2=0;
				foreach my $f2(keys%{$hCluster{$cluid2}}){
					$flag2++ unless(exists $matchs{$f2}{$a[0]});
				}
				if($flag2==0){
					$list{$a[0]}=$cluid2;
					$hCluster{$cluid2}{$a[0]}=1;
				}else{
					$count++;
					$list{$a[0]}=$count;
					$hCluster{$count}{$a[0]}=1;
				}
			}	
		}
		close IN;
		open OUT,">$opts{o1}/single_cluster/$aname2.cluster2";
		my %exi=();
		foreach my $key1(sort{$a<=>$b}keys%hCluster){
			print OUT "$key1\t";
			foreach my $key2(keys%{$hCluster{$key1}}){
				$exi{$key2}=1;
				foreach my $r1(keys%{$genes{$key2}}){
					print OUT "$r1,";
				}
			}
			print OUT "\n";
		}
		foreach my $key3(keys%genes){
			next if (exists $exi{$key3});
			$count++;
			print OUT "$count\t";
			foreach my $key4(keys%{$genes{$key3}}){
				print OUT "$key4,";
			}
			print OUT"\n";
		}
		close OUT;
	}
	caculate_family1($addname);
	#perl caculate_family.pl -i ancient.cluster2 -d DMel_PYes_MMus_BFlo_NVec_AQue.cluster2 -o try
	sub caculate_family1{
		my $cll=$_[0];
		open I,"<$opts{o1}/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene.families\n"); 
		my $line = "";
		my %ancient_spes_num = ();
		my $total_ancient = 0;
		while($line = <I>){
			my %species_status = ();
			$total_ancient++;
			chomp $line;
			my @lines = split /\t/,$line;
			my $entries = $lines[1];
			my @entries = split /\,/,$entries;
			foreach my $entry (@entries){
				my @entry_infos = split /\_/,$entry;
				my $species = $entry_infos[0];
				$species_status{$species} = 1;
			}
			foreach my $species (keys %species_status){
				$ancient_spes_num{$species}++;
			}
		}
		close I;

		open D,"<$opts{o1}/single_cluster/$cll.cluster2" or die;
		my $line1 = "";
		my %ancient_spes_num1 = ();
		my $total_ancient1 = 0;
		while($line1 = <D>){
			my %species_status1 = ();
			$total_ancient1++;
			chomp $line1;
			my @lines1 = split /\t/,$line1;
			my $entries1 = $lines1[1];
			my @entries1 = split /\,/,$entries1;
			foreach my $entry1 (@entries1){
				my @entry_infos1 = split /\_/,$entry1;
				my $species1 = $entry_infos1[0];
				$species_status1{$species1} = 1;
			}
			foreach my $species1 (keys %species_status1){
				$ancient_spes_num1{$species1}++;
			}
		}
		close D;
		open O,">$opts{o1}/gene_families.situation";
		print O "$total_ancient ancient families in total.\n";
		print O "Species\tA:ancient_gene_families_of_each_species\tB:total_ancient_gene_families\tC:total_gene_families_of_each_species\tA/Bpercentage\tA/Cpercentage\n";
		foreach my $species (sort keys %ancient_spes_num){
			my $percentage_ancient = $ancient_spes_num{$species} / $total_ancient;
			my $percentage_ancient2 = $ancient_spes_num{$species} / $ancient_spes_num1{$species};
			print O "$species\t$ancient_spes_num{$species}\t$total_ancient\t$ancient_spes_num1{$species}\t";
			printf O "%0.3f\t", $percentage_ancient;
			printf O "%0.3f\n", $percentage_ancient2;
		}
		close O;
		open C,">$opts{o1}/family_situation.table";
		print C "Species\tTotal_ancient_gene_families\tTotal_gene_families\tType1\tType2\n";
		foreach my $species (sort keys %ancient_spes_num){
			my $percentage_ancient = $ancient_spes_num{$species} / $total_ancient;
			my $percentage_ancient2 = $ancient_spes_num{$species} / $ancient_spes_num1{$species};
			print C "$species\t";
			printf C "%0.3f\t", $percentage_ancient;
			printf C "%0.3f\t", $percentage_ancient2;
			print C "total_ancient_gene_families\t";
			print C "total_gene_families\n";
		}
		close C;

	}
	#choose_gene_in_ancient_gene_families
	#@gr=split/\_/,$addname;
	my %hg;
	open I,"<$opts{o1}/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene.families\n"); 
	while(my $a = <I>){
		my @gr1=();my @single_gene_id=();
		chomp $a;
		@gr1 = split /\t/,$a;
		@single_gene_id = split /,/,$gr1[1];
		foreach my $gi (@single_gene_id){
			my @gr2 =();
			@gr2 = split /\_/,$gi;
			my $sp1=$gr2[0];
			if (exists $hg{$sp1}) {
				$hg{$sp1}=$hg{$sp1}."\t".$gi;
			}
			else{
				$hg{$sp1}=$gi;
			}
		}
	}
	close I;

	system("mkdir $opts{o1}/genes-in-ancient_gene_families");

	foreach my $i (sort keys %hg) {
		my $grr=$hg{$i};
		my @gr3=();my %h_seq=();
		@gr3 = split /\t/,$grr;
		open I2,"<$dir1/$i.pep";
		while(my $b = <I2>){
			chomp $b;
			$b=~s/\r//g;
			if ($b=~/>(\S+)/) { 
				$iddd=$1; 
			} 
			else { 
				$h_seq{$iddd} .= $b; 
			} 
		}
		close I2;
		open I,"<$opts{o1}/ids_match/$i.ids.match";
		my %h_xp=();
		while(my $a = <I>){
			chomp $a;my @xp=();
			@xp=split/\t/,$a;
			$h_xp{$xp[1]}=$xp[0];
		}
		close I;
		open O,">$opts{o1}/genes-in-ancient_gene_families/$i.ancient_gene_sequences";
		foreach my $m (@gr3) {
			if (exists $h_xp{$m}) {
				my $mmm=$h_xp{$m};
				if (exists $h_seq{$mmm}) {
					print O ">$mmm\n$h_seq{$mmm}\n";
				}
			}
		}
		close O;
	}



	system("mkdir $opts{o1}/ancient_gene_families.analysis");
	print "The [$opts{o1}/ancient_gene_families.analysis] directory has been created!\n";

	system("mv $opts{o1}/single_cluster/$addname.cluster2 $opts{o1}/ancient_gene_families.analysis/all_gene.families");
	system("mv $opts{o1}/ancient_gene.families $opts{o1}/ancient_gene_families.analysis/ancient_gene.families");
	system("mv $opts{o1}/gene_families.situation $opts{o1}/ancient_gene_families.analysis/gene_families.situation");
	system("mv $opts{o1}/ancient_gene_families.analysis/gene_families.situation $opts{o1}/ancient_gene_families.analysis/gene_families_summary.table");
	system "mv $opts{o1}/genes-in-ancient_gene_families $opts{o1}/ancient_gene_families.analysis/ancient_gene_families.genes";
	system "mv $opts{o1}/family_situation.table $opts{o1}/ancient_gene_families.analysis/family_situation.table";
}

if ($opts{a} eq "diamond") {
	my $dir1=$opts{i1}."/peps";
	my $spe_num=0;
	opendir P,$dir1;
	while (my $c=readdir P) {
		if ( $c=~/^(\S+).pep$/ ) {
			my $file=$dir1."/$c";
			testDB_script1_get_formatted_input_pep_file($1,$file); 
			$spe_num=$spe_num+1;
		}
	}
	close P;

	##sub testDB_script1_get_formatted_input_pep_file.pl
	sub testDB_script1_get_formatted_input_pep_file {
		my $spe=$_[0];
		my $filename=$_[1];
		my $n=0;
		open IN,"<$filename" or die("Could not open $filename\n"); 
		open OUT1,"> $opts{o1}/formatted_peps/$spe.pep.fa";
		open OUT2,"> $opts{o1}/ids_match/$spe.ids.match";
		while (my $a=<IN>) {
			chomp $a;
			if ($a=~/>(\S+)/) {
				$n++;
				my $rid=$spe."_$spe"."_amo$n";
				print OUT1 ">$rid\n";
				print OUT2 "$1\t$rid\n";
			}
			else {print OUT1 "$a\n";}
		}
		close IN;
		close OUT1;
		close OUT2;
	}
	###################################################
	system "cat $opts{o1}/formatted_peps/*fa > $opts{o1}/formatted_peps/all.fa";
	$fileExist = -e "$opts{o1}/formatted_peps/all.fa.database.dmnd";
	if ( $fileExist ) {
	}
	else {
	system "diamond makedb --in $opts{o1}/formatted_peps/all.fa -d $opts{o1}/formatted_peps/all.fa.database";
	}


	my $dir2=$opts{o1}."/formatted_peps";
	opendir F,$dir2;
	while (my $d=readdir F) {
		if ( $d=~/^(\S+).pep.fa$/ ) {
			my $file2=$dir2."/$d";
			system "diamond blastp --query $opts{o1}/formatted_peps/$1.pep.fa --db $opts{o1}/formatted_peps/all.fa.database --outfmt 6 --out $opts{o1}/blastp/$1.blastp  --evalue $opts{e} --threads $opts{t}";
			system "cut -f1,2,12 $opts{o1}/blastp/$1.blastp | sort -k 1 -nk 3 > $opts{o1}/score/$1.score";
			#system"less $opts{o1}/blastp/$1.blastp |awk '{print \$1\"\t\"\$2\"\t\"\$12}'|sort -k 1 -nk 3 > $opts{o1}/score/$1.score";
			cluster_out($1); 
		}
	}
	close F;
	##sub cluster_out.pl
	sub cluster_out{
		my $n=0;my %h_sub;
		my $spename=$_[0];
		open IN,"<$opts{o1}/score/$spename.score";
		while (my $a=<IN>){
			chomp $a;
			my @items=split/\t/,$a;
			$h_sub{$items[0]}=1;
			$h_sub{$items[1]}=1;
		}
		open OUT,">$opts{o1}/single_cluster/$spename.cluster2";
		foreach my $key (sort keys %h_sub) {
			if ($key=~/^$spename\_/) {
				$n=$n+1;
				print OUT "$n\t$key\n";
			}
		 }
		close IN;
		close OUT;
	}
	######################3
	my @b;my $pair;my @single_spe;my $spe1;my $spe2;my $addname;

	system"mkdir $opts{o1}/cluster1";
	print "The [$opts{o1}/cluster1] directory has been created!\n";
	
	system"mkdir $opts{o1}/table";
	print "The [$opts{o1}/table] directory has been created!\n";

	open T5,"<$opts{i1}/tree.nwk" or die("Could not open $opts{i1}/tree.nwk\n"); 
	while (my $a=<T5>){
		chomp $a;
		$a=~s/;//g;
		while ($a=~/,/) {
			@b=$a=~/(\w+,\w+)/g;
			foreach $pair (@b) {
				@single_spe=split/,/,$pair;
				$spe1=$single_spe[0];
				$spe2=$single_spe[1];
				$addname=$spe1."_".$spe2;
				system"cat $opts{o1}/score/$spe1.score $opts{o1}/score/$spe2.score > $opts{o1}/score/$addname.score";
				MacroSystenic_firstCluster($spe1,$spe2,$addname);
				count_inter_outer_scores($addname);
				MacroSystenic_step2($addname);
			}
			@b=();
			$a=~s/\((\w+),(\w+)\)/$1\_$2/g;
		}
		my $path_cate=$opts{i1}."/species_classification.table";
		get_ancietS($addname,$path_cate);
	}
	close T5;
	##sub get_ancietS_2.pl
	sub get_ancietS{
		$nn_num=1;$nnn1=1;$nnn2=2;$nnn3=3;
		my $cl=$_[0];
		my $pathdy=$_[1];
		my @it;my $aa;my $bb;my %h1;
		open T,"<","$pathdy" or die;
		while(my $a=<T>){
			if ($a!~/^Species/) {
				chomp $a;
				@it=split/\t/,$a;
				$aa=$it[0];
				$bb=$it[1];
				$h1{$aa}=$bb;
				$h_name{$bb}=$aa;
			}
		}
		foreach $ia (keys %h_name) {
			$hnn{$nn_num}=$ia;
			$nn_num=$nn_num+1;
		}
		open I,"<$opts{o1}/single_cluster/$cl.cluster2" or die;
		open OUT,">$opts{o1}/ancient_gene.families";
		while(<I>){
			chomp;
			my %count=();
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f1(@b){
				my @nms=split/_/,$f1;
				my $p=$nms[0];
				my $cla=$h1{$p};
				#if (exists $count{$cla}) {
				#	$count{$cla}=$count{$cla}+1;
				#}
				#else{$count{$cla}=0;}
				$count{$cla}++;
			}
			my $out=$count{$hnn{$nnn1}};my $yuankou=$count{$hnn{$nnn2}};my $houkou=$count{$hnn{$nnn3}};
			my $nq=2;
			if(($out>=$nq && $yuankou>=$nq) || ($out>=$nq && $houkou>=$nq) || ($yuankou>=$nq && $houkou>=$nq)){
				print OUT "$_\n";
			}	
		}
		close I;close T;
		close OUT;
	}
	##sub count_inter_outer_scores.pl
	sub count_inter_outer_scores{
		my $aname=$_[0];
		my %list=();
		my %outMaxScore=();
		my %interMaxScore=();
		open IN,"<","$opts{o1}/cluster1/$aname.cluster1" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f(@b){
				$list{$f}=$a[0];
			}
		}
		close IN;
		open IN,"<","$opts{o1}/score/$aname.score" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			next if ($a[1] eq $a[0]);
			if(exists $list{$a[0]} && exists $list{$a[1]}){
				my $c1=$list{$a[0]};
				my $c2=$list{$a[1]};
				if($c1!=$c2){
					if(!exists $interMaxScore{$c1}{$c2}){
						$interMaxScore{$c1}{$c2}=$a[2];
					}else{
						$interMaxScore{$c1}{$c2}=$a[2] if ($a[2]>$interMaxScore{$c1}{$c2});
					}
				}
			}
			if(exists $list{$a[0]} && !exists $list{$a[1]}){
				my $c3=$list{$a[0]};
				if(!exists $outMaxScore{$c3}){
					$outMaxScore{$c3}=$a[2];
				}else{
					$outMaxScore{$c3}=$a[2] if ($a[2]>$outMaxScore{$c3});
				}
			}
			if(!exists $list{$a[0]} && exists $list{$a[1]}){
				my $c4=$list{$a[1]};
				if(!exists $outMaxScore{$c4}){
					$outMaxScore{$c4}=$a[2];
				}else{
					$outMaxScore{$c4}=$a[2] if ($a[2]>$outMaxScore{$c4});
				}
			}	 
		} 
		close IN;
		my %final;
		foreach my $k1(keys%interMaxScore){
			foreach my $k2(keys%{$interMaxScore{$k1}}){
				push @{$final{'hh'}},[$interMaxScore{$k1}{$k2},$k1,$k2];
			}   
		}
		my %pres=();
		open OUT,">$opts{o1}/table/$aname.cluster1.table";
		foreach my $key1(keys%final){
			@{$final{$key1}}=sort{$b->[0] <=> $a->[0]}@{$final{$key1}};
			for(my $i=0;$i<@{$final{$key1}};$i++){
				my ($oc1,$oc2)=(0,0);
				my $scoress=${$final{$key1}}[$i][0];
				my $id1=${$final{$key1}}[$i][1];
				my $id2=${$final{$key1}}[$i][2];
				next if(exists $pres{$id1}{$id2} || exists $pres{$id2}{$id1});
				if(exists $outMaxScore{$id1}){
					$oc1=$outMaxScore{$id1};
				}
				if(exists $outMaxScore{$id2}){
					$oc2=$outMaxScore{$id2};
				}
				my $maxoc=max($oc1,$oc2);
				if($scoress>=$maxoc){
					print OUT "$id1\t$id2\t$scoress\t$maxoc\n";
					$pres{$id1}{$id2}=1;
				}
			}
		}
		close OUT;
	}
	##sub MacroSystenic_firstCluster.pl
	sub MacroSystenic_firstCluster{
		my $spea=$_[0];
		my $speb=$_[1];
		my $adname=$_[2];
		my %genelist1=();
		my %genelist2=();
		my %hCluster1=();
		my %hCluster2=();
	#====================================
		open IN,"<","$opts{o1}/single_cluster/$spea.cluster2" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $n1 (0..$#b) {
				$genelist1{$b[$n1]}=$a[0];
				$hCluster1{$a[0]}{$b[$n1]}=1;
			}
		}
		close IN;
	#====================================
		open IN,"<","$opts{o1}/single_cluster/$speb.cluster2" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $n2 (0..$#b) {
				$genelist2{$b[$n2]}=$a[0];
				$hCluster2{$a[0]}{$b[$n2]}=1;
			}
		}
		close IN;
	#=======================================
		my %ScoreStore1=();
		my %ScoreStore2=();
		open IN,"<","$opts{o1}/score/$adname.score" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			if(exists $genelist1{$a[0]} && exists $genelist2{$a[1]}){
				push @{$ScoreStore1{$genelist1{$a[0]}}},[$a[2],$genelist2{$a[1]}];
			}
			if(exists $genelist2{$a[0]} && exists $genelist1{$a[1]}){
				push @{$ScoreStore2{$genelist2{$a[0]}}},[$a[2],$genelist1{$a[1]}];
			}
		}
		close IN;
	#==========================================
		my %MaxScore1=();
		my %MaxScore2=();
		foreach my $key1 (keys%ScoreStore1) {
			@{$ScoreStore1{$key1}}=sort{$b->[0]<=>$a->[0]}@{$ScoreStore1{$key1}};
			$MaxScore1{$key1}=${$ScoreStore1{$key1}}[0][1];
		}
		foreach my $key2 (keys%ScoreStore2) {
			 @{$ScoreStore2{$key2}}=sort{$b->[0]<=>$a->[0]}@{$ScoreStore2{$key2}};
			$MaxScore2{$key2}=${$ScoreStore2{$key2}}[0][1];
		}
		%ScoreStore1=();
		%ScoreStore2=();
	#=================================================
		my $count=0;
		my %final=();
		my %repeat1=();
		my %repeat2=();
		foreach my $key3 (sort{$a<=>$b}keys%MaxScore1) {
			my $matchkey=$MaxScore1{$key3};
			if(exists $MaxScore2{$matchkey}){
				if($MaxScore2{$matchkey}==$key3){
					$repeat1{$key3}=1;
					$repeat2{$matchkey}=1;
					$count++;
					foreach my $k1 (sort keys%{$hCluster1{$key3}}) {
						$final{$count}{$k1}=1;
					}
					foreach my $k2 (sort keys%{$hCluster2{$matchkey}}) {
						$final{$count}{$k2}=1;
				   }
				}
			}
		}
	#========================================================
		open OUT,">$opts{o1}/cluster1/$adname.cluster1";
		foreach my $f1 (sort{$a<=>$b}keys%final) {
			print OUT "$f1\t";
			foreach my $f2 (sort keys%{$final{$f1}}) {
				print OUT "$f2,";
			}
			print OUT "\n";
		}
		foreach my $f3 (sort{$a<=>$b}keys%hCluster1) {
			if(!exists $repeat1{$f3}){
				$count++;
				print OUT "$count\t";
				foreach my $f4 (sort keys%{$hCluster1{$f3}}) {
					print OUT "$f4,";
				}
				print OUT "\n";
			}
		}  
		foreach my $f5 (sort{$a<=>$b}keys%hCluster2) {
			if(!exists $repeat2{$f5}){
				$count++;
				print OUT "$count\t";
				foreach my $f6 (sort keys%{$hCluster2{$f5}}) {
					print OUT "$f6,";
				}
				print OUT "\n";
			}
		}
		close OUT;
	}
	##sub MacroSystenic_step2.pl
	sub MacroSystenic_step2{
		my $aname2=$_[0];
		my %genes=();
		open IN,"<","$opts{o1}/cluster1/$aname2.cluster1" or die;
		while(<IN>){
			chomp;
			s/,$//;
			my @a=split /\t/;
			my @b=split /,/,$a[1];
			foreach my $f(@b){
				$genes{$a[0]}{$f}=1;
			}
		}
		my %matchs=();
		open IN,"<","$opts{o1}/table/$aname2.cluster1.table" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			$matchs{$a[0]}{$a[1]}=1;
			$matchs{$a[1]}{$a[0]}=1;
		}
		close IN;
		my %hCluster=();
		my %list=();
		my $count=0;
		open IN,"<","$opts{o1}/table/$aname2.cluster1.table" or die;
		while(<IN>){
			chomp;
			my @a=split /\t/;
			if(!exists $list{$a[0]} && !exists $list{$a[1]}){
				$count++;
				$list{$a[0]}=$count;
				$list{$a[1]}=$count;
				$hCluster{$count}{$a[0]}=1;
				$hCluster{$count}{$a[1]}=1;
			}
			if(exists $list{$a[0]} && !exists $list{$a[1]}){
				my $cluid1=$list{$a[0]};
				my $flag1=0;
				foreach my $f1(keys%{$hCluster{$cluid1}}){
					$flag1++ unless(exists $matchs{$f1}{$a[1]});
				}
				if($flag1==0){
					$list{$a[1]}=$cluid1;
					$hCluster{$cluid1}{$a[1]}=1;
				}else{
					$count++;
					$list{$a[1]}=$count;
					$hCluster{$count}{$a[1]}=1;
				}
			}
			if(!exists $list{$a[0]} && exists $list{$a[1]}){
				my $cluid2=$list{$a[1]};
				my $flag2=0;
				foreach my $f2(keys%{$hCluster{$cluid2}}){
					$flag2++ unless(exists $matchs{$f2}{$a[0]});
				}
				if($flag2==0){
					$list{$a[0]}=$cluid2;
					$hCluster{$cluid2}{$a[0]}=1;
				}else{
					$count++;
					$list{$a[0]}=$count;
					$hCluster{$count}{$a[0]}=1;
				}
			}	
		}
		close IN;
		open OUT,">$opts{o1}/single_cluster/$aname2.cluster2";
		my %exi=();
		foreach my $key1(sort{$a<=>$b}keys%hCluster){
			print OUT "$key1\t";
			foreach my $key2(keys%{$hCluster{$key1}}){
				$exi{$key2}=1;
				foreach my $r1(keys%{$genes{$key2}}){
					print OUT "$r1,";
				}
			}
			print OUT "\n";
		}
		foreach my $key3(keys%genes){
			next if (exists $exi{$key3});
			$count++;
			print OUT "$count\t";
			foreach my $key4(keys%{$genes{$key3}}){
				print OUT "$key4,";
			}
			print OUT"\n";
		}
		close OUT;
	}
	caculate_family($addname);
	#perl caculate_family.pl -i ancient.cluster2 -d DMel_PYes_MMus_BFlo_NVec_AQue.cluster2 -o try
	sub caculate_family{
		my $cll=$_[0];
		open I,"<$opts{o1}/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene.families\n"); 
		my $line = "";
		my %ancient_spes_num = ();
		my $total_ancient = 0;
		while($line = <I>){
			my %species_status = ();
			$total_ancient++;
			chomp $line;
			my @lines = split /\t/,$line;
			my $entries = $lines[1];
			my @entries = split /\,/,$entries;
			foreach my $entry (@entries){
				my @entry_infos = split /\_/,$entry;
				my $species = $entry_infos[0];
				$species_status{$species} = 1;
			}
			foreach my $species (keys %species_status){
				$ancient_spes_num{$species}++;
			}
		}
		close I;

		open D,"<$opts{o1}/single_cluster/$cll.cluster2" or die;
		my $line1 = "";
		my %ancient_spes_num1 = ();
		my $total_ancient1 = 0;
		while($line1 = <D>){
			my %species_status1 = ();
			$total_ancient1++;
			chomp $line1;
			my @lines1 = split /\t/,$line1;
			my $entries1 = $lines1[1];
			my @entries1 = split /\,/,$entries1;
			foreach my $entry1 (@entries1){
				my @entry_infos1 = split /\_/,$entry1;
				my $species1 = $entry_infos1[0];
				$species_status1{$species1} = 1;
			}
			foreach my $species1 (keys %species_status1){
				$ancient_spes_num1{$species1}++;
			}
		}
		close D;
		open O,">$opts{o1}/gene_families.situation";
		print O "$total_ancient ancient families in total.\n";
		print O "Species\tA:ancient_gene_families_of_each_species\tB:total_ancient_gene_families\tC:total_gene_families_of_each_species\tA/Bpercentage\tA/Cpercentage\n";
		foreach my $species (sort keys %ancient_spes_num){
			my $percentage_ancient = $ancient_spes_num{$species} / $total_ancient;
			my $percentage_ancient2 = $ancient_spes_num{$species} / $ancient_spes_num1{$species};
			print O "$species\t$ancient_spes_num{$species}\t$total_ancient\t$ancient_spes_num1{$species}\t";
			printf O "%0.3f\t", $percentage_ancient;
			printf O "%0.3f\n", $percentage_ancient2;
		}
		close O;
		open C,">$opts{o1}/family_situation.table";
		print C "Species\tTotal_ancient_gene_families\tTotal_gene_families\tType1\tType2\n";
		foreach my $species (sort keys %ancient_spes_num){
			my $percentage_ancient = $ancient_spes_num{$species} / $total_ancient;
			my $percentage_ancient2 = $ancient_spes_num{$species} / $ancient_spes_num1{$species};
			print C "$species\t";
			printf C "%0.3f\t", $percentage_ancient;
			printf C "%0.3f\t", $percentage_ancient2;
			print C "total_ancient_gene_families\t";
			print C "total_gene_families\n";
		}
		close C;

	}
	#choose_gene_in_ancient_gene_families
	#@gr=split/\_/,$addname;
	my %hg;
	open I,"<$opts{o1}/ancient_gene.families" or die("Could not open $opts{o1}/ancient_gene.families\n"); 
	while(my $a = <I>){
		my @gr1=();my @single_gene_id=();
		chomp $a;
		@gr1 = split /\t/,$a;
		@single_gene_id = split /,/,$gr1[1];
		foreach my $gi (@single_gene_id){
			my @gr2 =();
			@gr2 = split /\_/,$gi;
			my $sp1=$gr2[0];
			if (exists $hg{$sp1}) {
				$hg{$sp1}=$hg{$sp1}."\t".$gi;
			}
			else{
				$hg{$sp1}=$gi;
			}
		}
	}
	close I;

	system("mkdir $opts{o1}/genes-in-ancient_gene_families");

	foreach my $i (sort keys %hg) {
		my $grr=$hg{$i};
		my @gr3=();my %h_seq=();
		@gr3 = split /\t/,$grr;
		open I2,"<$dir1/$i.pep";
		while(my $b = <I2>){
			chomp $b;
			$b=~s/\r//g;
			if ($b=~/>(\S+)/) { 
				$iddd=$1; 
			} 
			else { 
				$h_seq{$iddd} .= $b; 
			} 
		}
		close I2;
		open I,"<$opts{o1}/ids_match/$i.ids.match";
		my %h_xp=();
		while(my $a = <I>){
			chomp $a;my @xp=();
			@xp=split/\t/,$a;
			$h_xp{$xp[1]}=$xp[0];
		}
		close I;
		open O,">$opts{o1}/genes-in-ancient_gene_families/$i.ancient_gene_sequences";
		foreach my $m (@gr3) {
			if (exists $h_xp{$m}) {
				my $mmm=$h_xp{$m};
				if (exists $h_seq{$mmm}) {
					print O ">$mmm\n$h_seq{$mmm}\n";
				}
			}
		}
		close O;
	}



	system("mkdir $opts{o1}/ancient_gene_families.analysis");
	print "The [$opts{o1}/ancient_gene_families.analysis] directory has been created!\n";

	system("mv $opts{o1}/single_cluster/$addname.cluster2 $opts{o1}/ancient_gene_families.analysis/all_gene.families");
	system("mv $opts{o1}/ancient_gene.families $opts{o1}/ancient_gene_families.analysis/ancient_gene.families");
	system("mv $opts{o1}/gene_families.situation $opts{o1}/ancient_gene_families.analysis/gene_families.situation");
	system("mv $opts{o1}/ancient_gene_families.analysis/gene_families.situation $opts{o1}/ancient_gene_families.analysis/gene_families_summary.table");
	system "mv $opts{o1}/genes-in-ancient_gene_families $opts{o1}/ancient_gene_families.analysis/ancient_gene_families.genes";
	system "mv $opts{o1}/family_situation.table $opts{o1}/ancient_gene_families.analysis/family_situation.table";
}

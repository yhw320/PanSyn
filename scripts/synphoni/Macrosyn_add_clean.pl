#!/usr/bin/perl -w
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i1=s","m=s","o1=s","h|help");
if (!(defined $opts{i1} and defined $opts{m} and defined $opts{o1})) {
		die "************************************************\n
Options:
	-i1	The parameter refers to the same parameter [-i1] used in the previously executed script [Macrosyn1.pl]
	-m	Enter the abbreviation for the name of the species that represents the ancestral genome (example: NVec)
	-o1	The parameter refers to the same parameter [-o1] used in the previously executed script [Macrosyn1.pl]
	Optional:
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
Options:
	-i1	The parameter refers to the same parameter [-i1] used in the previously executed script [Macrosyn1.pl]
	-m	Enter the abbreviation for the name of the species that represents the ancestral genome (example: NVec)
	-o1	The parameter refers to the same parameter [-o1] used in the previously executed script [Macrosyn1.pl]
	Optional:
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
####
####
if ($opts{o1} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o1} =~ s/$slash$//;
}
####

my $spename1=$opts{m};
my $dir1=$opts{i1}."/peps";
opendir P,$dir1;
while (my $c=readdir P) {
	if ( $c=~/^(\S+).pep$/) {
		if ($c!~$spename1) {
			my $spename2=$1;
			my $gff1="Ancestors_database/".$spename1."/$spename1.gff";
			my $pathway_gffm="$opts{i1}/$gff1";
			my $pathway_gffn="$opts{i1}/gffs/$spename2.gff";
			my $pathway=$opts{o1};
			my $block_name="Ancestors_database/".$spename1."/$spename1.block";
			my $pathway_b="$opts{i1}/$block_name";
			my $nv_cg=$spename1."_".$spename2;
			#system"mkdir $pathway/$nv_cg";
			system"cat $pathway/blastp/$spename1.blastp $pathway/blastp/$spename2.blastp > $pathway/$nv_cg.blastp";
			get_pairs_from_cluster($nv_cg,$spename1,$spename2);
			#system "perl scripts/get_pairs_from_cluster.2.pl single_cluster/ancient.cluster2 $nv_cg/$nv_cg.blastp $spename1 $spename2 > $nv_cg/$nv_cg.pairs";
			synteny_match_id($spename1,$spename2,$nv_cg);
			#system "perl scripts/AJ-synteny-match-id2.pl -i ids_match/$spename1.ids.match -t ids_match/$spename2.ids.match -q $nv_cg/$nv_cg.pairs -o $nv_cg/$nv_cg.pairs2";
			system "rm $pathway/$nv_cg.blastp";
			sub get_pairs_from_cluster {
				my $one=$_[0];
				my $p1=$_[1];
				my $p2=$_[2];
				my %matchs=();
				open IN,"<","$pathway/$one.blastp" or die;
				while(<IN>){
					chomp;
					my @a=split /\t/;
					next unless (($a[0]=~/^$p1\_/ && $a[1]=~/^$p2\_/) || ($a[0]=~/^$p2\_/ && $a[1]=~/^$p1\_/));
					$matchs{$a[0]}{$a[1]}=$a[11];
				}
				close IN;
				my %list=();
				open IN,"<$pathway/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $pathway/ancient_gene_families.analysis/ancient_gene.families.\n"); 
				open O,"> $pathway/$one.pairs";#NV_GG.pairs
				while(<IN>){
					chomp;
					s/,$//;
					my @a=split /\t/;
					my @b=split /,/,$a[1];
					my @t1=();
					my @t2=();
					foreach my $f1(@b){
						push @t1,$f1 if ($f1=~/^$p1\_/);
						push @t2,$f1 if ($f1=~/^$p2\_/);
					}
					next unless (scalar(@t1)>0 && scalar(@t2)>0);
					my %tmp;
					foreach my $r1(@t1){
						foreach my $r2(@t2){
							if (exists $matchs{$r1}{$r2}){
								push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r1}{$r2}];
							}else{
								push @{$tmp{'hh'}},[$r1,$r2,0]
							}
							if(exists $matchs{$r2}{$r1}){
								push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r2}{$r1}];
							}else{
								push @{$tmp{'hh'}},[$r1,$r2,0]
							}
						}
					} 
					foreach my $k2(keys%tmp){
						@{$tmp{$k2}}=sort{$b->[2]<=>$a->[2]}@{$tmp{$k2}};
						for(my $i=0;$i<@{$tmp{$k2}};$i++){
							my $g1=${$tmp{$k2}}[$i][0];
							my $g2=${$tmp{$k2}}[$i][1];
							if(!exists $list{$g1} && !exists $list{$g2}){
								print O "$g1\t$g2\n";
								$list{$g1}=1;
								$list{$g2}=1;
							}
						}
					}	
				}
				close IN;
				close O;
			}
			#sub synteny-match-id2.pl
			sub synteny_match_id{
				my %h1;my %h2;my $spe_id11=();my $spe_id22=();
				my $spe1=$_[0];
				my $spe2=$_[1];
				my $name=$_[2];
				open I,"< $pathway/ids_match/$spe1.ids.match";#NV.id
				while (my $a=<I>) {
					my $spe_id1=$a;
					$spe_id11.=$spe_id1;
				}
				open T,"< $pathway/ids_match/$spe2.ids.match";#CG.id
				while (my $b=<T>) {
					my  $spe_id2=$b;
					$spe_id22.=$spe_id2;
				}
				my $all_id=$spe_id11.$spe_id22;
				my @kk=split/\n/,$all_id;
				foreach $iii(@kk){
					my @itms=split/\t/,$iii;
					$h1{$itms[1]}=$itms[0];
				}
				open Q,"< $pathway/$name.pairs";#NV_GG.pairs
				open O,"> $pathway/$name.pairs2";#NV_GG.pairs2
				while (my $a=<Q>) {
					chomp $a;
					my @itms=split/\t/,$a;
					if (exists $h1{$itms[0]} and exists $h1{$itms[1]}) {
						print O "$h1{$itms[0]}\t$h1{$itms[1]}\n";
					}
					else {print O "error\t$a\n";}
				}
				close I;close T;close Q;close O;
			}

		}
	}
}
close P;
######################

#######################


if (-d "$opts{o1}/add_species") {
	print "The $opts{o1}/add_species directory already exists!\nThe old directory has been deleted.\n";
	system "rm -r $opts{o1}/add_species";
}


system "mkdir $opts{o1}/add_species";
print "The [$opts{o1}/add_species] directory has been created!\n";

system "mkdir $opts{o1}/add_species/$opts{m}";
print "The [$opts{o1}/add_species/$opts{m}] directory has been created!\n";

system "cat $opts{o1}/*.pairs2 > $opts{o1}/add_species/$opts{m}/all.pairs2";
system "rm $opts{o1}/*.pairs2";
system "rm $opts{o1}/*.pairs";
system "cat $opts{i1}/peps/*.pep > $opts{o1}/add_species/$opts{m}/all.peps";
system "sed -i '\/^[^>]\/s\/>\/\\n>\/g' $opts{o1}/add_species/$opts{m}/all.peps";

#system q(sed -i '/^[^>]/s/>/\\n>/g' $opts{o1}/add_species/$opts{m}/all.peps);

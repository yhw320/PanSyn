#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i1=s","m=s","n=s","o1=s","o3=s","e=s","h|help");
if (!(defined $opts{i1} and defined $opts{m} and defined $opts{n} and defined $opts{o1} and defined $opts{o3})) {
		die "************************************************\n
	-i1	The parameter refers to the same parameter [-i1] used in the previously executed script [Macrosyn1.pl]
	-m	Enter the abbreviation for the name of the species representing the ancestral genome (example: NVec)
	-n	Enter the abbreviation for the name of the interested species (example: HSap)
	-o1	The parameter refers to the same parameter [-o1] used in the previously executed script [Macrosyn1.pl]
	-o3	Full path to the new [outputDir3] directory containing output files
	Optional:
	-e	Protein alignment evalue (default:1e-5)
	-h|-help Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-i1	The parameter refers to the same parameter [-i1] used in the previously executed script [Macrosyn1.pl]
	-m	Enter the abbreviation for the name of the species representing the ancestral genome (example: NVec)
	-n	Enter the abbreviation for the name of the interested species (example: HSap)
	-o1	The parameter refers to the same parameter [-o1] used in the previously executed script [Macrosyn1.pl]
	-o3	Full path to the new [outputDir3] directory containing output files
	Optional:
	-e	Protein alignment evalue (default:1e-5)
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
my $slash;
####
if ($opts{o1} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o1} =~ s/$slash$//;
}
####
my $slash;
####
if ($opts{o3} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o3} =~ s/$slash$//;
}
####

if (!(defined $opts{e})) {
	$opts{e}=1e-5;
}
my $spename1=$opts{m};
my $spename2=$opts{n};
my $pathway=$opts{o1};
my $nv_cg=$spename1."_".$spename2;

system"cat $pathway/blastp/$spename1.blastp $pathway/blastp/$spename2.blastp > $opts{o3}/$nv_cg.blastp";
get_pairs_from_cluster($nv_cg,$spename1,$spename2);

synteny_match_id($spename1,$spename2,$nv_cg);

sub get_pairs_from_cluster {
    my $one=$_[0];
    my $p1=$_[1];
    my $p2=$_[2];
    my %matchs=();
    open IN,"<","$opts{o3}/$one.blastp" or die;
    while(<IN>){
	    chomp;
	    my @a=split /\t/;
	    next unless (($a[0]=~/^$p1\_/ && $a[1]=~/^$p2\_/) || ($a[0]=~/^$p2\_/ && $a[1]=~/^$p1\_/));
	    $matchs{$a[0]}{$a[1]}=$a[10];
    }
    close IN;
    my %list=();my %h_hs_mm=();my %h_mm_hs=();
	open IN,"<$pathway/ancient_gene_families.analysis/ancient_gene.families" or die("Could not open $pathway/ancient_gene_families.analysis/ancient_gene.families.\n"); 
    open O,"> $opts{o3}/$one.pairs.score";#NV_GG.pairs
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
			    if (exists $matchs{$r1} and exists $matchs{$r1}{$r2}){
					if ($matchs{$r1}{$r2}<=$opts{e}) {
						$h_hs_mm{$r1}{$r2}=$matchs{$r1}{$r2};
					}				
				    #push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r1}{$r2}];
			    }#else{
				 #  push @{$tmp{'hh'}},[$r1,$r2,0]
			    #}
			    if(exists $matchs{$r2} and exists $matchs{$r2}{$r1}){
					if ($matchs{$r2}{$r1}<=$opts{e}) {
						$h_mm_hs{$r2}{$r1}=$matchs{$r2}{$r1};
					}
				    #push @{$tmp{'hh'}},[$r1,$r2,$matchs{$r2}{$r1}];
			    }#else{
				  #  push @{$tmp{'hh'}},[$r1,$r2,0]
			    #}
		    }
	    } 
		foreach my $g1 (keys %h_hs_mm) {
			foreach my $g2 (keys %{$h_hs_mm{$g1}}) {
				#print O "$g1\t$g2\t$h_hs_mm{$g1}{$g2}\n";
				$h_f{$g1}{$g2}=$h_hs_mm{$g1}{$g2};
			}
		}
		foreach my $g1 (keys %h_mm_hs) {
			foreach my $g2 (keys %{$h_mm_hs{$g1}}) {
				if (exists $h_hs_mm{$g2} and exists $h_hs_mm{$g2}{$g1}) {
				}
				else{
					#print O "$g1\t$g2\t$h_mm_hs{$g2}{$g1}\n";
					$h_f{$g2}{$g1}=$h_mm_hs{$g1}{$g2};
				}
			}
		}

    }
    close IN;
	foreach my $i1 (keys %h_f) {
		foreach my $i2 (keys %{$h_f{$i1}}) {
			print O "$i1\t$i2\t$h_f{$i1}{$i2}\n"
		}
	}
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

    open Q,"< $opts{o3}/$name.pairs.score";#NV_GG.pairs
    open O,"> $opts{o3}/$name.pairs2.score";#NV_GG.pairs2
    while (my $a=<Q>) {
	    chomp $a;
	    my @itms=split/\t/,$a;
	    if (exists $h1{$itms[0]} and exists $h1{$itms[1]}) {
		    print O "$h1{$itms[0]}\t$h1{$itms[1]}\t$itms[2]\n";
	    }
	    else {print O "error\t$a\n";}
    }
    close I;close T;close Q;close O;
}

#!/usr/bin/perl
use Getopt::Long;
my %opts;

GetOptions(\%opts,"pansyn=s","orth=s","n=s","o=s","b=s","h|help");
if (!( defined $opts{orth} and defined $opts{o} and defined $opts{pansyn})) {
		die "************************************************\n
	-orth	Full path to the [Orthogroups.txt] file
	-o	Full path of the [outputDir_S4] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
	-n	The maxmum number of intervening spacer genes allowed between adjacent genes in a block (default: 5)
	-b	The minimum number of genes contained within a microsynteny block (default: 3)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-orth	Full path to the [Orthogroups.txt] file
	-o	Full path of the [outputDir_S4] directory
	-pansyn Full path the [scripts] provided by PanSyn
	Optional:
	-n	The maxmum number of intervening spacer genes allowed between adjacent genes in a block (default: 5)
	-b	The minimum number of genes contained within a microsynteny block (default: 3)
	-h|-help	Print this help page
		*************************************************\n";
}

my $slash; 



####
if ($opts{o} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{o} =~ s/$slash$//;
}
####
####
if ($opts{pansyn} =~ /(\/)$/) {
    # 存储捕获的结果
    $slash = $1;

    # 删除末尾的 /
    $opts{pansyn} =~ s/$slash$//;
}
####

if (!(defined $opts{n})) {
	$opts{n}=5;
}
if (!(defined $opts{b})) {
	$opts{b}=3;
}


####################################
#system "mkdir $opts{o}/final_result";

system "perl $opts{pansyn}/orthoFinderToOrthogroup.pl $opts{orth} >$opts{o}/new_orthogroups.txt";
#chdir $opts{o};
system "perl $opts{pansyn}/prepMicroSynt.pl $opts{o} \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') $opts{n} $opts{o}/new_orthogroups.txt";

opendir P,$opts{o};
#chdir "$opts{o}/final_result";

while (my $c=readdir P) {
	if ( $c=~/^b.(\S+).fa.command$/ ) {
		my $file=$opts{o}."/$c";
		system "perl run_commands.pl -i $file";
	}
}
close P;

system "perl $opts{pansyn}/makeClusters.pl $opts{o} \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') .$opts{n}.blocks $opts{b} 0.3 0.5 > $opts{o}/$opts{n}.blocks.$opts{b}.syn.synt";

system "python $opts{pansyn}/correct_blocks_coordinates.py $opts{o}/$opts{n}.blocks.$opts{b}.syn.synt \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') > $opts{o}/$opts{n}.blocks.$opts{b}.syn_corrected.synt";



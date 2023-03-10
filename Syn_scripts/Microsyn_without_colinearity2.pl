#!/usr/bin/perl
use Getopt::Long;
my %opts;
GetOptions(\%opts,"orth=s","s=s","n=s","o=s","h|help");
if (!( defined $opts{orth} and defined $opts{s} and defined $opts{o})) {
		die "************************************************\n
	-orth	Full path to the [Orthogroups.txt] file (generated by Orthofinder)
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-n	Intervening genes number (default: 5)
	-h|-help	Print this help page
		*************************************************\n";
}
if (defined $opts{h} or defined $opts{help}) {
		die "************************************************\n
	-orth	Full path to the [Orthogroups.txt] file (generated by Orthofinder)
	-s	Full path to the [syn_scripts] file PanSyn provided
	-o	Full path to the [outputDir] folder containing output files
	Optional:
	-n	Intervening genes number (default: 5)
	-h|-help	Print this help page
		*************************************************\n";
}
if (!(defined $opts{n})) {
	$opts{n}=5;
}

####################################
#system "mkdir $opts{o}/final_result";

system "perl $opts{s}/orthoFinderToOrthogroup.pl $opts{orth} >$opts{o}/new_orthogroups.txt";
chdir $opts{o};
system "perl $opts{s}/prepMicroSynt.pl \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') $opts{n} $opts{o}/new_orthogroups.txt";

opendir P,$opts{o};
#chdir "$opts{o}/final_result";

while (my $c=readdir P) {
	if ( $c=~/^b.(\S+).fa.sh$/ ) {
		my $file=$opts{o}."/$c";
		system "$file";
	}
}
close P;

system "perl $opts{s}/makeClusters.pl \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') .5.blocks 3 0.3 0.5 > 5.blocks.3.syn.synt";

system "python3 $opts{s}/correct_blocks_coordinates.py 5.blocks.3.syn.synt \$(ls $opts{o}/gffs/*.chrom | xargs | sed -e 's/ /,/g') > 5.blocks.3.syn_corrected.synt";

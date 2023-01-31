#!/usr/bin/perl -w
use File::Basename;
use strict;

my @p = sort split(/\,/,$ARGV[0]);
my @s = split(/\,/,$ARGV[1]);

print "$ARGV[0]\n";
print "$ARGV[1]\n";

my $clusdir=$ARGV[2];
my $maxsize=100;
if ($ARGV[3]) { $maxsize=$ARGV[3] }

my $cd=`pwd`;
chomp $cd;

my $script = dirname(__FILE__) . "/findSyntBlocks_CLUS_ORTH3.pl";

my $c=0;
for my $size (@s) {

for my $i (0..($#p-1)) {
 for my $j (($i+1)..$#p) {
  my $p1=$p[$i];
  my $p2=$p[$j];
 my ($p1n)=$p1=~/([^\/]*)\.chrom/;
 my ($p2n)=$p2=~/([^\/]*)\.chrom/;

  #if ($p[$j]>$p[$i]) {
  # $p1=$p[$j];
  # $p2=$p[$i];
  #}
  $c++;
  open(O,">b.$c.fa.sh");
 print O "#!/bin/bash\n";
 print O "#\$ -S /bin/bash\n";
 print O "cd $cd\n";
 print O "perl $script $p1 $p2 $clusdir $size $maxsize > $p1n-$p2n.$size.blocks\n";
 #print O "perl /proj/Simakov/scripts/MICROSYNT/findSyntBlocks_CLUS_ORTH2.pl $p1 $p2 $clusdir $size $chromrand $maxsize $pmap > $p1-$p2.$size.randblocks\n";
 close O;
 `chmod u+x b.$c.fa.sh`;
}
}
}

print STDERR " submit: sbatch --array=1-$c --constraint=array-1core job.sh\n";
open(O,">job.sh");
print O "#!/bin/sh\ncd $cd\n./b.\$\{SLURM_ARRAY_TASK_ID\}.fa.sh\n";
close O;


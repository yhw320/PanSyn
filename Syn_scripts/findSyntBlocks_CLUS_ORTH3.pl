#!/usr/bin/perl -w
use strict;

if ($#ARGV==-1) {
 die(" Usage: [p1] [p2] [clus] [Nmax] [max size per proteome]\n");
}


my $maxsize=$ARGV[4];
my %blocks=();

my $p1=$ARGV[0];
my $p2=$ARGV[1];

my %prot=();
open(I,"<$p1");
while (<I>) {
 chomp;
 my @tmp = split(/\s+/,$_);
 $prot{$tmp[1]}=$p1;
}
close I;

open(I,"<$p2");
while (<I>) {
 chomp;
 my @tmp = split(/\s+/,$_);
 $prot{$tmp[1]}=$p2;
}
close I;


print STDERR "Prot.map loaded\n";


my %clus=();
open(I,"<$ARGV[2]");
while (<I>) {
 chomp;
 my @tmp=split(/\t/,$_);
 my %tmp2=();
 for my $i (2..$#tmp) {
  if (not exists $prot{$tmp[$i]}) { #print STDERR " $tmp[$i] not exist!!!"; 
   next }
  $tmp2{$prot{$tmp[$i]}}++;
 }
 my $proc=1;
 for my $i (keys %tmp2) {
  if ($tmp2{$i}>$maxsize) { $proc=0 }
 }
 if ((exists $tmp2{$p1})&&(exists $tmp2{$p2})&&($proc==1)) {
  for my $i (2..$#tmp) { $clus{$tmp[$i]}=1 }
 }

}

my %spec1=load_chr($p1);
my %spec2=load_chr($p2);
my $NMAX=$ARGV[3];

print STDERR "Analyzing blocks\n";
open(I,"<$ARGV[2]");
my $c=0;
my %seen=();
while (<I>) {
 chomp;
 my @tmp=split(/\t/,$_);
 my @seq1=();my @seq2=();
 my %tmp2=();
 for my $i (2..$#tmp) {
  if (exists $spec1{$tmp[$i]}) {push @seq1,$tmp[$i]}
  if (exists $spec2{$tmp[$i]}) {push @seq2,$tmp[$i]}
  if (not exists $prot{$tmp[$i]}) { 
  # print STDERR " $tmp[$i] not exist!!!"; 
  next }
  $tmp2{$prot{$tmp[$i]}}++;
 }
 my $proc=1;
 for my $i (keys %tmp2) {
  if ($tmp2{$i}>$maxsize) { $proc=0 }
 }


if (($#seq1>=0)&&($#seq2>=0)&&($proc==1)) {
for my $seqx (@seq1) {
for my $seqy (@seq2) {
 my $xChr=${$spec1{$seqx}}[0];
 my $xNum=${$spec1{$seqx}}[1];
 my $yChr=${$spec2{$seqy}}[0];
 my $yNum=${$spec2{$seqy}}[1];
 #print "BLOCKS ".keys(%blocks).". Now: $tmp[0]:$xChr, $xNum \t $tmp[1]:$yChr, $yNum\n";
 my %add_h=(); 
 for my $i (keys(%blocks)) {
  my @t1=@{$blocks{$i}{$p1}};
  my @t2=@{$blocks{$i}{$p2}};
  my $t1Chr=$t1[0][1];
  my $t2Chr=$t2[0][1];
  my $ok1=0;my $ok2=0;
  if (($t1Chr eq $xChr)&&($t2Chr eq $yChr)) {
  for my $x (0..$#t1) {
   if (($xNum>=($t1[$x][2]-$NMAX))&&($xNum<=($t1[$x][2]+$NMAX))) {
    $ok1=1;
   }
  }
  for my $y (0..$#t2) {
   if (($yNum>=($t2[$y][2]-$NMAX))&&($yNum<=($t2[$y][2]+$NMAX)) ) {
    $ok2=1;
   }
  }
  if (($ok1==1)&&($ok2==1)) { $add_h{$i}=1; }
  }
 }
 my @add = keys %add_h;
 if ($#add>=0) {
  if ($#add==0) {
  } else {
   my %tmpbl=%blocks;
   for my $x (1..$#add) {
    #print STDERR "Merging $add[$x] to $add[0] ...\n";
    my @bl1=@{$blocks{$add[$x]}{$p1}};
    my @bl2=@{$blocks{$add[$x]}{$p2}};
    for my $z (0..$#bl1) {
     if (not exists $seen{$add[0]}{$bl1[$z][0]}) {
     push @{$blocks{$add[0]}{$p1}}, [ ($bl1[$z][0],$bl1[$z][1],$bl1[$z][2]) ]; $seen{$add[0]}{$bl1[$z][0]}=1 }
    }
    for my $z (0..$#bl2) {
     if (not exists $seen{$add[0]}{$bl2[$z][0]}) {
     push @{$blocks{$add[0]}{$p2}}, [ ($bl2[$z][0],$bl2[$z][1],$bl2[$z][2]) ]; $seen{$add[0]}{$bl2[$z][0]}=1 }
    } 
    delete $blocks{$add[$x]};
   }
  }
  if (not exists $seen{$add[0]}{$seqx}) {
  push @{$blocks{$add[0]}{$p1}}, [($seqx,$xChr,$xNum)]; $seen{$add[0]}{$seqx}=1 }
  if (not exists $seen{$add[0]}{$seqy}) {
  push @{$blocks{$add[0]}{$p2}}, [($seqy,$yChr,$yNum)]; $seen{$add[0]}{$seqy}=1 }
 } else {
  $c++;
  if (not exists $seen{$c}{$seqx}) {
  push @{$blocks{$c}{$p1}}, [($seqx,$xChr,$xNum)]; $seen{$c}{$seqx}=1}
  if (not exists $seen{$c}{$seqy}) {
  push @{$blocks{$c}{$p2}}, [($seqy,$yChr,$yNum)]; $seen{$c}{$seqy}=1}
 }

}
}
}
}
close I;

print STDERR "Printing blocks ... \n";
for my $i (keys(%blocks)) {
 if (($#{$blocks{$i}{$p1}}>=1)&&($#{$blocks{$i}{$p2}}>=1)) {
 my ($p1n)=$p1=~/([^\/]*)\.chrom/;
 my ($p2n)=$p2=~/([^\/]*)\.chrom/;
 my $gn1=$#{$blocks{$i}{$p1}}+1;
 my $gn2=$#{$blocks{$i}{$p2}}+1;
 my $gn=$gn1; if ($gn2<$gn) { $gn=$gn2 }
 print "$i\t$gn\t".$p1n."\t";
 my @t1=@{$blocks{$i}{$p1}};
 my @t2=@{$blocks{$i}{$p2}};
 print "$t1[0][1]\t$p2n\t$t2[0][1]\t";
 my @tmp = ();
 for my $x (@t1) {
  push @tmp, ${$x}[0];
 }
 @tmp = sort { ${$spec1{$a}}[1] <=> ${$spec1{$b}}[1] } @tmp;
 print "".join("\t",@tmp)."\t";
 @tmp=();
 for my $x (@t2) {
  push @tmp, ${$x}[0];
 }
 @tmp = sort { ${$spec2{$a}}[1] <=> ${$spec2{$b}}[1] } @tmp;
 print "".join("\t",@tmp)."\n";
 }
}


sub load_chr {
 my $p = $_[0];
 my %counts=();
 open(I,"<$p");
 my %d=();
 while (<I>) {
  chomp;
  my @tmp = split(/\s+/,$_);
#  if ($tmp[0] eq $p) {
#   if (exists $clus{$tmp[1]}) {
    $d{$tmp[2]}{$tmp[1]}=$tmp[4];
    $counts{$tmp[2]}++;
#   }
#  }
 }
 close I;
 my %o=();
 for my $i (keys(%d)) {
  my @k=sort{ $d{$i}{$a} <=> $d{$i}{$b} } keys %{$d{$i}};
  for my $j (0..$#k) {
   #print "#$i\t$k[$j]\t".($j+1)."\n";
   @{$o{$k[$j]}}=($i,($j+1));
  }
 }
 my @cs=sort {$counts{$a} <=> $counts{$b} } keys %counts;
 my @cs_a=();
 my $sum=0;
 for my $xx (@cs) { push @cs_a, $counts{$xx}; $sum+=$counts{$xx}; }
 print STDERR "Orthol/chr: $p\tN=".($#cs+1)."\tmedian=".$cs_a[$#cs_a/2]."\tmean=".($sum/($#cs+1))."\n";
 return %o;
}



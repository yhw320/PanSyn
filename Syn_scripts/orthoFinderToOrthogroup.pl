#! /usr/bin/perl -w 
use strict;

if (@ARGV!=1) {
  die("Usage: $0 ORTHOGROUPSTXT\n".
    "Convert an OrthoFinder Orthogroups.txt file into an ortho file\n");
}

open(FH,$ARGV[0]);

while(<FH>) {
  next if(/^#/);
  my @i = split;
  my $id = shift(@i);
  $id =~ s/:$//;
  @i = map { (split(m/\s+/,$_))[0] } @i;
  print "$id\t".@i."\t".join("\t",@i)."\n";
}

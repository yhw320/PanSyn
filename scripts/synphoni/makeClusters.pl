#!/usr/bin/perl -w
use strict;

if ($#ARGV==-1) {
 die(" Usage: [proteomes] [block extension, e.g. .5.blocks] [min size] [min overlap] [min overlap between spec] [cond file] [[cutoff]]\n");
}

my $path=$ARGV[0];
my $cutoff=1;
if ($ARGV[7]) { $cutoff=$ARGV[7] }
my $MINSPOL=$ARGV[5];
my @cn=('red','green','magenta','cyan','brown','gold','beige','silver','gray','maroon','teal','navy','purple','indianRed','coral','DarkOrange','LightYellow','SteelBlue','CadetBlue','PaleGreen','springGreen','wheat','tan','RosyBrown','Peru','Chocolate','LawnGreen','Crimson','FireBrick','MediumVioletRed','PaleVioletRed');
my @colcode=( [255,0,0], [0,255,0], [0,0,255] ,[255,255,0],[0,255,255],[255,255,255],[255,0,255],[190,190,190],[125,125,125] ,[125,0,0],[0,125,125],[0,0,125],[125,0,125],[205,92,92],[255,177,80],[255,140,0],[255,255,224],[70,130,180],[95,158,160],[152,251,152],[0,255,127],[245,222,179],[210,180,140],[188,143,143],[205,133,63],[210,105,30],[124,252,0],[220,20,60],[178,34,34],[199,21,133],[219,112,147]);

my $MINOVERLAP=$ARGV[4];

my %cond=();
if ($ARGV[6]) {
if (-s $ARGV[6]) {
open(I,"<$ARGV[6]");
while (<I>) {
 chomp;
 my @tmp = split(/\t/,$_);
 for my $cc (1..($#tmp-1)) {
  my @in=split(/\,/,$tmp[$cc]);
  for my $x (@in) { $cond{$tmp[0]}{1}{$cc}{$x}=1 }
 }
 my @out=split(/\,/,$tmp[$#tmp]);
 for my $x (@out) {
 $cond{$tmp[0]}{2}{$x}=1
 }
}
close I;
}
}

my @p = split(/\,/,$ARGV[1]);
my %col=();
for my $i (0..$#p) {
 $col{$p[$i]}=$i;
}
my $ext=$ARGV[2];
my $ms=$ARGV[3];

my %chrom=();
my %spNames=();
my %prot0=();
for my $xx (@p) {
open(I,"<$xx");
while (<I>) {
 chomp;
 my @tmp = split(/\s+/,$_);
 @{$chrom{$tmp[1]}}=($tmp[2],$tmp[4]);
 my ($nn)=$xx=~/([^\/]*)\.chrom/;
 $spNames{$nn}=$nn;
 $prot0{$tmp[1]}=$nn;
}
close I;
}

print STDERR "Chrom, spnames, prot map loaded\n";

print STDERR "Loading pairwise block data ... \n";
my %clus=();
my %clus_a=();
my $curid=0;
my %oc=();
my %prot=();
my %init_graph=();
my @ss=();
my $spec;
for my $i (0..($#p-1)) {
 for my $j (($i+1)..$#p) {
  my ($p1n)=$p[$i]=~/([^\/]*)\.chrom/;
  my ($p2n)=$p[$j]=~/([^\/]*)\.chrom/;

  my $f = "$path/$p1n-$p2n$ext";
  if (!(-e "$f")) {
   $f="$path/$p2n-$p1n$ext";
  }
  if (!(-e $f)) { print STDERR "$f not found\n" }
  open(I,"<$f");
  while (<I>) {
   chomp;
   my @tmp = split(/\t/,$_);
   my %ctmp=();
   for my $x (6..$#tmp) {
    if (exists $prot0{$tmp[$x]}) {
    push @{$ctmp{$prot0{$tmp[$x]}}}, $tmp[$x];
        }
   }
   my @ctmpa=sort{$#{$ctmp{$a}}<=>$#{$ctmp{$b}}} keys %ctmp;
   if ($#ctmpa>=2) { print STDERR "More than 2 species are present: ".join(",",@ctmpa)."in $f\n"; next; }
   if ($#ctmpa<1) { print STDERR "Only one species present: ".join(",",@ctmpa)." in $f\n"; next; }
   my $n=$#{$ctmp{$ctmpa[0]}}+1;
   my @s1s=@{$ctmp{$ctmpa[0]}}; #sp 1 seq in that block
   my @s2s=@{$ctmp{$ctmpa[1]}}; #sp 2 seq in that block
   my $s1=$prot0{$s1s[0]};
   my $s2=$prot0{$s2s[0]};
   if ($n>=$ms) {
    for my $x (@s1s) {
     for my $y (@s2s) {
      $init_graph{$x}{$y}=1; #connect seq
      $init_graph{$y}{$x}=1;
     }
    }

    %oc=();
    for my $x (@s1s) {
     if (exists $clus{$x}) {
      for my $kk (keys %{$clus{$x}}) { $oc{$kk}++ }
     }
    }
    @ss=@s1s;
    $spec=$s1;
    my $nci1=merge();

    %oc=();
    for my $x (@s2s) {
     if (exists $clus{$x}) {
      for my $kk (keys %{$clus{$x}}) { $oc{$kk}++ }
     }
    }
    @ss=@s2s;
    $spec=$s2;
    my $nci2=merge();
   }
  }
  close I;
 }
}

print STDERR "Making graph ... \n";
my %graph=();
my %sca=();
#for my $x (keys(%init_graph)) {
# for my $y (keys(%{$init_graph{$x}})) {
#  for my $xx (keys %{$clus{$x}}) {
#  for my $yy (keys %{$clus{$y}}) {
#  $graph{$xx}{$yy}=1;
#  $graph{$yy}{$xx}=1;
#  $sca{$prot{$xx}."_".${$chrom{$x}}[0]}{$xx}=1;
#  $sca{$prot{$yy}."_".${$chrom{$y}}[0]}{$yy}=1;
#  }
#  }
# }
#}

print STDERR "Total : ".keys(%clus_a)." clusters\n";

for my $xx (keys %clus_a) {
 for my $yy (keys %clus_a) {
  if ((not exists $graph{$xx}{$yy})&&($xx<$yy)) {
   my @seq1=@{$clus_a{$xx}};
   my @seq2=@{$clus_a{$yy}};
   my $initol=0; my $inittotal=0;
   my %lcx=(); my %lcy=();
   for my $x (@seq1) {
    for my $y (@seq2) {
     if (exists $init_graph{$x}{$y}) { $initol++; $lcx{$x}++; $lcy{$y}++; }
      #$sca{$prot{$xx}."_".${$chrom{$x}}[0]}{$xx}=1;
      #$sca{$prot{$yy}."_".${$chrom{$y}}[0]}{$yy}=1;
     $inittotal++;
    }
   }
   if ((keys(%lcx)>=(($#seq1+1)*$MINSPOL))&&(keys(%lcy)>=(($#seq2+1)*$MINSPOL))) {
    $graph{$xx}{$yy}=1;
    $graph{$yy}{$xx}=1;
   }
  }
 }
}

print STDERR "Preparing graph ... \n";
my %bs=();
my %pos_g=();
my %sums=();
for my $i (keys(%graph)) {
  my @seq=@{$clus_a{$i}};
  my $sp=$prot{$i};
  my %p=();
  for my $j (@seq) {
   $p{$j}=${$chrom{$j}}[1];
  }
  my @sp=sort {$p{$a} <=> $p{$b}} keys %p;
  $bs{$sp}{$i}=($p{$sp[$#sp]}-$p{$sp[0]});
  $pos_g{$i}=$p{$sp[0]};
  $sums{$sp}+=$bs{$sp}{$i};
}

#my %av=();
for my $i (@p) {
 #$av{$i}=$sums{$i}/keys(%{$bs{$i}});
}

my $skip=1;
my ($fn)=$ARGV[2]=~/\.(.*)/;
if ($skip==0) {
print STDERR "Printing .graph and .clans ... \n";
open(O,">$path/$fn.graph");
open(P,">$path/$fn.clans");
print P "sequences=".keys(%graph)."\n";
#print P "<param>\n</param>\n";
print P "<rotmtx>\n1.0;0.0;0.0;\n0.0;1.0;0.0;\n0.0;0.0;1.0\n</rotmtx>\n";
print O "graph G {\n";
print O " overlap=false\n";
my %idMap=();
my $idtmp=0;
my %tmpGroup=();
print P "<seq>\n";
for my $i (keys(%graph)) {
 print O " $i [label=\"$i \(".($#{$clus_a{$i}}+1)."\)\",style=filled, color=$cn[$col{$prot{$i}}]  , fillcolor = $cn[$col{$prot{$i}}] ];\n";
 print P ">$spNames{$prot{$i}}::Block $i\nA\n";
 push @{$tmpGroup{$spNames{$prot{$i}}}}, $idtmp;
 $idMap{$i}=$idtmp;
 $idtmp++;
}
print P "</seq>\n";
print P "<seqgroups>\n";
my $idtmp2=0;
for my $i (keys (%tmpGroup)) {
 print P "name=$i\n";
 print P "type=0\n";
 print P "size=5\n";
 print P "color=$colcode[$idtmp2][0];$colcode[$idtmp2][1];$colcode[$idtmp2][2]\n";
 print P "numbers=".join(";",@{$tmpGroup{$i}})."\n";
 $idtmp2++;
}
print P "</seqgroups>\n";
print P "<pos>\n";
$idtmp=0;
for my $i (keys (%graph)) {
 print P "$idtmp 0 0 0\n";
 $idtmp++;
}
print P "</pos>\n";
print P "<hsp>\n";
for my $i (keys(%graph)) {
 my @id=keys %{$graph{$i}};
 for my $j (@id) {
   if ($i<$j) {
   #print O "  c".$i."r".int($bs{$prot{$i}}{$i}/$av{$prot{$i}}*100)." -- c".$j."r".int($bs{$prot{$j}}{$j}/$av{$prot{$j}}*100).";\n";
   print O "  $i=3 -- $j [len=2,style=bold];\n";
   #determine overlap
   my @seq1=@{$clus_a{$i}};
   my @seq2=@{$clus_a{$j}};
   my $overlap1=0;
   my $overlap2=0;
   for my $xx (@seq1) {
    my $overlap_var=0;
    for my $xxx (keys %{$init_graph{$xx}}) {
     if (exists $clus{$xxx}{$j}) { $overlap_var=1 }
    }
    if ($overlap_var==1) { $overlap1++ }
   }
   for my $xx (@seq2) {
    my $overlap_var=0;
    for my $xxx (keys %{$init_graph{$xx}}) {
     if (exists $clus{$xxx}{$i}) { $overlap_var=1 }
    }
    if ($overlap_var==1) { $overlap2++ }
   }
   my $ratio=$overlap1/($#seq1+1);
   if (($overlap1/($#seq1+1))==0) { print STDERR "$i to $j o1: ($overlap1/($#seq1+1))\n"}
   if (($overlap2/($#seq2+1))==0) { print STDERR "$i to $j o2: ($overlap2/($#seq2+1))\n"}
   if ($ratio>($overlap2/($#seq2+1))) { $ratio=$overlap2/($#seq2+1) }
   $ratio=100*$ratio;
   #end determine overlap
   #print P "$idMap{$i} $idMap{$j}:1E-".(int($ratio))."\n";
   print P "$idMap{$i} $idMap{$j}:1E-20\n";
   }
 }
}
for my $i (keys(%sca)) {
 if (keys(%{$sca{$i}})>1) {
 print O "/* $i*/\n";
 my @tmp=keys(%{$sca{$i}});
 my @seqid= sort {$pos_g{$a} <=> $pos_g{$b}} @tmp;
 print O "  ".join(" -- ",@seqid)." [len=1,style=bold];\n";
 for my $j (0..($#seqid-1)) {
  #print P "$idMap{$seqid[$j]} $idMap{$seqid[$j+1]}:1E-10\n";
 }
 }
}
print O "}\n";
close O;
print P "</hsp>\n";
close P;
}

print STDERR "Printing clusters ... \n";
my $clus_id=0;
#build clusters
open(C,">$path/$fn.$ms.syn.clusters");
my %seenclus=();
for my $i (keys %clus_a) {
 if (exists $seenclus{$i}) {
 } else {
  my %tocheck=%{$graph{$i}};
  my %checked=();
  $clus_id++;
  my %clus_out=();
  $clus_out{$i}=$clus_id;
  while (keys(%tocheck)>0) {
   for my $x (keys %tocheck) {
    if (not exists $checked{$x}) {
     $clus_out{$x}=1;
     for my $y (keys %{$graph{$x}}) {
      $tocheck{$y}=1;
     }
     $checked{$x}=1;
    }
    delete $tocheck{$x};
   }
  }
  for my $x (keys %clus_out) { $seenclus{$x}=$clus_id}
  print C "$clus_id\t".join("\t",keys(%clus_out))."\n";
 }
}
close C;

for my $i (keys(%clus_a)) {
 #determine all proteomes in the syntenic supercluster
 my %unseencon=%{$graph{$i}};
 my %supspec=();
 $supspec{$prot{$i}}=1;
 my %seencon=();
 $seencon{$i}=1;
 while (keys(%unseencon)>0) {
  for my $x (keys(%unseencon)) {
   if (not exists $seencon{$x}) {
   $supspec{$prot{$x}}++;
   $seencon{$x}=1;
   for my $y (keys(%{$graph{$x}})) {
    $unseencon{$y}=1;
   }
   }
   delete($unseencon{$x});
  }
 }
 my @classif;
 for my $c (keys(%cond)) {
 my %cin=(); my $out=0;
 my $incon=keys(%{$cond{$c}{1}});
 for my $x (keys(%supspec)) {
  for my $xx (1..$incon) {
  if (exists $cond{$c}{1}{$xx}{$x}) {$cin{$xx}++}
  }
  if (exists $cond{$c}{2}{$x}) {$out++}
 }
 my $in=1;
 my $miningroup=0;
 for my $xx (1..$incon) { if (not exists $cin{$xx}) { } else { $in++; if ($cin{$xx}>=$cutoff) { $miningroup++ } } }
 if ((keys(%cin)>=2)&&($out==0)&&($miningroup>=2)) { push @classif, $c }
 }
 if ($#classif==-1) { push @classif, "-" }
 #end
 my @seq=@{$clus_a{$i}};
 my @sort_seq=sort { ${$chrom{$a}}[1] <=> ${$chrom{$b}}[1] } @seq;
 if (!(${$chrom{$sort_seq[$#sort_seq]}}[0] eq ${$chrom{$sort_seq[0]}}[0])) {
  print STDERR "WTF not same chrom: $sort_seq[$#sort_seq] != $sort_seq[0] :: ${$chrom{$sort_seq[$#sort_seq]}}[0] vs ${$chrom{$sort_seq[0]}}[0]\n";
 }
 my $sizent=${$chrom{$sort_seq[$#sort_seq]}}[1]-${$chrom{$sort_seq[0]}}[1];
 my $sp=$prot{$i};
 my @maps=();
 my %prots=();
 my $con=keys(%{$graph{$i}});
 for my $x (keys(%{$graph{$i}})) {
  if (exists $clus_a{$x}) { push @maps, "$x\($spNames{$prot{$x}}\)"; $prots{$spNames{$prot{$x}}}=1; }
  else { push @maps, "!".$x }
 }
 my $locmin=1000000000000;
 my $locmax=0;
 my $scf="";
 for my $zzz (@seq) {
  $scf=${$chrom{$zzz}}[0];
  if (${$chrom{$zzz}}[1]<$locmin) { $locmin=${$chrom{$zzz}}[1] }
  if (${$chrom{$zzz}}[1]>$locmax) { $locmax=${$chrom{$zzz}}[1] }
 }
 print "$i\t$spNames{$sp}\t$con\t".join(",",@maps)."\t".keys(%prots)."\t".join(",",keys(%supspec))."\t".join(",",@classif)."\t$scf:$locmin..$locmax\t$sizent\t".join(",",@seq)."\n";
}



sub merge {
my $oid=-1;
my @oldclus=();
for my $kk (keys %oc) {
 if ($oc{$kk}>=$MINOVERLAP) {
  push @oldclus, $kk;
 }
}
if ($#oldclus==-1) {
 $curid++;
 for my $x (@ss) {
  #if (not exists $clus{$x}) {push @{$clus_a{$curid}}, $x }
  push @{$clus_a{$curid}}, $x;
  $clus{$x}{$curid}=1;
 }
 $oid=$curid;
} else {
 $oid=$oldclus[0];
 for my $x (@ss) {
  if (not exists $clus{$x}{$oldclus[0]}) { push @{$clus_a{$oldclus[0]}},$x } else { }
  $clus{$x}{$oldclus[0]}=1;
 }
 if ($#oldclus>0) {
  for my $tmpi (1..$#oldclus) {
   for my $x (@{$clus_a{$oldclus[$tmpi]}}) {
    if (not exists $clus{$x}{$oldclus[0]}) {
    push @{$clus_a{$oldclus[0]}},$x;
    $clus{$x}{$oldclus[0]}=1;
    }
    if (exists $clus{$x}{$oldclus[$tmpi]} ) {
    delete $clus{$x}{$oldclus[$tmpi]};
    }
   }
   delete $clus_a{$oldclus[$tmpi]};
  }
 }
}
if ($oid==-1) {die("Oid=$oid\n") }
$prot{$oid}=$spec;
return $oid;
}
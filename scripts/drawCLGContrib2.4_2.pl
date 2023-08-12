#!/usr/bin/perl -w
use strict;
use warnings;
use SVG;
use Color::Mix;
use Color::Rgb;

if ($#ARGV==-1) {
 die(" Usage: [msynt,chrorder,[[param=value:param=value] 
	msynt is required
	chrorder can be empty
	other parameters are optional and plot various data, separated by \":\" (individual sub-paramaters separated by \",\"

	E.g., perl $0 bfl.msynt:bfl.chr.order[,horiz,selpat]:type=scale:type=alg,file=bil.fam,width=30 human.msynt::type=scale:type=alg,file=bil.fam,width=30

	Main parameter:
	type = segment (default), heatmap, alg, cent, scale

	Data file (if required):
	file = coordinate or data file

	General plotting:
	width = width of track, pixel
	chrcol = chromosome column, 1
        stcol = start column, 2
        spcol = stop column, 3
	valcol = value column, 4 (if -1/-/+/+1 the considered as direction)

	For segment:
	colors = list of colors (separated by -), blue (default)

	For heatmap:
	maxval = maximum value (for heatmap), 1
	minval = minimum value (for heatmap), 0
		(value - minval)/maxval

	For alg:
	window = window size, default 10
	slide = slide size, default 5
	maxbreak = max distance (bp) between orthologs, default 5000000

	For scale:
	ticks = tick separation, bp (default 5000000)
	ticklab = 1 for yes, 0 for no (default 0)
	scalefirst = 1 for scale only for first chromosome, 0 for all chromosomes (chromosome, default 1)

 v. Hannah & Oleg

  \n");
}

my %defcolor=();

sub get_distinct_colors {
  use POSIX 'ceil';

  my $n = shift;
  my $discrete = ceil($n ** (1/3));
  my @vals = map 1 - (($_-1) / $discrete), 1 .. $discrete;
  my @colors;
  my ($r, $g, $b) = (0,0,0);

  for my $i (1 .. $n) {
    push @colors, [@vals[$r,$g,$b]];
    if (++$b == $discrete) {
      if (++$g == $discrete) {
        $r = ($r + 1) % $discrete;
        $g = 0;
      }
      $b = 0;
    }
  }   
  return \@colors;
} 

sub get_distinct_colors2 {
 my @colors=(

[0.8941176,0.1019608,0.10980392],
[0.5529412,0.2980392,0.41568627],
[0.2156863,0.4941176,0.72156863],
[0.2549020,0.5882353,0.50588235],
[0.3019608,0.6862745,0.29019608],
[0.4470588,0.4941176,0.46274510],
[0.5960784,0.3058824,0.63921569],
[0.7960784,0.4000000,0.31764706],
[1.0000000,0.4980392,0.00000000],
[1.0000000,0.7490196,0.09803922],
[1.0000000,1.0000000,0.20000000],
[0.8235294,0.6666667,0.17647059],
[0.6509804,0.3372549,0.15686275],
[0.8078431,0.4196078,0.45098039],
[0.9686275,0.5058824,0.74901961],
[0.7843137,0.5529412,0.67450980],
[0.6000000,0.6000000,0.60000000],


  [ 0.600, 0.060, 0.060 ],
 [ 0.800, 0.320, 0.320],
 [1.000, 0.700, 0.700],
 [ 0.600, 0.330, 0.060],
 [0.800, 0.560, 0.320],
 [1.000, 0.850, 0.700],
 [0.420, 0.600, 0.060],
 [0.640, 0.800, 0.320],
 [0.900, 1.000, 0.700],
 [0.060, 0.420, 0.600],
 [0.320, 0.640, 0.800],
 [0.700, 0.900, 1.000],
 [0.150, 0.060, 0.600],
 [0.400, 0.320, 0.800],
 [0.750, 0.700, 1.000],
 [0.600, 0.600, 0.600],
 [0.200, 0.200, 0.200],
 [1.000,0.4313,0.7058],
 [1.000,0.2313,0.4058]
 );
 return \@colors;
}

 my $scheme2=Color::Mix->new;   
 #my @s2 = $scheme2->analogous('f00000', 20, 20);
 my @s2=@{get_distinct_colors2()};


my $wid=20;
my $hei=400;
my $spacer=25;

my $svg = SVG->new(
    width  => 6000,
    height => ($hei+100)*($#ARGV+1)+300,
);


my %b=();
my %bs=();
my $max=0;

my $xloc=100;
my $yloc=100;

my %SEGCOL=();
my %HITSEGPOS=();
my %ORTHOMAP=();

my $VERT=1;

for my $com (@ARGV[0..$#ARGV]) {
my @data = split /\:/,$com;
%b=();
%bs=();
$max=0;
plot_chr(\@data);
if ($VERT==1) {
	$xloc=15*$#data;
$yloc+=$hei+100;
} else {
	$xloc=15;
	$yloc+=300;
}
}

print $svg->xmlify;




sub plot_chr {
 my $FILE=${$_[0]}[0];
 my $ordbool="";
 my $totalwid=0;
 my %cent=();
 my $selpat=".*";
 if (${$_[0]}[1]) { 
	 my @tmpxxx=split /\,/,${$_[0]}[1];
	 for my $xxx (@tmpxxx) {
		 if ($xxx=~/horiz/) { $VERT=0 }
		 if ($xxx=~/selpat=(.*)/) { $selpat=$1 }
		 if ($xxx=~/ordbool=(.*)/) { $ordbool=$1 }
		 if ($xxx=~/algcolor=(.*)/) { 
			 open(I,"<$1");
			 while (<I>) {
				 chomp;
				 my @tmp = split /\t/;
				 $tmp[5]=~s/\"//g;
				 $defcolor{$tmp[0]}=$tmp[5];
			 }
			 $defcolor{"NOALG"}="#e1e1cd";
			 close I;
			 }
	 }
 }
 print STDERR " vertical = $VERT\n";
 my %segbound=();
 my %segopt=();
 my %famasn=();
 my %famasn_og=();

#INIT ADDFEAT
if ($#{$_[0]}>=2) {
 for my $addfeat (2..$#{$_[0]}) {
 my @addaa=split /\,/,${$_[0]}[$addfeat];
 $segopt{$addfeat}{"type"}="segment";
 @{$segopt{$addfeat}{"color"}}=("blue","green");
 $segopt{$addfeat}{"width"}=$wid/2;
 $segopt{$addfeat}{"chrcol"}=1;
 $segopt{$addfeat}{"stcol"}=2;
 $segopt{$addfeat}{"spcol"}=3;
 $segopt{$addfeat}{"valcol"}=4;
 $segopt{$addfeat}{"maxval"}=1;
 $segopt{$addfeat}{"minval"}=0;
 $segopt{$addfeat}{"window"}=10;
 $segopt{$addfeat}{"slide"}=5;
 $segopt{$addfeat}{"scalefirst"}=1;
 $segopt{$addfeat}{"maxbreak"}=5000000;
 $segopt{$addfeat}{"file"}="NA";
 #@{$segopt{$addfeat}{"ticks"}}=(0,0.5,1);
 $segopt{$addfeat}{"ticks"}=5000000;
 $segopt{$addfeat}{"ticklab"}=0;
 $segopt{$addfeat}{"coloralg"}=0;
 for my $aa (@addaa[0..$#addaa]) {
  if ($aa=~/^file=(.*)/) {
   $segopt{$addfeat}{"file"}=$1
  }
  if (($aa=~/^color=(.*)/)||($aa=~/^colors=(.*)/)) {
   @{$segopt{$addfeat}{"color"}}=split /\-/, $1;
  }
  if ($aa=~/^width=(.*)/) {
   $segopt{$addfeat}{"width"}=$1;
  }
  if ($aa=~/^scalefirst=(.*)/) {
   $segopt{$addfeat}{"scalefirst"}=$1;
  }
  if ($aa=~/^maxbreak=(.*)/) {
   $segopt{$addfeat}{"maxbreak"}=$1;
  }
  if ($aa=~/^type=(.*)/) {
   $segopt{$addfeat}{"type"}=$1;
  }
  if ($aa=~/^chrcol=(.*)/) {
    $segopt{$addfeat}{"chrcol"}=$1;
  }
  if ($aa=~/^stcol=(.*)/) {
    $segopt{$addfeat}{"stcol"}=$1;
  }
  if ($aa=~/^spcol=(.*)/) {
    $segopt{$addfeat}{"spcol"}=$1;
  }
  if ($aa=~/^valcol=(.*)/) {
    $segopt{$addfeat}{"valcol"}=$1;
  }
  if ($aa=~/^maxval=(.*)/) {
    $segopt{$addfeat}{"maxval"}=$1;
  }
  if ($aa=~/^minval=(.*)/) {
    $segopt{$addfeat}{"minval"}=$1;
  }
  if ($aa=~/^window=(.*)/) {
    $segopt{$addfeat}{"window"}=$1;
  }
  if ($aa=~/^slide=(.*)/) {
    $segopt{$addfeat}{"slide"}=$1;
  }
  if ($aa=~/^ticks=(.*)/) {
    $segopt{$addfeat}{"ticks"}=$1;
  }
  if ($aa=~/^ticklab=(.*)/) {
    $segopt{$addfeat}{"ticklab"}=$1;
  }
  if ($aa=~/^connect=(.*)/) {
    $segopt{$addfeat}{"connect"}=$1;
  }
  if ($aa=~/^highlight=(.*)/) {
    $segopt{$addfeat}{"highlight"}=$1;
  }
  if ($aa=~/^coloralg=(.*)/) {
    $segopt{$addfeat}{"coloralg"}=$1;
  }




 }
 if (($segopt{$addfeat}{"type"} eq "segment")||($segopt{$addfeat}{"type"} eq "heatmap")) {
  open(I,"<$segopt{$addfeat}{'file'}");
  print STDERR " $segopt{$addfeat}{'file'} \n";
  while (<I>) {
   chomp;
   my @tmp = split /\s+/;
   my ($sid)=$tmp[$segopt{$addfeat}{"chrcol"}-1]=~/^([^\.]*)/;
   push @{$segbound{$sid}{$addfeat}},[@tmp];
  }
  close I;
 }

 if (($segopt{$addfeat}{"type"} eq "alg")) {
  open(I,"<$segopt{$addfeat}{'file'}");
  print STDERR " $segopt{$addfeat}{'file'} \n";
  while (<I>) {
   chomp;
   my @tmp = split /\s+/;
   for my $xxxx (@tmp[2..$#tmp]) { $famasn{$xxxx}=$tmp[1]; $famasn_og{$xxxx}=$tmp[0]; }
  }
  close I;
 }
 if (($segopt{$addfeat}{"type"} eq "cent")) {
  open(I,"<$segopt{$addfeat}{'file'}");
  print STDERR " $segopt{$addfeat}{'file'} \n";
  my $ccid=0;
  while (<I>) {
   chomp;
   my @tmp = split /\s+/;
   $ccid++;
   $cent{$tmp[0]}{$ccid}{1}=$tmp[1];
   $cent{$tmp[0]}{$ccid}{2}=$tmp[2];
  }
  close I;
 }

 print STDERR "$totalwid :: ".$segopt{$addfeat}{"width"}."\n";
$totalwid+=$segopt{$addfeat}{"width"};
 }
}

my %bpos=();
my %allclg=();
open(I,"<$FILE");
while (<I>) {
 chomp;
 my @tmp = split /\s+/;
 my $seg=1;
 my $chr=$tmp[0];
 if ($tmp[0]=~/([^\.]*)\.(.*)/) {
  $chr=$1;
  $seg=$2;
 }
 if ($chr=~/^$selpat/) { } else { next }
 if (exists $famasn{$tmp[1]}) { push @{$bpos{$chr}{$seg}}, [($famasn{$tmp[1]},$tmp[2],$famasn_og{$tmp[1]})]; $allclg{$famasn{$tmp[1]}}=1 } #position in chr / famasn = CLG
 if (not exists $b{$chr}{$seg}{1}) { $b{$chr}{$seg}{1}=1 } else { if ($tmp[2]<$b{$chr}{$seg}{1}) { $b{$chr}{$seg}{1}=$tmp[2] } }
 if (not exists $b{$chr}{$seg}{2}) { $b{$chr}{$seg}{2}=$tmp[2] } else { if ($tmp[2]>$b{$chr}{$seg}{2}) { $b{$chr}{$seg}{2}=$tmp[2] } }
 if (not exists $bs{$chr}) { $bs{$chr}=$tmp[2] } else { if ($tmp[2]>$bs{$chr}) { $bs{$chr}=$tmp[2] } } # chr size
 if ($tmp[2]>$max) { $max=$tmp[2] }
}
close I;

my @chr=sort { 
 my $alpha=0; 
 my $ac=$a; 
 if ($a=~/^\D{1,}(\d+)/) { $ac=$1; } 
 if ($a=~/^\D{1,}$/) {$alpha=1 }
 if ($a=~/^(\d+)\D{1,}/) { $ac=$1; }
 my $bc=$b; 
 if ($b=~/^\D{1,}(\d+)/) { $bc=$1; } 
 if ($b=~/^(\d+)\D{1,}/) { $bc=$1; }
 if ($b=~/^\D{1,}$/) {  if ($alpha==1) { $alpha=3 } else { $alpha=2 } }
 if ($alpha==0) {  return $ac <=> $bc } else { if ($alpha==1) { return 1<=>0 }; if ($alpha==2) { return 0<=>1 }; if ($alpha==3) { return $ac cmp $bc };  }
} keys %bs;

print STDERR "ordbool = $ordbool ; selpat = $selpat\n";
if (($ordbool eq "")||($ordbool eq $selpat)) { } else {
 %bs=();
 my $count=0;
 my %tmpcountx=();
 open(I,"<$ordbool");
 while (<I>) {
  chomp;
  my @tmp = split /\s+/;
  $count++;
  if ($tmp[0]=~/SPACER/) { $tmpcountx{$tmp[0]}++; $tmp[0]=$tmp[0].".$tmpcountx{$tmp[0]}" }
  print STDERR " tmp = $tmp[0]\n";
  $bs{$tmp[0]}=$count;
  if ($tmp[1]) { $b{$tmp[0]}{1}{1}=1; $b{$tmp[0]}{1}{2}=$tmp[1]; }
 }
 @chr=sort {$bs{$a} <=> $bs{$b}} keys %bs;
}

my $chrcount=0;
for my $x (@chr) {
	$chrcount++;
 print STDERR " CHR = $x ... ";
 if ($x=~/SPACER(\d+)/) { $xloc+=$1; next }
 if ($VERT==1) {
 $svg->text( x=>$xloc-10, y=>$yloc-20,transform=>'rotate(-90,'.($xloc-15).','.($yloc-25).') transform(100,0) scale(-1,1)')->cdata("$x");
 } else {
	 $svg->text( x=>$xloc-10, y=>$yloc-$totalwid+10,transform=>'rotate(-90,'.($xloc-15).','.($yloc-25).') transform(100,0) scale(-1,1)')->cdata("$x");
 }
 my $maxsp=0;
 my $lastpos=0;
 #my $minsp=10000000;
 for my $y (keys %{$b{$x}}) {
  my $sp=int $b{$x}{$y}{2}/$max*$hei;
  #if ($sp<$minsp) { $minsp=$sp }
  if ($sp>$maxsp) { $maxsp=$sp; $lastpos=$b{$x}{$y}{2}; }
  for my $z (0..$#{$bpos{$x}{$y}}) {
	  my $famxxx=${$bpos{$x}{$y}}[$z][2];
	  my $posxxx=$xloc+int ${$bpos{$x}{$y}}[$z][1]/$max*$hei;
	  push @{$ORTHOMAP{$FILE}{$famxxx}}, [ ($posxxx, $yloc+25) ];
  }
 }
 print STDERR " maxsp=$maxsp ";
 #$maxsp=$maxsp-$minsp+1;
 
my $xoff=-$wid;
#PLOTTING through ADDFEAT
for my $addfeat (2..$#{$_[0]}) {
 my @tmpcola=@{$segopt{$addfeat}{"color"}};
 if ($segopt{$addfeat}{"type"} eq "heatmap") {
  @tmpcola=@{get_distinct_colors(10)};
 }
 my $tmpwid=$segopt{$addfeat}{"width"};
 if (exists $segbound{$x}{$addfeat}) {
  #print STDERR "$x exists!\n";
  my $tmpcol=0; my $dir=0;
  for my $xxx (0..$#{$segbound{$x}{$addfeat}}) {
  if (${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1]) { 
   if ((${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1] eq "+")||(${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1] eq "1")) {
   $tmpcol=0; $dir=1;
   } else {
   if ((${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1] eq "-")||(${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1] eq "-1")) {
    $tmpcol=1; $dir=1;
   } else {
    if (${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1]=~/\d+\,\d+\,\d+/) {
    $tmpcol="rgb(".${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1].")"; } 
   }
   }
  }

  my $color="";
  if ($segopt{$addfeat}{"type"} eq "segment") {
   if ($tmpcol=~/^\d+$/) { 
    $color=$tmpcola[$tmpcol];
   } else {
    $color=$tmpcol;
   }
  }
  if ($segopt{$addfeat}{"type"} eq "heatmap") {
   $tmpcol=(${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"valcol"}-1]-$segopt{$addfeat}{"minval"})/$segopt{$addfeat}{"maxval"};
  @tmpcola=("#08306b","#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b","#7f0000");
  @tmpcola=("#ebebff","#d8d8ff","#c4c4ff","#b1b1ff","#9d9dff","#8989ff","#7676ff","#6262ff","#4e4eff","#3b3bff","#2727ff","#1414ff");
   $tmpcol=int 12*$tmpcol;
   if ($tmpcol>$#tmpcola) { $tmpcol=$#tmpcola }
   if ($tmpcol<0) { $tmpcol=0 }
   # my @cola=@{$tmpcola[$tmpcol]};
  #$cola[0]=int(255*$cola[0]);
  #$cola[1]=int(255*$cola[1]);
  #$cola[2]=int(255*$cola[2]);
  #$color=join ",", @cola;
  $color="$tmpcola[$tmpcol]";
  }

  my $cury=int ${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"stcol"}-1]/$max*$hei;
  my $cury2=int ${$segbound{$x}{$addfeat}}[$xxx][$segopt{$addfeat}{"spcol"}-1]/$max*$hei;
  my $yh=$cury2-$cury+1;
  my $TPMX=$xloc-$tmpwid-$xoff;
  my $TPMY=$yloc+$cury;
  my $TPMW=$tmpwid;
  my $TPMH=$yh;
  if ($VERT==0) {
   $TPMX=$xloc+$cury;
   $TPMY=$yloc-$tmpwid-$xoff;
   $TPMW=$yh;
   $TPMH=$tmpwid;
  }
  $svg->rectangle(
   x=>$TPMX,
   y=>$TPMY,
   width=>$TPMW,
   height=>$TPMH,
   style=> {
    'fill' => $color,
    'stroke'         => $color,
    'fill-opacity' => 1,
   }
  );
  if ($dir==0) {
   $tmpcol++;
   if ($tmpcol==($#tmpcola+1)) { $tmpcol=0 }
  }
  }

  if ($chrcount==($#chr+1)) {
   $svg->text( x=>$xloc+$maxsp, y=>$yloc-$tmpwid-$xoff+10,transform=>'rotate(-90,'.($xloc+$hei-50).','.($yloc-$tmpwid-$xoff+10).') transform(100,0) scale(-1,1)')->cdata($segopt{$addfeat}{"file"});
  }

 }
 
 if ($segopt{$addfeat}{"type"} eq "scale") {
  if (($chrcount>1)&&($segopt{$addfeat}{"scalefirst"}==1)) { next }
  my $tickpos=0;
  if ($lastpos==0) { next }
  while (($tickpos/$lastpos)<=1) {
   my $ticktext=$tickpos;
   my $ticksuf="";
   if ($tickpos=~/(\d+)(\d{3,3})/) { $ticktext=$1; $ticksuf="K" }
   if ($tickpos=~/(\d+)(\d{6,6})/) { $ticktext=$1; $ticksuf="M" }
   $svg->line(
   x1=>$xloc-$xoff-$tmpwid,
   y1=>$yloc+$tickpos/$lastpos*$maxsp,
   x2=>$xloc-$xoff,
   y2=>$yloc+$tickpos/$lastpos*$maxsp,
   style=> {
    'fill' => "rgb(0,0,0)",
    'stroke' => "rgb(0,0,0)",
   }
   );
   if ($segopt{$addfeat}{"ticklab"}==1) {
   $svg->text( x=>$xloc-$xoff-$tmpwid, y=>$yloc+$tickpos/$lastpos*$maxsp-5)->cdata("$ticktext$ticksuf");
   }
   $tickpos+=$segopt{$addfeat}{"ticks"};
  }
 }

 if ($segopt{$addfeat}{"type"} eq "orthomap") {
	 print STDERR " orthomap";
	 if ($segopt{$addfeat}{"connect"}) {
		 my $connect=$segopt{$addfeat}{"connect"};
		 if (exists $ORTHOMAP{$connect}) {
		  for my $y (keys %{$bpos{$x}}) {
			  for my $xxxx (0..$#{$bpos{$x}{$y}}) {
				  my $orthocol="190,190,190";
				  if ($segopt{$addfeat}{"highlight"}) {
                                   if (${$bpos{$x}{$y}}[$xxxx][0] eq $segopt{$addfeat}{"highlight"}) { 
				     my $coltmp=$defcolor{${$bpos{$x}{$y}}[$xxxx][0]};
                                         $coltmp=~s/\#//g;
                                        my @rgb = map $_ / 255, unpack 'C*', pack 'H*', $coltmp;
                                        $rgb[0]=int(255*$rgb[0]);
                                        $rgb[1]=int(255*$rgb[1]);
                                        $rgb[2]=int(255*$rgb[2]);
                                        $orthocol=join(",",@rgb);
				   }
                          	  }
				  if ($segopt{$addfeat}{"coloralg"} eq "1") {
					 my $coltmp=$defcolor{${$bpos{$x}{$y}}[$xxxx][0]};
					 $coltmp=~s/\#//g;
 					my @rgb = map $_ / 255, unpack 'C*', pack 'H*', $coltmp;
 					$rgb[0]=int(255*$rgb[0]);
 					$rgb[1]=int(255*$rgb[1]);
 					$rgb[2]=int(255*$rgb[2]);
 					$orthocol=join(",",@rgb);
				  }
				  my $famxxx=${$bpos{$x}{$y}}[$xxxx][2];
				  my $posxxx=$xloc + int ${$bpos{$x}{$y}}[$xxxx][1]/$max*$hei;
				  for my $zz (0..$#{$ORTHOMAP{$connect}{$famxxx}}) {
					  my $conx=${$ORTHOMAP{$connect}{$famxxx}}[$zz][0];
					  my $cony=${$ORTHOMAP{$connect}{$famxxx}}[$zz][1];
  $svg->line(
   x1=>$posxxx,
   y1=>$yloc-3*$tmpwid-$xoff,
   x2=>$conx,
   y2=>$cony,  style=> { 'fill' => "rgb($orthocol)",  'stroke' => "rgb($orthocol)" } ); 
				  }
			  }
		  }
		 }
 	}
 }

if ($segopt{$addfeat}{"type"} eq "alg") {
  my @aaalg=sort keys %allclg;


  my $TPMX=$xloc-$tmpwid-$xoff;
  my $TPMY=$yloc;
  my $TPMW=$tmpwid;
  my $TPMH=$maxsp;
  if ($VERT==0) {
   $TPMX=$xloc;
   $TPMY=$yloc-$tmpwid-$xoff;
   $TPMW=$maxsp;
   $TPMH=$tmpwid;
  }

  $svg->rectangle(
   x=>$TPMX,
   y=>$TPMY,
   width=>$TPMW,
   height=>$TPMH,
   style=> {
    'fill' => "rgb(225,225,205)",
    'stroke' => "rgb(225,225,205)",
   }
  );


 my %halg=();
 for my $xalg (0..$#aaalg) { $halg{$xalg}=$xalg; }
  for my $y (keys %{$bpos{$x}}) {
   my @leftline=();
   my @rightline=();
   my @yline=();
   my $win=$segopt{$addfeat}{"window"};
   my $sl=$segopt{$addfeat}{"slide"};
   my $curpos=0;
   if ($#{$bpos{$x}{$y}}<=$win) { $win=$#{$bpos{$x}{$y}} }
   if ($win<1) { next }
   while (($curpos+$win-1)<=$#{$bpos{$x}{$y}}) {
	   if ($VERT==1) {
    push @leftline, ($xloc-$xoff-$tmpwid);
    push @rightline, ($xloc-$xoff-$tmpwid);
    push @yline,$yloc+ ${$bpos{$x}{$y}}[$curpos+int($win/2)][1]/$max*$hei;
    } else {
	    push @leftline, ($yloc-$xoff-$tmpwid);
    push @rightline, ($yloc-$xoff-$tmpwid);
     push @yline,$xloc+ ${$bpos{$x}{$y}}[$curpos+int($win/2)][1]/$max*$hei;
    }
    $curpos+=$sl;
   }
   my %posdisc=();
   my $maxbreak=$segopt{$addfeat}{"maxbreak"};
   #my $stw=${$bpos{$x}{$y}}[$curpos-1][1];
   #my $spw=${$bpos{$x}{$y}}[$curpos+$win-1][1];
   for my $xx (0..$#aaalg) {
   $curpos=0;
   my $wincount=-1;
   if ($#{$bpos{$x}{$y}}<=$win) { $win=$#{$bpos{$x}{$y}} }
   while (($curpos+$win-1)<=$#{$bpos{$x}{$y}}) {
    my %counts=();
    $wincount++;
    for my $z ($curpos..($curpos+$win-1)) {
     $counts{${$bpos{$x}{$y}}[$z][0]}++;
     if ($z>$curpos) {
      my $prevpos0=${$bpos{$x}{$y}}[$z-1][1];
      my $curpos0=${$bpos{$x}{$y}}[$z][1];
      if (($curpos0-$prevpos0+1)>$maxbreak) { $posdisc{$wincount}=1; }
     }
    }
    my $c=0; if (exists $counts{$aaalg[$xx]}) { $c=$counts{$aaalg[$xx]} }
   # print STDERR "$x\t$stw\t$spw\t$aaalg[$xx]\t$wincount\t$c\n";
    $rightline[$wincount]+= $c/($win+1)*($tmpwid+1);
    $curpos+=$sl;
   }
   
#print STDERR " ".join(",",@rightline)."\n";
#print STDERR " ".join(",",@leftline)."\n";
#print STDERR " ".join(",",@yline)."\n";
my @discaa=sort {$a <=> $b } keys %posdisc;

my $base=0;
while ($base<=$#yline) {
 my $ext=$base;
 if (exists $posdisc{$base}) { $base++; next }
 while ((not exists $posdisc{$ext})&&($ext<=$#yline)) {
  $ext++;
 }
 $ext--;
#print STDERR "at $x: $aaalg[$xx]: $base to $ext\n";

my $xv=[@rightline[$base..$ext], reverse @leftline[$base..$ext] ];
my $yv=[@yline[$base..$ext] , reverse @yline[$base..$ext] ];
if ($VERT==0) {
 $yv=[@rightline[$base..$ext], reverse @leftline[$base..$ext] ];
 $xv=[@yline[$base..$ext] , reverse @yline[$base..$ext] ];
}
#print STDERR " xv=".join ",",@{$xv}." yv=$yv\n";
my $points = $svg->get_path(
    x     =>  $xv,
    y     =>  $yv,
    -type =>'polygon'
);
 
my $colid=-1;
my $col="";
if (keys(%defcolor)>1) {
	#	print STDERR " $xx -> $aaalg[$xx] -> $defcolor{$aaalg[$xx]}\n ";
 $col=$defcolor{$aaalg[$xx]};
 $col=~s/\#//g;
 my @rgb = map $_ / 255, unpack 'C*', pack 'H*', $col;
 $rgb[0]=int(255*$rgb[0]);
 $rgb[1]=int(255*$rgb[1]);
 $rgb[2]=int(255*$rgb[2]);
 $col=join(",",@rgb);
} else {
 if ($aaalg[$xx]=~/(\d+)/) {
  ($colid)=$aaalg[$xx]=~/(\d+)/; $colid--; } else {
  my ($tmpid)=$aaalg[$xx]=~/(\w)$/;
  my $tmpid3=-1; my %mapcol=();
  for my $tmpid2 ('A'..'Z') { $tmpid3++; $mapcol{$tmpid2}=$tmpid3; }
  $colid=$mapcol{$tmpid};
 }
 my @cola=@{$s2[$colid]};
  $cola[0]=int(255*$cola[0]);
  $cola[1]=int(255*$cola[1]);
  $cola[2]=int(255*$cola[2]);
  my $col=join ",", @cola;
}
my $poly = $svg->polygon(
    %$points,
    style=> {
    'fill' => "rgb($col)",
    'stroke'         => "rgb($col)",
   }
);
 $base=$ext+1;
}
 @leftline=@rightline;




   }
  }
 #}
 }

  $xoff+=$tmpwid+5;
 }



 for my $y (keys %{$b{$x}}) {
  my $st=int $b{$x}{$y}{1}/$max*$hei;
  $st++;
  my $sp=int $b{$x}{$y}{2}/$max*$hei;
  my $cw=$wid;
  my $curx=$xloc+($wid-$cw)/2;
  my $cury=$yloc+$st;
  my $curw=$cw;
  my $curh=($sp-$st+2);

  my $col="";
  $col="225,225,225";
 
  @{$HITSEGPOS{$x}{$y}}=($curx,$cury,$curw,$curh);
  if ($col eq "225,225,225") { } else {
  $svg->rectangle(
   x=>$curx,
   y=>$cury,
   width=>$curw,
   height=>$curh,
   style=> {
    'fill' => "rgb($col)",
    'stroke'         => "rgb($col)",
   }
  );
 }
 } 
 print STDERR "\n";
 if ($VERT==1) {
 $xloc+=$wid+$spacer+$xoff;
 } else {
	 $xloc+=$maxsp+30;
 }
}
}

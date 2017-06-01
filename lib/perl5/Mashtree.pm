#!/usr/bin/env perl
package Mashtree;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use List::Util qw/shuffle/;

use threads;
use threads::shared;

use lib dirname($INC{"Mashtree.pm"});
use Bio::Matrix::IO;

our @EXPORT_OK = qw(
           logmsg openFastq _truncateFilename distancesToPhylip createTreeFromPhylip sortNames
           @fastqExt @fastaExt @bamExt @vcfExt @richseqExt
           $MASHTREE_VERSION
         );

local $0=basename $0;

######
# CONSTANTS

our $VERSION = "0.14";
our $MASHTREE_VERSION=$VERSION;
our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fas .fa);
our @bamExt=qw(.sorted.bam .bam);
our @vcfExt=qw(.vcf.gz .vcf);
# Richseq extensions were obtained mostly from bioperl under
# the genbank, embl, and swissprot entries, under
# the source for Bio::SeqIO
our @richseqExt=qw(.gb .gbank .genbank .gbk .gbs .gbf .embl .ebl .emb .dat .swiss .sp);

# Helpful things
my $fhStick :shared;  # A thread can only open a fastq file if it has the talking stick.

#################################################
### COMMON SUBS/TOOLS (not object subroutines) ##
#################################################
# Redefine how scripts die
$SIG{'__DIE__'} = sub {
  local $0=basename($0);
  my $e = $_[0] || ""; 
  my $callerSub=(caller(1))[3] || (caller(0))[3] || "UnknownSub";

  $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; 
  die("$0: $callerSub: $e"); 
};
# Centralized logmsg
#sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}
sub logmsg {
  local $0=basename $0;
  my $parentSub=(caller(1))[3] || (caller(0))[3];
  $parentSub=~s/^main:://;

  # Find the thread ID and stringify it
  my $tid=threads->tid;
  $tid=($tid) ? "(TID$tid)" : "";

  my $msg="$0: $parentSub$tid: @_\n";

  print STDERR $msg;
}

# Opens a fastq file in a thread-safe way.
sub openFastq{
  my($fastq,$settings)=@_;

  my $fh;

  lock($fhStick);

  my @fastqExt=qw(.fastq.gz .fastq .fq.gz .fq);
  my($name,$dir,$ext)=fileparse($fastq,@fastqExt);
  if($ext =~/\.gz$/){
    open($fh,"zcat $fastq | ") or die "ERROR: could not open $fastq for reading!: $!";
  } else {
    open($fh,"<",$fastq) or die "ERROR: could not open $fastq for reading!: $!";
  }
  return $fh;
}

# Removes fastq extension, removes directory name,
# truncates to a length, and adds right-padding.
sub _truncateFilename{
  my($file,$settings)=@_;
  $$settings{truncLength}||=255;
  my $name=basename($file,@fastqExt);
  $name=substr($name,0,$$settings{truncLength}); 
  $name.=" " x ($$settings{truncLength}-length($name)); 
  return $name;
}


# 1. Read the mash distances
# 2. Create a phylip file
sub distancesToPhylip{
  my($distances,$outdir,$settings)=@_;

  my $phylip = "$outdir/distances.phylip"; 
  # NOTE: need to regenerate the combined distances each time
  # because I need to allow variation in the input samples.
  #return $phylip if(-e $phylip);

  # The way phylip is, I need to know the genome names
  # a priori
  my %name;
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";
  while(<MASHDIST>){
    next if(/^#/);
    my($name)=split(/\t/,$_);
    $name=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
    $name{_truncateFilename($name,$settings)}=1;
  }
  close MASHDIST;
  my @name=sortNames([keys(%name)],$settings);
  # Index the array
  my $columnIndex=0;
  for(@name){
    $name{$_}=$columnIndex++;
  }

  # Load up the matrix object
  logmsg "Reading the distances file at $distances";
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";
  my $query="UNKNOWN"; # Default ID in case anything goes wrong
  my @m;
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      $query=_truncateFilename($1,$settings);
    } else {
      my ($reference,$distance)=split(/\t/,$_);
      $reference=_truncateFilename($reference,$settings);
      $distance=sprintf("%0.8f",$distance);
      $m[$name{$query}][$name{$reference}]=$distance;
      $m[$name{$reference}][$name{$query}]=$distance;
    }
  }
  close MASHDIST;
  #my $matrixObj=Bio::Matrix::Generic->new(-rownames=>\@name,-colnames=>\@name,-values=>\@m);

  # taking this method from write_matrix in http://cpansearch.perl.org/src/CJFIELDS/BioPerl-1.6.924/Bio/Matrix/IO/phylip.pm
  my $str;
  $str.=(" " x 4) . scalar(@name)."\n";
  for(my $i=0;$i<@name;$i++){
    $str.=$name[$i];
    my $count=0;
    for(my $j=0;$j<@name;$j++){
      if($count < $#name){
        $str.=$m[$i][$j]. "  ";
      } else {
        $str.=$m[$i][$j];
      }
      $count++;
    }
    $str.="\n";
  }
  open(PHYLIP,">",$phylip) or die "ERROR: could not write to $phylip: $!";
  print PHYLIP $str;
  close PHYLIP;
  return $phylip;
}

sub sortNames{
  my($name,$settings)=@_;
  my @sorted;
  if($$settings{'sort-order'} =~ /^(abc|alphabet)$/){
    @sorted=sort { $a cmp $b } @$name;
  } elsif($$settings{'sort-order'}=~/^rand(om)?/){
    @sorted=shuffle(@$name);
  } elsif($$settings{'sort-order'} eq 'input-order'){
    @sorted=@$name;
  } else {
    die "ERROR: I don't understand sort-order $$settings{'sort-order'}";
  }
  return @sorted;
}

# Create tree file with BioPerl
sub createTreeFromPhylip{
  my($phylip,$outdir,$settings)=@_;

  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $matrix   = Bio::Matrix::IO->new(-format=>"phylip", -file=>$phylip)->next_matrix;
  my $treeObj = $dfactory->make_tree($matrix);
  open(TREE,">","$outdir/tree.dnd") or die "ERROR: could not open $outdir/tree.dnd: $!";
  print TREE $treeObj->as_text("newick");
  print TREE "\n";
  close TREE;

  return $treeObj;

}

1; # gotta love how we we return 1 in modules. TRUTH!!!


#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Uses Mash and BioPerl to create a NJ tree based on distances.
# Run this script with -h for help and usage.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;
use File::Basename qw/basename dirname fileparse/;
use File::Copy qw/copy/;

use threads;
use Thread::Queue;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Mashtree qw/logmsg @fastqExt @fastaExt/;
use Bio::Tree::DistanceFactory;
use Bio::Matrix::IO;
use Bio::Tree::Statistics;

local $0=basename $0;
my $fhStick :shared;  # A thread can only open a fastq file if it has the talking stick.

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outmatrix=s tempdir=s numcpus=i genomesize=i mindepth=i truncLength=i warn-on-duplicate kmerlength=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{truncLength}||=250;  # how long a genome name is
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{kmerlength}||=21;
  logmsg "Temporary directory will be $$settings{tempdir}";

  # Mash-specific options
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}||=5;

  die usage() if($$settings{help});

  # "reads" are either fasta assemblies or fastq reads
  my @reads=@ARGV;
  die usage() if(@reads < 2);

  # Check for prereq executables.
  for my $exe(qw(mash)){
    system("$exe -h > /dev/null 2>&1");
    die "ERROR: could not find $exe in your PATH" if $?;
  }

  logmsg "$0 on ".scalar(@reads)." files";

  my $sketches=sketchAll(\@reads,"$$settings{tempdir}",$settings);

  my $distances=mashDistance($sketches,$$settings{tempdir},$settings);

  my $phylip = distancesToPhylip($distances,$$settings{tempdir},$settings);

  my $treeObj = createTree($phylip,$$settings{tempdir},$settings);

  print $treeObj->as_text('newick');
  print "\n";
  
  return 0;
}

# Run mash sketch on everything, multithreaded.
sub sketchAll{
  my($reads,$sketchDir,$settings)=@_;

  mkdir $sketchDir;

  my $readsQ=Thread::Queue->new(@$reads);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashSketch,$sketchDir,$readsQ,$settings);
  }
  
  $readsQ->enqueue(undef) for(@thr);

  my @mshList;
  for(@thr){
    my $mashfiles=$_->join;
    for my $file(@$mashfiles){
      push(@mshList,$file);
    }
  }

  return \@mshList;
}

# Individual mash sketch
sub mashSketch{
  my($sketchDir,$Q,$settings)=@_;

  my @msh;
  while(defined(my $fastq=$Q->dequeue)){
    my($fileName,$filePath,$fileExt)=fileparse($fastq,@fastqExt,@fastaExt);

    # Do different things depending on fastq vs fasta
    my $sketchXopts="";
    if(grep {$_ eq $fileExt} @fastqExt){
      $sketchXopts.="-m $$settings{mindepth} -g $$settings{genomesize} ";
    } elsif(grep {$_ eq $fileExt} @fastaExt) {
      $sketchXopts.=" ";
    } else {
      logmsg "WARNING: I could not understand what kind of file this is by its extension ($fileExt): $fastq";
    }
      
    logmsg "Sketching $fastq";
    my $outPrefix="$sketchDir/".basename($fastq);

    # See if the user already mashed this file locally
    if(-e "$fastq.msh"){
      logmsg "Found locally mashed file $fastq.msh. I will use it.";
      copy("$fastq.msh","$outPrefix.msh");
    }

    if(-e "$outPrefix.msh"){
      logmsg "WARNING: ".basename($fastq)." was already mashed. You need unique filenames for this script. This file will be skipped: $fastq";
    } elsif(-s $fastq < 1){
      logmsg "WARNING: $fastq is a zero byte file. Skipping.";
    } else {
      system("mash sketch -k $$settings{kmerlength} -s 10000 $sketchXopts -o $outPrefix $fastq > /dev/null 2>&1");
      die if $?;
    }

    push(@msh,"$outPrefix.msh");
  }

  return \@msh;
}

# Parallelized mash distance
sub mashDistance{
  my($mshList,$outdir,$settings)=@_;

  # Make a temporary file with one line per mash file.
  # Helps with not running into the max number of command line args.
  my $mshListFilename="$outdir/mshList.txt";
  open(my $mshListFh,">",$mshListFilename) or die "ERROR: could not write to $mshListFilename: $!";
  print $mshListFh $_."\n" for(@$mshList);
  close $mshListFh;

  my $mshQueue=Thread::Queue->new(@$mshList);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$mshQueue,$mshListFilename,$settings);
  }

  $mshQueue->enqueue(undef) for(@thr);

  my $distfile="$outdir/distances.tsv";
  open(DIST,">",$distfile) or die "ERROR: could not open $distfile for writing: $!";
  for(@thr){
    my $distfiles=$_->join;
    for my $file(@$distfiles){
      # Print the contents of each dist file to the
      # main dist file.
      open(ONEDISTFILE,"<",$file) or die "ERROR: could not open $file for reading: $!";
      while(<ONEDISTFILE>){
        print DIST $_;
      }
      close ONEDISTFILE;
    }
  }
  close DIST;

  return $distfile;
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshQueue,$mshList,$settings)=@_;
  my @dist;
  while(defined(my $msh=$mshQueue->dequeue)){
    my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh -l $mshList > $outfile");
    die if $?;

    push(@dist,$outfile);
  }

  return \@dist;
}

# 1. Read the mash distances
# 2. Create a phylip file
sub distancesToPhylip{
  my($distances,$outdir,$settings)=@_;

  my $phylip = "$outdir/distances.phylip"; 
  return $phylip if(-e $phylip);

  logmsg "Reading the distances file at $distances";
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";

  my $id="UNKNOWN"; # Default ID in case anything goes wrong
  my %m; #matrix for distances
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      $id=_truncateFilename($1,$settings);
    } else {
      my @F=split(/\t/,$_);
      $F[0]=_truncateFilename($F[0],$settings);
      $m{$id}{$F[0]}=sprintf("%0.6f",$F[1]);
    }
  }
  close MASHDIST;

  # Create the phylip file.
  # Make the text first so that we can edit it a bit.
  # TODO I should probably make the matrix the bioperl way.
  logmsg "Creating the distance matrix file for fneighbor.";
  my %seenTruncName;
  my $phylipText="";
  my @genome=sort{$a cmp $b} keys(%m);
  for(my $i=0;$i<@genome;$i++){ 
    my $name=_truncateFilename($genome[$i],$settings);
    $phylipText.="$name  "; 
    if($seenTruncName{$name}++){
      
    }
    for(my $j=0;$j<@genome;$j++){
      $phylipText.=$m{$genome[$i]}{$genome[$j]}."  ";
    }
    $phylipText.= "\n";
  }
  $phylipText=~s/  $//gm;

  # Make the phylip file.
  open(PHYLIP,">",$phylip) or die "ERROR: could not open $phylip for writing: $!";
  print PHYLIP "    ".scalar(@genome)."\n";
  print PHYLIP $phylipText;
  close PHYLIP;

  return $phylip;
}

# Create tree file with BioPerl
sub createTree{
  my($phylip,$outdir,$settings)=@_;

  logmsg "Creating a NJ tree with BioPerl";
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $matrix   = Bio::Matrix::IO->new(-format=>"phylip", -file=>$phylip)->next_matrix;
  my $treeObj = $dfactory->make_tree($matrix);
  open(TREE,">","$outdir/tree.dnd") or die "ERROR: could not open $outdir/tree.dnd: $!";
  print TREE $treeObj->as_text("newick");
  print TREE "\n";
  close TREE;

  return $treeObj;

}

#######
# Utils
#######

# Removes fastq extension, removes directory name,
# truncates to a length, and adds right-padding.
sub _truncateFilename{
  my($file,$settings)=@_;
  my $name=basename($file,@fastqExt);
  $name=substr($name,0,$$settings{truncLength}); 
  $name.=" " x ($$settings{truncLength}-length($name)); 
  return $name;
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

# Lifted from List::MoreUtils
sub uniq (@)
{
    my %seen = ();
    my $k;
    my $seen_undef;
    grep { defined $_ ? not $seen{ $k = $_ }++ : not $seen_undef++ } @_;
}


sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz *.fasta > tree.dnd
  NOTE: fasta files are read as assembly files; fastq files
        are read as raw reads. Fastq file can be gzipped.
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.
  --truncLength        250  How many characters to keep in a filename
  --warn-on-duplicate       Warn instead of die when a duplicate
                            genome name is found

  MASH SKETCH OPTIONS
  --genomesize         5000000
  --mindepth           5     
  --kmerlength         21
  "
}

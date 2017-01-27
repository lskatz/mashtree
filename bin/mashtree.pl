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
use Mashtree qw/logmsg @fastqExt @fastaExt _truncateFilename distancesToPhylip createTreeFromPhylip $MASHTREE_VERSION/;
use Mashtree::Db;
use Bio::Tree::DistanceFactory;
use Bio::Matrix::IO;
use Bio::Tree::Statistics;

my %delta :shared=(); # change in amplitude for peak detection, for each fastq
my $scriptDir=dirname $0;
my $dbhLock :shared;  # Use this as a lock so that only one thread writes to the db at a time
local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outmatrix=s tempdir=s numcpus=i genomesize=i mindepth|min-depth=i truncLength=i kmerlength=i sort-order=s sketch-size=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{truncLength}||=250;  # how long a genome name is
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{'sort-order'}||="ABC";

  # Mash-specific options
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}//=5;
  $$settings{kmerlength}||=21;
  $$settings{'sketch-size'}||=10000;

  # Make some settings lowercase
  for(qw(sort-order)){
    $$settings{$_}=lc($$settings{$_});
  }

  die usage() if($$settings{help});

  # "reads" are either fasta assemblies or fastq reads
  my @reads=@ARGV;
  die usage() if(@reads < 2);

  # Check for prereq executables.
  for my $exe(qw(mash)){
    system("$exe -h > /dev/null 2>&1");
    die "ERROR: could not find $exe in your PATH" if $?;
  }

  logmsg "Temporary directory will be $$settings{tempdir}";
  logmsg "$0 on ".scalar(@reads)." files";

  my $sketches=sketchAll(\@reads,"$$settings{tempdir}",$settings);

  my $phylip = mashDistance($sketches,$$settings{tempdir},$settings);

  logmsg "Creating a NJ tree with BioPerl";
  my $treeObj = createTreeFromPhylip($phylip,$$settings{tempdir},$settings);

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
      my $minDepth=determineMinimumDepth($fastq,$$settings{mindepth},$$settings{kmerlength},$settings);
      $sketchXopts.="-m $minDepth -g $$settings{genomesize} ";
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
      system("mash sketch -k $$settings{kmerlength} -s $$settings{'sketch-size'} $sketchXopts -o $outPrefix $fastq  1>&2");
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

  # Instatiate the database and create the table before the threads get to it
  my $mashtreeDbFilename="$outdir/distances.sqlite";
  my $mashtreeDb=Mashtree::Db->new($mashtreeDbFilename);

  my $mshQueue=Thread::Queue->new(@$mshList);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$mshQueue,$mshListFilename,$mashtreeDbFilename,$settings);
  }

  $mshQueue->enqueue(undef) for(@thr);

  logmsg "Joining $$settings{numcpus} threads";
  for(@thr){
    my $distfiles=$_->join;
  }

  my $phylip = "$outdir/distances.phylip";
  logmsg "Converting to phylip format into $phylip";
  open(my $phylipFh, ">", $phylip) or die "ERROR: could not write to $phylip: $!";
  print $phylipFh $mashtreeDb->toString("phylip");
  close $phylipFh;

  return $phylip;
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshQueue,$mshList,$mashtreeDbFilename,$settings)=@_;
  my @distFile;

  my $mashtreeDb=Mashtree::Db->new($mashtreeDbFilename);
  while(defined(my $msh=$mshQueue->dequeue)){
    my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh -l $mshList > $outfile");
    die "ERROR with 'mash dist -t $msh -l $mshList'" if $?;

    push(@distFile,$outfile);
  }

  # Disconnect right away, before the lock, if no records
  # are going to be inserted into the database
  if(@distFile < 1){
    $mashtreeDb->disconnect();
  }

  lock($dbhLock);
  for my $dist(@distFile){
    $mashtreeDb->addDistances($dist);
  }
}

sub determineMinimumDepth{
  my($fastq,$mindepth,$kmerlength,$settings)=@_;

  $delta{$fastq}//=100;
  my $defaultDepth=1; # if no valley is detected

  return $mindepth if($mindepth > 0);

  my $basename=basename($fastq,@fastqExt);
  
  # Run the min abundance finder to find the valleys
  my $minAbundanceTempdir="$$settings{tempdir}/$basename.minAbundance.tmp";
  mkdir $minAbundanceTempdir;
  my $histFile="$$settings{tempdir}/$basename.hist.tsv";
  my @valley=`$scriptDir/min_abundance_finder.pl $fastq --kmer $kmerlength --tempdir $minAbundanceTempdir --valleys --nopeaks --hist $histFile --delta $delta{$fastq}`;
  die "ERROR with min_abundance_finder.pl on $fastq: $!" if($?);
  chomp(@valley);
  # Some cleanup of large files
  unlink $_ for(glob("$minAbundanceTempdir/*"));
  rmdir $minAbundanceTempdir;

  # If there is no valley, return a default value
  if(!defined $valley[1]){
    $delta{$fastq}=int($delta{$fastq}/2);
    if($delta{$fastq} > 10){
      logmsg "Trying again to determine a min depth with delta==$delta{$fastq} on $fastq";
      return determineMinimumDepth(@_);
    }
    logmsg "WARNING: no valleys were found! Reporting minimum kmer coverage as $defaultDepth.";
    return $defaultDepth;
  }
  
  # Discard the header but keep the first line
  my($minKmerCount, $countOfCounts)=split(/\t/,$valley[1]);
  $minKmerCount=1 if(!defined($minKmerCount) || $minKmerCount < 1);

  logmsg "Setting the min depth as $minKmerCount for $fastq (delta==$delta{$fastq})";

  return $minKmerCount;
}

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz *.fasta > tree.dnd
  NOTE: fasta files are read as assembly files; fastq files
        are read as raw reads. Fastq file can be gzipped.
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.

  TREE OPTIONS
  --truncLength        250  How many characters to keep in a filename
  --sort-order         ABC  For neighbor-joining, the sort order can
                            make a difference. Options include:
                            ABC (alphabetical), random, input-order

  MASH SKETCH OPTIONS
  --genomesize         5000000
  --mindepth           5    If mindepth is zero, then it will be
                            chosen in a smart but slower method,
                            to discard lower-abundance kmers.
  --kmerlength         21
  --sketch-size        10000

  Mashtree version $MASHTREE_VERSION
  "
}


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
use File::Copy qw/cp mv/;
use List::Util qw/shuffle/;

use threads;
use Thread::Queue;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Mashtree qw/logmsg @fastqExt @fastaExt createTreeFromPhylip/;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Tree::DistanceFactory;
use Bio::Tree::Statistics;
use Bio::Matrix::IO;

local $0=basename $0;
my $writeStick :shared;  # Only one thread can write at a time

exit main();

sub main{
  my $settings={};
  my @wrapperOptions=qw(help outmatrix=s distance-matrix tempdir=s numcpus=i save-space reps=i);
  my @mashOptions=qw(kmerlength=i truncLength=i genomesize=i mindepth=i warn-on-duplicate sketch-size=i);
  GetOptions($settings,@wrapperOptions,@mashOptions) or die $!;
  $$settings{numcpus}||=1;
  $$settings{reps}||=0;
  $$settings{kmerlength}||=21;
  $$settings{genomesize}||=5000000;
  $$settings{'sketch-size'}||=10000;
  $$settings{mindepth}||=10;
  $$settings{truncLength}||=250;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});

  logmsg "Temporary directory will be $$settings{tempdir}";

  die usage() if($$settings{help});

  # "reads" are either fasta assemblies or fastq reads
  my @reads=@ARGV;
  die usage() if(@reads < 2);
  
  # Figure out options for mashtree.pl
  my $mashtreexopts="--numcpus $$settings{numcpus} ";
  for my $getopt(@mashOptions){
    my $option=$getopt;
    $option=~s/=.+$//; # remove equals sign and after

    # Treat boolean and options with values differently.
    # If the option never had equals sign in it, it's boolean
    # and the value should be stripped.
    my $value=$$settings{$option};
    if($getopt eq $option){
      $value="";
    }
    $mashtreexopts.="--$option $value " if($$settings{$option});
  }
  
  # Some filenames we'll expect
  my $observeddir="$$settings{tempdir}/observed";
  my $obsDistances="$observeddir/distances.phylip";
  my $observedTree="$observeddir/tree.dnd";

  # Make the observed directory and run Mash
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree.pl --tempdir $observeddir $mashtreexopts @reads > $observedTree.tmp && mv $observedTree.tmp $observedTree");
  die if $?;

  # Read the distance matrix file for shuffling in the bootstrap step.
  my $matrixObj=Bio::Matrix::IO->new(-format=>"phylip", -file=>$obsDistances)->next_matrix;

  logmsg "Running $$settings{reps} replicates";
  my @bsTree;
  for(my $i=1;$i<=$$settings{reps}; $i++){
    # While we are running replicates, if the user enters ctrl-C,
    # Mashtree will simply stop running replicates and will
    # move onto the bootstrapped tree.
    local $SIG{INT}=sub{$$settings{reps}=$i; logmsg "Caught ^C. Only ran $i reps."; return 0;};

    my $tempdir="$$settings{tempdir}/rep$i";
    mkdir($tempdir);
    if($i % 100 == 0){
      logmsg "Mashtree replicate $i - $tempdir";
    }

    if(!-e "$tempdir/tree.dnd"){
      # Test the neighbor-joining algorithm by shuffling the input order
      # Make a new shuffled matrix
      my $shuffledMatrixObj=$matrixObj;
      @{$shuffledMatrixObj->{_names}} = shuffle(@{$shuffledMatrixObj->{_names}});
      my $matrixOut=Bio::Matrix::IO->new(-format=>'phylip',-file=>">$tempdir/distances.phylip");
      $matrixOut->write_matrix($shuffledMatrixObj);
      $matrixOut->close; # need to flush the buffer to the file for the next step to work
      
      # Make a tree from the shuffled matrix
      createTreeFromPhylip("$tempdir/distances.phylip",$tempdir,$settings);
    }

    push(@bsTree,Bio::TreeIO->new(-file=>"$tempdir/tree.dnd")->next_tree);
  }
  
  # Combine trees into a bootstrapped tree and write it 
  # to an output file. Then print it to stdout.
  logmsg "Adding bootstraps to tree";
  my $biostat=Bio::Tree::Statistics->new;
  my $guideTree=Bio::TreeIO->new(-file=>"$observeddir/tree.dnd")->next_tree;
  my $bsTree=$biostat->assess_bootstrap(\@bsTree,$guideTree);
  for my $node($bsTree->get_nodes){
    next if($node->is_Leaf);
    my $id=$node->bootstrap||$node->id||0;
    $node->id($id);
  }
  open(my $treeFh,">","$$settings{tempdir}/bstree.dnd") or die "ERROR: could not write to $$settings{tempdir}/bstree.dnd: $!";
  print $treeFh $bsTree->as_text('newick');
  print $treeFh "\n";
  close $treeFh;

  system("cat $$settings{tempdir}/bstree.dnd"); die if $?;

  if($$settings{'distance-matrix'}){
    cp("$$settings{tempdir}/distances.tsv",$$settings{'distance-matrix'});
  }
  
  return 0;
}

#######
# Utils
#######

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz *.fasta > tree.dnd
  NOTE: fasta files are read as assembly files; fastq files
        are read as raw reads. Fastq file can be gzipped.
  MASH TREE OPTIONS
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.
  --truncLength        250  How many characters to keep in a filename
  --warn-on-duplicate       Warn instead of die when a duplicate
                            genome name is found
  MISC OPTIONS
  --distance-matrix    ''   Output file for distance matrix
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping.
  --save-space              Save space in the temporary directory
                            where possible
  MASH SKETCH OPTIONS
  --genomesize         5000000
  --mindepth           10
  --kmerlength         21
  --sketch-size        10000
  "
}


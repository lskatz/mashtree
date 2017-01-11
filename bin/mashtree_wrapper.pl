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
  my @wrapperOptions=qw(help distance-matrix tempdir=s reps=i);
  GetOptions($settings,@wrapperOptions) or die $!;
  $$settings{reps}||=0;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});

  logmsg "Temporary directory will be $$settings{tempdir}";

  die usage() if($$settings{help});
  die usage() if(@ARGV < 1);
  
  ## Catch some options that are not allowed to be passed
  # Tempdir: All mashtree.pl temporary directories will be under the
  # wrapper's tempdir.
  if(grep(/^\-+tempdir$/,@ARGV) || grep(/^\-+t$/,@ARGV)){
    die "ERROR: tempdir was specified for mashtree.pl";
  }

  my $mashOptions=join(" ",@ARGV);
  
  # Some filenames we'll expect
  my $observeddir="$$settings{tempdir}/observed";
  my $obsDistances="$observeddir/distances.phylip";
  my $observedTree="$observeddir/tree.dnd";

  # Make the observed directory and run Mash
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree.pl --tempdir $observeddir $mashOptions > $observedTree.tmp && mv $observedTree.tmp $observedTree");
  die if $?;

  # Read the distance matrix file for shuffling in the bootstrap step.
  my $matrixObj=Bio::Matrix::IO->new(-format=>"phylip", -file=>$obsDistances)->next_matrix;

  logmsg "Running $$settings{reps} replicates";
  my @bsTree;
  # TODO: multithread
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

    # Test the neighbor-joining algorithm by shuffling the input order
    # Make a new shuffled matrix
    my $shuffledMatrixObj=$matrixObj;
    @{$shuffledMatrixObj->{_names}} = shuffle(@{$shuffledMatrixObj->{_names}});
    my $matrixOut=Bio::Matrix::IO->new(-format=>'phylip',-file=>">$tempdir/distances.phylip");
    $matrixOut->write_matrix($shuffledMatrixObj);
    $matrixOut->close; # need to flush the buffer to the file for the next step to work
    
    # Make a tree from the shuffled matrix
    createTreeFromPhylip("$tempdir/distances.phylip",$tempdir,$settings);

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
  my $usage="$0: a wrapper around mashtree.pl.
  Usage: $0 [options] [-- mashtree.pl options] *.fastq.gz *.fasta > tree.dnd
  --distance-matrix    ''   Output file for distance matrix
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping.
  
  MASHTREE.PL OPTIONS:\n".
  # Print the mashtree options starting with numcpus,
  # skipping the tempdir option.
  `mashtree.pl --help 2>&1 | grep -A 999 numcpus | grep -v ^Stopped`;

  return $usage;
}


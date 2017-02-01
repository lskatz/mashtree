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
use Mashtree::Db;
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
  my @wrapperOptions=qw(help outmatrix=s tempdir=s reps=i numcpus=i);
  GetOptions($settings,@wrapperOptions) or die $!;
  $$settings{reps}||=0;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  $$settings{numcpus}||=1;

  logmsg "Temporary directory will be $$settings{tempdir}";

  die usage() if($$settings{help});
  die usage() if(@ARGV < 1);
  
  ## Catch some options that are not allowed to be passed
  # Tempdir: All mashtree.pl temporary directories will be under the
  # wrapper's tempdir.
  if(grep(/^\-+tempdir$/,@ARGV) || grep(/^\-+t$/,@ARGV)){
    die "ERROR: tempdir was specified for mashtree.pl";
  }
  # Numcpus: this needs to be specified in the wrapper and will
  # appropriately be transferred to the mashtree.pl script
  if(grep(/^\-+numcpus$/,@ARGV) || grep(/^\-+n$/,@ARGV)){
    die "ERROR: numcpus was specified for mashtree.pl";
  }
  # Outmatrix: the wrapper script needs to control where
  # the matrix goes because it can only have the outmatrix
  # for the observed run and not the replicates for speed's
  # sake.
  if(grep(/^\-+outmatrix$/,@ARGV) || grep(/^\-+o$/,@ARGV)){
    die "ERROR: outmatrix was specified for mashtree.pl";
  }

  my $mashOptions=join(" ",@ARGV);
  
  # Some filenames we'll expect
  my $observeddir="$$settings{tempdir}/observed";
  my $obsDistances="$observeddir/distances.phylip";
  my $observedTree="$observeddir/tree.dnd";
  my $outmatrix="$observeddir/distances.tsv";

  # Make the observed directory and run Mash
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree.pl --outmatrix $outmatrix --tempdir $observeddir --numcpus $$settings{numcpus} $mashOptions > $observedTree.tmp && mv $observedTree.tmp $observedTree");
  die if $?;

  # Multithreaded reps
  my $repQueue=Thread::Queue->new(1..$$settings{reps});
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&repWorker, "$observeddir/distances.sqlite", $repQueue, $settings);
    $repQueue->enqueue(undef);
  }

  my @bsTree;
  for(@thr){
    my $treeArr=$_->join;
    for(@$treeArr){
      push(@bsTree,Bio::TreeIO->new(-file=>$_)->next_tree);
    }
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

  if($$settings{'outmatrix'}){
    cp($outmatrix,$$settings{'outmatrix'});
  }
  
  return 0;
}

sub repWorker{
  my($dbFile,$repQueue,$settings)=@_;

  my $mashtreeDb=Mashtree::Db->new($dbFile);

  my @bsTree;
  while(defined(my $rep=$repQueue->dequeue())){
    my $tempdir="$$settings{tempdir}/rep$rep";
    mkdir($tempdir);
    if($rep % 100 == 0){
      logmsg "Mashtree replicate $rep - $tempdir";
    }

    # Test the neighbor-joining algorithm by shuffling the input order
    # Make a new shuffled matrix
    my $phylipString=$mashtreeDb->toString("phylip","rand");
    open(my $phylipFh,">","$tempdir/distances.phylip") or die "ERROR: could not open $tempdir/distances.phylip: $!";
    print $phylipFh $phylipString;
    close $phylipFh;

    # Make a tree from the shuffled matrix
    createTreeFromPhylip("$tempdir/distances.phylip",$tempdir,$settings);
    push(@bsTree,"$tempdir/tree.dnd");
  }
  return \@bsTree;
}

#######
# Utils
#######

sub usage{
  my $usage="$0: a wrapper around mashtree.pl.
  Usage: $0 [options] [-- mashtree.pl options] *.fastq.gz *.fasta > tree.dnd
  --outmatrix          ''   Output file for distance matrix
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping.
  --numcpus            1    This will be passed to mashtree.pl and will
                            be used to multithread reps.
  
  MASHTREE.PL OPTIONS:\n".
  # Print the mashtree options starting with numcpus,
  # skipping the tempdir option.
  `mashtree.pl --help 2>&1 | grep -A 999 "TREE OPTIONS" | grep -v ^Stopped`;

  return $usage;
}


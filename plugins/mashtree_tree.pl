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
use POSIX qw/floor/;
use List::Util qw/min max/;
use Scalar::Util qw/looks_like_number/;

use threads;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg @fastqExt @fastaExt @mshExt @richseqExt _truncateFilename createTreeFromPhylip $MASHTREE_VERSION/;
use Mashtree::Db;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);

  my ($db) = @ARGV;
  die usage() if(!$db || $$settings{help});
  
  # Get the distances data into phylip format
  my $mashDb  = Mashtree::Db->new($db);
  my $phylip  = $mashDb->toString('', 'phylip');
  my $phylipFile = "$$settings{tempdir}/distances.phylip";
  open(my $fh, ">", $phylipFile) or die "ERROR: could not write to $phylipFile: $!";
  print $fh $phylip;
  close $fh;

  # Get the distances into a tree using bioperl
  logmsg "Creating tree with BioPerl";
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $matrix   = Bio::Matrix::IO->new(-format=>"phylip", -file=>$phylipFile)->next_matrix;
  my $treeObj = $dfactory->make_tree($matrix);
  print $treeObj->as_text("newick");
  print "\n";

  return 0;
}


sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 [options] mash.sqlite
  --tempdir            ''   
  --numcpus            1    This script uses Perl threads.
  "
}


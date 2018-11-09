#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>

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

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg @fastqExt @fastaExt @mshExt @richseqExt _truncateFilename createTreeFromPhylip $MASHTREE_VERSION/;
use Mashtree::Db;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outtree=s outmatrix=s tempdir=s numcpus=i genomesize=i mindepth|min-depth=i truncLength=i kmerlength=i sort-order=s sketch-size=i version)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);

  my($db)=@ARGV;

  my $mashdb = Mashtree::Db->new($db);
  print Dumper {$mashdb->toString('',"tsv")};

  return 0;
}

sub usage{
  "$0: dump the database
  Usage: $0 [options] db.sqlite
  "
}


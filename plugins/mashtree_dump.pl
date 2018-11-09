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
  GetOptions($settings,qw(help tempdir=s numcpus=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("MASHTREE_DUMP.XXXXXX",CLEANUP=>1,TMPDIR=>1);

  my($db)=@ARGV;
  die usage() if($$settings{help} || !$db);

  my $mashdb = Mashtree::Db->new($db);
  my %distances = $mashdb->toString('',"tsv");
  my @genome = keys(%distances);
  while(my($g1,$distances)=each(%distances)){
    print $g1;
    for my $g2(@genome){
      my $dist;
      if($g1 eq $g2){
        $dist = 0;
      } else {
        $dist = $distances{$g1}{$g2} || $distances{$g2}{$g1} || ".";
      }
      print "\t$dist";
    }
    print "\n";
  }

  return 0;
}

sub usage{
  "$0: dump the database
  Usage: $0 [options] db.sqlite
  "
}


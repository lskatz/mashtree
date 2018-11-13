#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;
use File::Basename qw/basename dirname fileparse/;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg $MASHTREE_VERSION/;
use Mashtree::Db;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("MASHTREE_INIT.XXXXXX",CLEANUP=>1,TMPDIR=>1);

  my($db)=@ARGV;
  die usage() if($$settings{help} || !$db);

  my $mashdb = Mashtree::Db->new($db);

  return 0;
}

sub usage{
  "$0: initialize the database
  Usage: $0 [options] db.sqlite
  "
}


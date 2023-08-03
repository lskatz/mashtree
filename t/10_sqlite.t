#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;
use Data::Dumper;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use lib "$RealBin/../lib/perl5/x86_64-linux-thread-multi/";

use Test::More tests => 5;

use_ok 'Mashtree';
use_ok 'Mashtree::Db';
use Mashtree qw/_truncateFilename/;

$ENV{PATH}="./bin:$ENV{PATH}";
$ENV{PATH}=~s/quicktree//gi; # remove quicktree for now because it produces a diff ordering of trees

# Test for the sqlite3 executable
system("which sqlite3");
ok($?==0, "SQLite executable") or
  diag("Could not find sqlite. Make sure that the executable sqlite3 is set up properly");

my $tempdir=tempdir("testMash.XXXXXX",CLEANUP=>1,TMPDIR=>1);

my $mashDistancesFile="$tempdir/testdistances.txt";
my $sqliteFile="$tempdir/testdb.sqlite";

# Create a large distances file to see if it can handle large inserts
open(my $fh, ">", $mashDistancesFile) or die "ERROR: could not write to $mashDistancesFile: $!";
print $fh '#query test'."\n";
for(my $i=0;$i<20000;$i++){
  print $fh join("\t","genome$i",$i)."\n";
}
close $fh;

# Make the large inserts and test the hashsum of the db

my $db=Mashtree::Db->new($sqliteFile);
my $numInserted=$db->addDistances($mashDistancesFile);
ok($numInserted == 20000, "Added 20k distances to the database");

subtest 'Testing for specific distances' => sub {
  plan tests=>9;
  for my $expectedDistance(0,10,100,256,987,1234,1432){
    my $dist = $db->findDistance("genome$expectedDistance","test");
    is $expectedDistance,$dist, "Expected distance: $expectedDistance";
  }

  my $distances = $db->findDistances("test");
  my $expectedNumDistances=20000;
  is $expectedNumDistances, scalar(keys(%$distances)), "Expected number of distances with Mashtree::Db::findDistances(): $expectedNumDistances";

  is $$distances{"genome654"}, 654, "Specific distance from Mashtree::Db::findDistances()==654";
};


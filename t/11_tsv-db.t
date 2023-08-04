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

use Test::More tests => 4;

use_ok 'Mashtree';
use_ok 'Mashtree::Db';
use Mashtree qw/_truncateFilename/;

$ENV{PATH}="./bin:$ENV{PATH}";
$ENV{PATH}=~s/quicktree//gi; # remove quicktree for now because it produces a diff ordering of trees

my $tempdir=tempdir("testMashDb.XXXXXX",CLEANUP=>1,TMPDIR=>1);

my $mashDistancesFile="$tempdir/testdistances.txt";
my $dbFile="$tempdir/testdb.db.tsv";

# Create a large distances file to see if it can handle large inserts
open(my $fh, ">", $mashDistancesFile) or die "ERROR: could not write to $mashDistancesFile: $!";
# Each genome is labeled as genomeX and is X distance from genome 'test'
print $fh '#query test'."\n";
for(my $i=0;$i<20000;$i++){
  print $fh join("\t","genome$i",$i)."\n";
}
close $fh;

# Make the large inserts and test the hashsum of the db

my $db=Mashtree::Db->new($dbFile);
my $numInserted=$db->addDistances($mashDistancesFile);
ok($numInserted == 20000, "Added 20k distances to the database (actually added: $numInserted)");

subtest 'Testing for specific distances' => sub {
  plan tests=>16;
  # Reload the database internally
  $db->readDatabase;

  # Test for these expected distances one at a time
  for my $expectedDistance(0,10,100,256,987,1234,1432){
    my $dist = $db->findDistance("genome$expectedDistance","test");
    is $expectedDistance,$dist, "Expected distance: $expectedDistance";

    my $distRev = $db->findDistance("test","genome$expectedDistance",);
    is $expectedDistance,$distRev, "Expected distance in reverse: $expectedDistance";
  }

  my $distances = $db->findDistances("test");
  my $expectedNumDistances=20000;
  is $expectedNumDistances, scalar(keys(%$distances)), "Expected number of distances with Mashtree::Db::findDistances(): $expectedNumDistances";

  is $$distances{"genome654"}, 654, "Specific distance from Mashtree::Db::findDistances()==654";
};


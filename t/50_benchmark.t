#!/usr/bin/env perl

use strict;
use warnings;
use Benchmark ':all';
use Test::More tests=>4;
use File::Temp qw/tempdir/;
use Data::Dumper;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use lib "$RealBin/../lib";
use Mashtree::Db;

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

my $db = Mashtree::Db->new($dbFile);
my $numInserted=$db->addDistances($mashDistancesFile);
ok($numInserted == 20000, "Added 20k distances to the database (actually added: $numInserted)");

my $distsHash = $db->readDatabase;
is(keys(%$distsHash), 1, "Keys of dists file total 1");
is(keys(%{ $$distsHash{test} }), 20000, "Number of other genomes compared against 'test' genome is 20k.");

unlink($dbFile);
$db = Mashtree::Db->new($dbFile);
$numInserted=$db->addDistancesFromHash($distsHash);
ok($numInserted == 20000, "Added 20k distances to the database (actually added: $numInserted)");

timethese(30, {
  'addDistancesFromMashFile' => sub{
    unlink($dbFile);
    my $db = Mashtree::Db->new($dbFile);
    $db->addDistances($mashDistancesFile);
  },
  'addDistancesFromHash'     => sub{
    unlink($dbFile);
    my $db = Mashtree::Db->new($dbFile);
    $db->addDistancesFromHash($distsHash);
  },
  'readDatabase'             => sub{
    $db->readDatabase;
  },
});


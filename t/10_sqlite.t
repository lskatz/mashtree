#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;
use Digest::MD5 qw/md5_hex/;
use Data::Dumper;

use Test::More tests => 6;

use_ok 'Mashtree';
use_ok 'Mashtree::Db';
use Mashtree qw/_truncateFilename/;

$ENV{PATH}="./bin:$ENV{PATH}";
$ENV{PATH}=~s/quicktree//gi; # remove quicktree for now because it produces a diff ordering of trees

# Test for the sqlite3 executable
system("which sqlite3");
ok($?==0, "SQLite executable");

my $tempdir=tempdir("testMash.XXXXXX",CLEANUP=>1,TMPDIR=>1);

my $mashDistancesFile="$tempdir/testdistances.txt";
my $sqliteFile="$tempdir/testdb.sqlite";
my $distancesMd5sum='50ad6869eff5a99de5465f3f7e30723c';
my $dbMd5sum='2bf8b5866060828687485809286b3d64';

# Create a large distances file to see if it can handle large inserts
open(my $fh, ">", $mashDistancesFile) or die "ERROR: could not write to $mashDistancesFile: $!";
print $fh '#query test'."\n";
for(my $i=0;$i<20000;$i++){
  print $fh join("\t","genome$i",$i)."\n";
}
close $fh;
ok(md5_hex(`cat $mashDistancesFile`) eq $distancesMd5sum, "Make the distances file for input to the database");

# Make the large inserts and test the hashsum of the db

my $db=Mashtree::Db->new($sqliteFile);
$db->addDistances($mashDistancesFile);
ok(md5_hex(`cat $sqliteFile`) eq $dbMd5sum, "Added distances to the database");

my $distance=$db->findDistance(_truncateFilename("test"),_truncateFilename("genome7009"));
is($distance,7009,"Get distance value from database");

#system("md5sum 

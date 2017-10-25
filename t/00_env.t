#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 2;

# Is Mash installed?
my $mash_is_missing = system("mash >/dev/null 2>&1");
is $mash_is_missing, 0, "Found Mash in PATH";
if($mash_is_missing){
  die "ERROR: the executable Mash was not found in PATH";
}

# test for version
my $version=-1;
open(my $fh, 'mash |') or die "ERROR: could not run mash by itself: $!";
while(<$fh>){
  chomp;
  if(/version\s+(\d+(\.?\d+)?)/){
    $version=$1;
    diag "Found Mash version $version";
  }
}
close $fh;
if($version < 0){
  diag "No version was found when running mash";
}
ok($version >= 2, "Version should be >=2");


#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 3;

use lib './lib';
use_ok 'Mashtree';

# Is Mash installed?
my $mash_is_missing = system("mash >/dev/null 2>&1");
is $mash_is_missing, 0, "Found Mash in PATH";
if($mash_is_missing){
  die "ERROR: the executable Mash was not found in PATH";
}

# Is Quicktree installed?
my $quicktree_is_missing = system("quicktree > /dev/null 2>&1");
   $quicktree_is_missing = $quicktree_is_missing % 256;
is $quicktree_is_missing, 0, "Found quicktree in PATH";
if($quicktree_is_missing){
  die "ERROR: the executable quicktree was not found in PATH";
}


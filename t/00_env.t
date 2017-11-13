#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 2;

use lib './lib';
use_ok 'Mashtree';

# Is Mash installed?
my $mash_is_missing = system("mash >/dev/null 2>&1");
is $mash_is_missing, 0, "Found Mash in PATH";
if($mash_is_missing){
  die "ERROR: the executable Mash was not found in PATH";
}


#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 5;
use FindBin qw/$RealBin/;

use lib './lib';
use_ok 'Mashtree';
$ENV{PATH}="$RealBin/../bin:".$ENV{PATH};

# Is Mash installed?
my $mash_is_missing = system("mash > mash.log 2>&1");
is $mash_is_missing, 0, "Found Mash in PATH" or
  diag("The executable Mash was not found in PATH: ".`cat mash.log`);

END{unlink("mash.log");}

# Is Quicktree installed?
my $quicktree_is_missing = system("quicktree 2>&1 | grep -C 50 -i quicktree > mash.log 2>&1"); # overwrite the mash log
my $quicktree_is_missing_mod = $quicktree_is_missing >> 8; # shift the exit code down 8 bits
is $quicktree_is_missing_mod, 0, "Found quicktree in PATH" or 
  diag("The executable quicktree was not found in PATH: ".`cat mash.log`);

# Test out mashtree exe
my $version = `mashtree --version`;
my $exit_code = $? >> 8;
is $exit_code, 0, "Mashtree --version exit code: $exit_code";

# If mash, quicktree, or mashtree gave an exit code, bail out of the whole test
if($mash_is_missing || $quicktree_is_missing || $exit_code){
  BAIL_OUT("Prerequisite software was not found");
}

$version =~ s/Mashtree\s*//; # Mashtree trim
$version =~ s/^\s+|\s+$//;   # whitespace trim
my $found_nonversion = !! ($version=~/([^\.\d])/) + 0;
is($found_nonversion, 0, "Looking for version numbers in the format of 1.2.3 (Version returned was $version)");


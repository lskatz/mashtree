#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 5;
use FindBin qw/$RealBin/;
use Scalar::Util qw/looks_like_number/;

use lib './lib';
use_ok 'Mashtree';
$ENV{PATH}="$RealBin/../bin:".$ENV{PATH};

# Is Mash installed?
my $mash_is_missing = system("mash >/dev/null 2>&1");
is $mash_is_missing, 0, "Found Mash in PATH" or
  diag("The executable Mash was not found in PATH.");

# Is Quicktree installed?
my $quicktree_is_missing = system("quicktree > /dev/null 2>&1");
   $quicktree_is_missing = $quicktree_is_missing % 256;
is $quicktree_is_missing, 0, "Found quicktree in PATH" or 
  diag("The executable quicktree was not found in PATH.");

# Test out mashtree exe
my $version = `mashtree --version`;
my $exit_code = $? % 256;
is $exit_code, 0, "Mashtree --version exit code: $exit_code";

$version =~ s/[^\d\.]//g;
ok looks_like_number($version), "Version number returned: $version";


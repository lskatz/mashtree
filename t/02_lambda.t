#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 2;

use_ok 'Mashtree';

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="(sample3:0.00195,sample4:0.00205,(sample1:0.00205,sample2:0.00205):0.00010);";

# Test to see if the correct tree is made
my $mashtree=`mashtree.pl --numcpus 1 t/lambda/*.fastq.gz`;
chomp($mashtree);
is $mashtree, $correctMashtree, "Mashtree produced the expected tree"


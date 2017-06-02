#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib/perl5';
use File::Basename qw/dirname/;

use Test::More tests => 2;

use_ok 'Mashtree';

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="(CFSAN000189.gbk:0.00000,(CFSAN001112.ref:0.00001,CFSAN001140_1:0.00019):0.00001,((CFSAN000968.ref:0.00000,CFSAN001115.ref:0.00000):0.00001,(CFSAN000189.ref:0.00000,(CFSAN000191.ref:0.00002,(CFSAN000211.gbk:0.00045,CFSAN000961.gbk:0.00005):0.00003):0.00000):0.00001):0.00000);";

# Test to see if the correct tree is made
my $mashtree=`mashtree.pl --numcpus 1 t/filetypes/*`;
chomp($mashtree);
is $mashtree, $correctMashtree, "Mashtree produced the expected tree"


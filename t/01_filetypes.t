#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib/perl5';
use File::Basename qw/dirname/;

use Test::More tests => 2;

use_ok 'Mashtree';

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="(CFSAN000189.ref.fasta:0.00000,CFSAN001112.ref.fasta:0.00002,((CFSAN000968.ref.fasta:0.00001,(CFSAN001140_1:0.00010,(CFSAN000211.gbk:0.02103,(CFSAN000189.gbk:0.97421,CFSAN000961.gbk:0.02579):0.00302):0.02295):0.00009):0.00002,(CFSAN000191.ref.fasta:0.00000,CFSAN001115.ref.fasta:0.00000):0.00002):0.00000);";

# Test to see if the correct tree is made
my $mashtree=`mashtree.pl --numcpus 1 t/data/*`;
chomp($mashtree);
is $mashtree, $correctMashtree, "Mashtree produced the expected tree"


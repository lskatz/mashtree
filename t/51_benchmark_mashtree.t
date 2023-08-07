#!/usr/bin/env perl

use strict;
use warnings;
use Benchmark ':all';
use Test::More tests=>1;
use File::Temp qw/tempdir/;
use File::Which qw/which/;
use Data::Dumper;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use lib "$RealBin/../lib";

$ENV{PATH}="$RealBin/../bin:$ENV{PATH}";

my $exe = which('mashtree');
ok($exe ne "", "found a path to mashtree: $exe");

timethese(20, {
  'vanilla mashtree' => sub{
    my $dnd = `mashtree --numcpus 1 --genomesize 40000 $RealBin/lambda/*.fastq.gz 2>/dev/null`;
    die if $?;
  },
});
timethese(5, {
  'bootstrap mashtree' => sub{
    my $dnd = `mashtree_bootstrap.pl --reps 10 --numcpus 1 $RealBin/lambda/*.fastq.gz -- --genomesize 40000 2>/dev/null`;
    die if $?;
  },
  'jackknife mashtree' => sub{
    my $dnd = `mashtree_jackknife.pl --reps 10 --numcpus 1 $RealBin/lambda/*.fastq.gz -- --genomesize 40000 2>/dev/null`;
    die if $?;
  },
});


#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';
use File::Basename qw/dirname/;
use Getopt::Long qw/GetOptions/;

use Test::More tests => 2;

use_ok 'Mashtree';
use Mashtree qw/treeDist/;

$ENV{PATH}="./bin:$ENV{PATH}";

my $settings={};
GetOptions($settings,qw(numcpus=i)) or die $!;
$$settings{numcpus}||=1;

my $correctMashtree='((CFSAN000191.ref:0.0000277081,(CFSAN000189.ref:0.0000023811,CFSAN000189.gbk:0.0000000000):0.0000032768):0.0000050668,(CFSAN001115.ref:0.0000099241,CFSAN001112.ref:0.0000210585):0.0000032901,(((CFSAN001140_1:0.0112883298,CFSAN000211.gbk:0.0004680698):0.0000119053,CFSAN000968.ref:0.0000164555):0.0000085410,CFSAN000961.gbk:0.0000418070):0.0000036369);';
$correctMashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

# Test to see if the correct tree is made
my $mashtree=`mashtree --min-depth 0 --numcpus $$settings{numcpus} t/filetypes/*`;
chomp($mashtree);
$mashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

#print STDERR "$mashtree\n";

my $treedist=treeDist($correctMashtree, $mashtree);
ok $treedist < 2, "Correct min_abundance_filter tree (dist: $treedist)";


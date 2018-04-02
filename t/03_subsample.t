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

#my $correctMashtree="(CFSAN001112.ref:0.00001,CFSAN001115.ref:0.00000,((((CFSAN001140_1:0.16779,CFSAN000211.gbk:0.00021):0.00025,CFSAN000191.ref:0.00000):0.00000,CFSAN000189.ref:0.00000):0.00001,((CFSAN000968.ref:0.00000,CFSAN000961.gbk:0.00010):0.00000,CFSAN000189.gbk:0.00000):0.00001):0.00000);";
my $correctMashtree='(CFSAN000189.gbk:0.0000000000,((CFSAN001115.ref:0.0000027557,CFSAN000961.gbk:0.0000544902):0.0000051566,CFSAN001112.ref:0.0000163108):0.0000107139,((((CFSAN001140_1:0.1677920464,CFSAN000211.gbk:0.0002099480):0.0002478605,CFSAN000968.ref:0.0000000000):0.0000227821,CFSAN000191.ref:0.0000131458):0.0000084592,CFSAN000189.ref:0.0000010247):0.0000001422);';

# Test to see if the correct tree is made
my $mashtree=`mashtree --min-depth 0 --numcpus $$settings{numcpus} t/filetypes/*`;
chomp($mashtree);

#print STDERR "$mashtree\n";

my $treedist=treeDist($correctMashtree, $mashtree);
ok $treedist < 1, "Correct min_abundance_filter tree (dist: $treedist)";


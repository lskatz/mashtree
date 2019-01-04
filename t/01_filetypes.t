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

#my $correctMashtree="(CFSAN000189.gbk:0.00000,(CFSAN001112.ref:0.00001,CFSAN001140_1:0.00019):0.00001,((CFSAN000968.ref:0.00000,CFSAN001115.ref:0.00000):0.00001,(CFSAN000189.ref:0.00000,(CFSAN000191.ref:0.00002,(CFSAN000211.gbk:0.00045,CFSAN000961.gbk:0.00005):0.00003):0.00000):0.00001):0.00000);";
my $correctMashtree="(CFSAN001112.ref:0.0000197713,(CFSAN001115.ref:0.0000107516,((CFSAN000968.ref:0.0000186460,(CFSAN000961.gbk:0.0000252712,CFSAN000211.gbk:0.0004753628):0.0000146923):0.0000077003,(CFSAN000191.ref:0.0000276279,(CFSAN000189.ref:0.0000023811,CFSAN000189.gbk:0.0000000000):0.0000033570):0.0000037416):0.0000038797):0.0000011147,CFSAN001140_1:0.0002225407);";
$correctMashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

# Test to see if the correct tree is made
my $mashtree=`mashtree --numcpus $$settings{numcpus} t/filetypes/*`;
chomp($mashtree);
$mashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

ok treeDist($mashtree,$correctMashtree) < 11, "Produce correct tree from mixed file types";


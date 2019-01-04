#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 2;

use_ok 'Mashtree';
use Mashtree qw/treeDist/;

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="((sample2:0.0020443525,sample1:0.0021037373):0.0000540274,sample3:0.0019622177,sample4:0.0020673526);";

# Test to see if the correct tree is made
my $mashtree=`mashtree --file-of-files --numcpus 1 t/file-of-files.txt`;
chomp($mashtree);

is($mashtree, $correctMashtree, "File of files test on lambda");


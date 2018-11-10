#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 4;

use_ok 'Mashtree';
use Mashtree qw/treeDist/;

$ENV{PATH}="./bin:$ENV{PATH}";
$ENV{PATH}="./plugins:$ENV{PATH}";

my $db="mash.sqlite";
END{unlink($db);}

my $correctMashtree="(sample1:0.00210,sample2:0.00204,(sample3:0.00196,sample4:0.00207):0.00005);";

system("mashtree_init.pl $db");
is $?, 0, "Ran mashtree_init.pl";

system("mashtree_mash.pl --numcpus 1 $db t/lambda/sample*.fastq.gz");
is $?, 0, "Ran mashtree_init.pl";

system("mashtree_optimize.pl $db");
is $?, 0, "Ran mashtree_optimize.pl";

system("mashtree_dump.pl $db");
my $tree = `mashtree_tree.pl $db`; chomp($tree);
is $tree, $correctMashtree, "mashtree_tree.pl";



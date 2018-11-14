#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 6;

use_ok 'Mashtree';
use Mashtree qw/treeDist/;

$ENV{PATH}="./bin:$ENV{PATH}";
$ENV{PATH}="./plugins:$ENV{PATH}";

my $db="mash.sqlite";
END{unlink($db);}

my $correctMashtree="(sample1:0.00210,sample2:0.00204,(sample3:0.00196,sample4:0.00207):0.00005);";
my $estimatedTree="(sample2:0.00307,sample3:0.00098,(sample1:0.00313,sample4:0.00108):0.00100);";

system("mashtree_init.pl $db > /dev/null 2>&1");
is $?, 0, "Ran mashtree_init.pl";

# Make mash distances but stagger them so that there
# is missing data.
system("mashtree_mash.pl --numcpus 1 $db t/lambda/sample*.fastq.gz > /dev/null 2>&1");
is $?, 0, "Ran mashtree_mash.pl";

# Delete one distance to test for missing data
subtest "Delete distance between 1 and 2" => sub{
  plan tests => 2;
  system(qq(sqlite3 $db "DELETE FROM DISTANCE WHERE GENOME1='sample2' AND GENOME2='sample1';"));
  is $?, 0, "Deleted distance from sample1 and sample2";
  system(qq(sqlite3 $db "DELETE FROM DISTANCE WHERE GENOME1='sample1' AND GENOME2='sample2';"));
  is $?, 0, "Deleted distance from sample1 and sample2";
};

#system("mashtree_dump.pl $db");
system("mashtree_guessBlanks.pl $db");
is $?, 0, "Ran mashtree_guessBlanks.pl to fill in distance between 1 and 2";

#system("mashtree_dump.pl $db");
my $tree = `mashtree_tree.pl $db `; 
chomp($tree);
is $tree, $estimatedTree, "mashtree_tree on estimated data";



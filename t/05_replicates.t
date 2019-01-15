#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use File::Basename qw/dirname/;
use Bio::TreeIO;
use IO::String;

use Test::More tests => 2;

use_ok 'Mashtree';
use Mashtree;

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="((sample2:0.0020443525,sample1:0.0021037373)66:0.0000540274,sample3:0.0019622177,sample4:0.0020673526)83;";
$correctMashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

# Test to see if the correct tree is made
END{unlink "lambdadist.tsv";}
my $mashtree=`mashtree_wrapper.pl --reps 5 --numcpus 2 $RealBin/lambda/*.fastq.gz`;
chomp($mashtree);
$mashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

my $fh = IO::String->new($mashtree);
my $tree = Bio::TreeIO->new(-fh=>$fh, -format=>"newick")->next_tree;
subtest "Bootstrapping test" => sub{
  plan tests => 3;
  my @nodes = $tree->get_nodes;
  for my $node(grep {!$_->is_Leaf} @nodes){
    ok($node->id > 1, "Bootstrap is greater than zero: ".$node->id);
  }

  my $correctNodeString = "sample1 sample2 sample3 sample4";
  my $nodeString = join(" ", sort map{$_->id} grep { $_->is_Leaf} @nodes);
  is $correctNodeString, $nodeString, "Taxon names in the tree";
};
  



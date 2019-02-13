#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use File::Basename qw/dirname/;
use Bio::TreeIO;
use IO::String;
use Scalar::Util qw/looks_like_number/;

use Test::More tests => 5;

use_ok 'Mashtree';
use Mashtree;
use Mashtree::Db;

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="((sample2:0.0020443525,sample1:0.0021037373)66:0.0000540274,sample3:0.0019622177,sample4:0.0020673526)83;";
$correctMashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

# Test to see if the correct tree is made
END{unlink "lambdadist.tsv"; system("rm -rf $RealBin/lambda/jackknife.tmp jackknife.log");}
my $mashtree=`mashtree_jackknife.pl --tempdir $RealBin/lambda/jackknife.tmp --reps 100 --numcpus 2 $RealBin/lambda/*.fastq.gz 2>jackknife.log`;
if($?){
  BAIL_OUT("Mashtree exited with error:\n".`cat jackknife.log`);
}
my $passed = ok(defined($mashtree),"Mashtree_jackknife.pl ran and produced a string");
$mashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

my $fh = IO::String->new($mashtree);
my $tree = Bio::TreeIO->new(-fh=>$fh, -format=>"newick")->next_tree;
$passed = is(ref($tree),"Bio::Tree::Tree","Produced a BioPerl tree object");
if(!$passed){
  BAIL_OUT("Tree was not produced out of this string:\n$mashtree");
}

subtest "Parts of the tree file intact" => sub{
  plan tests => 3;
  my @nodes = $tree->get_nodes;
  my @expectedBootstrap = (100, 11);
  my $nodeCounter=0;
  for my $node(grep {!$_->is_Leaf} @nodes){
    ok(looks_like_number($node->id), "Bootstrap is a number: ".$node->id);
    note("Usually this bootstrap is around $expectedBootstrap[$nodeCounter], give or take 5%");
    $nodeCounter++;
  }

  my $correctNodeString = "sample1 sample2 sample3 sample4";
  my $nodeString = join(" ", sort map{$_->id} grep { $_->is_Leaf} @nodes);
  is $correctNodeString, $nodeString, "Taxon names in the tree: $nodeString";
};
  
# Test to validate distances on the first rep
subtest "Database of distances for rep1" => sub{
  plan tests=>20;
  my $db = Mashtree::Db->new("$RealBin/lambda/jackknife.tmp/rep1/distances.sqlite");
  my @dist = split(/\n/, $db->toString('','phylip'));
  chomp(@dist);
  shift(@dist);
  for(my $i=0;$i<@dist;$i++){
    my ($name, @dist) = split(/\s+/, $dist[$i]);
    for(my $j=0; $j<@dist; $j++){
      ok(looks_like_number($dist[$j]), "Distance is a number")
        or note "  Found $dist[$j]";
      if($j==$i){
        ok($dist[$j] == 0,"Distance against self is 0");
      }
    }
  }
};

#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 2;

use_ok 'Mashtree::Mash';

my $expectedTree="((t/lambda/sample3.fastq.gz:0.17960,(t/lambda/sample1.fastq.gz:0.17876,(t/lambda/sample4.fastq.gz:0.17601,t/lambda/sample2.fastq.gz:0.17584):0.00205):0.00122):0.154475,(t/filetypes/CFSAN000968.ref.fasta.gz:0.00082,(t/filetypes/CFSAN001112.ref.fasta.gz:0.00031,(t/filetypes/CFSAN000191.ref.fasta.gz:0.00043,t/filetypes/CFSAN000189.ref.fasta.gz:0.00006):0.00018):0.00043):0.154475):0.01;";


# Since the tree can be somewhat random due to very similar
# profiles, test 100 trees.
# See if at least one is the correct one. I would like
# to stable-sort the tree somehow in the future to avoid
# this issue.
subtest 'Mashtree::Mash tree' => sub{
  my $numAttempts=500;
  for my $testCounter(0..$numAttempts-1){
    # Every so many tries, redo mash
    if($testCounter % 42 == 0){
      mashSketch($testCounter);
    }
    my $msh = Mashtree::Mash->new(
      ["t/lambda/lambda.msh","t/filetypes/filetypes.msh"]
    );
    my $tree = $msh->tree;
    my $newick = $tree->as_text("newick");
    if($newick=~/^\s*$/){
      fail("Empty tree");
      done_testing();
      return;
    }
    if($expectedTree eq $newick){
      pass("Found the correct tree on test number $testCounter");
      done_testing();
      return;
    }
  }
  fail("Did not find the correct tree after $numAttempts tries.");
};

sub mashSketch{
  my($seed)=@_;
  # Sketch a bunch of files to prepare.
  # TODO if I make a module for mash sketch, replace this line
  system("mash sketch -S $seed t/lambda/*.fastq.gz -o t/lambda/lambda.msh");
  die "ERROR running mash sketch on the lambda dataset" if $?;
  END{system("rm -f t/lambda/lambda.msh")}
  system("mash sketch -S $seed t/filetypes/*.fasta.gz -o t/filetypes/filetypes.msh");
  die "ERROR running mash sketch on the 'filetypes' dataset" if $?;
  END{system("rm -f t/filetypes/filetypes.msh")}
}


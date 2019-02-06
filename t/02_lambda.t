#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib './lib';
use File::Basename qw/dirname/;

use Test::More tests => 3;

use_ok 'Mashtree';
use Mashtree qw/treeDist/;

$ENV{PATH}="./bin:$ENV{PATH}";

my $correctMashtree="(sample3:0.00195,sample4:0.00205,(sample1:0.00205,sample2:0.00205):0.00010);";
$correctMashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression

# Test to see if the correct tree is made
END{unlink "lambdadist.tsv";}
my $mashtree=`mashtree --outmatrix lambdadist.tsv --numcpus 1 t/lambda/*.fastq.gz`;
chomp($mashtree);
$mashtree=~s/(\d+\.)(\d+)/$1 . substr($2,0,4)/ge; # global and expression
my $dist=treeDist($mashtree,$correctMashtree);
is $dist , 0, "Lambda test set tree";

# Test for the correct distance matrix
my %matrix=(
          'sample4' => {
                         'sample4' => 0,
                         'sample2' => '0.00417555',
                         'sample1' => '0.0042153'
                       },
          'sample2' => {
                         'sample2' => 0,
                         'sample1' => '0.00414809'
                       },
          'sample3' => {
                         'sample4' => '0.00402957',
                         'sample2' => '0.00405078',
                         'sample3' => 0,
                         'sample1' => '0.0041298'
                       },
          'sample1' => {
                         'sample1' => 0
                       }
        );
# mirror the matrix
while(my($ref,$queryHash)=each(%matrix)){
  while(my($query,$dist)=each(%$queryHash)){
    $matrix{$query}{$ref}=$dist;
  }
}

subtest "Test matrix" => sub {
  plan tests => 16;
  open(MATRIX, "lambdadist.tsv") or die "ERROR: could not read lambdadist.tsv: $!";
  my $header=<MATRIX>;
  chomp($header);
  my (undef,@header)=split(/\t/,$header);
  while(my $distances=<MATRIX>){
    chomp($distances);
    my($label,@dist)=split /\t/,$distances;
    for(my $i=0;$i<@header;$i++){
      # Round to 8 digits to be compliant with Mash
      my $got      = sprintf("%0.8f", $dist[$i]);
      my $expected = sprintf("%0.8f", $matrix{$label}{$header[$i]});
      is $got, $expected, "Distance between $label and $header[$i] (should be $dist[$i])";
    }
  }
  close MATRIX;
};




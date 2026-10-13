#!/usr/bin/env perl

# Test that --sort-order is honored when making the
# phylip file that goes into neighbor joining.
# https://github.com/lskatz/mashtree/issues/95

use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use lib "$RealBin/../lib/perl5";
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use List::Util qw/shuffle/;

use Test::More tests => 4;

use Mashtree;
use Mashtree::Db;

$ENV{PATH}="./bin:$ENV{PATH}";

my $tempdir = tempdir(basename($0).".XXXXXX", TMPDIR=>1, CLEANUP=>1);

# Return the sample names of a phylip string, in order
sub phylipNames{
  my($phylip)=@_;
  my @row = split(/\n/, $phylip);
  shift(@row); # number of samples
  return [map{(split(/\s+/, $_))[0]} @row];
}

# Return the sample names of a distance matrix string, in order
sub matrixNames{
  my($matrix)=@_;
  my($header) = split(/\n/, $matrix);
  my(undef, @name) = split(/\t/, $header);
  return \@name;
}

subtest 'Mashtree::Db sort order' => sub {
  plan tests => 12;

  # Names that are not in alphabetical order
  my @input = qw(zeta beta alpha);
  my @abc   = sort @input;

  my $db = Mashtree::Db->new("$tempdir/sortorder.db.tsv", {significant_figures=>10});
  my %dist;
  for my $i(0..$#input){
    for my $j(0..$#input){
      $dist{$input[$i]}{$input[$j]} = abs($i - $j) / 10;
    }
  }
  $db->addDistancesFromHash(\%dist);

  is_deeply(phylipNames($db->toString(\@input, "phylip")), \@abc, "phylip: default is alphabetical");
  is_deeply(phylipNames($db->toString(\@input, "phylip", "ABC")), \@abc, "phylip: ABC");
  is_deeply(phylipNames($db->toString(\@input, "phylip", "input-order")), \@input, "phylip: input-order");
  is_deeply(matrixNames($db->toString(\@input, "matrix", "input-order")), \@input, "matrix: input-order");
  is_deeply(phylipNames($db->toString(undef, "phylip", "input-order")), \@abc, "phylip: input-order without a genome list is alphabetical");

  # Distances must follow the names when they are reordered
  my %phylipDist;
  for my $row((split(/\n/, $db->toString(\@input, "phylip", "input-order")))[1..@input]){
    my($name, @d) = split(/\s+/, $row);
    $phylipDist{$name}{$input[$_]} = $d[$_] + 0 for(0..$#d);
  }
  my %matrixDist;
  for my $row((split(/\n/, $db->toString(\@input, "matrix", "input-order")))[1..@input]){
    my($name, @d) = split(/\t/, $row);
    $matrixDist{$name}{$input[$_]} = $d[$_] + 0 for(0..$#d);
  }
  is_deeply(\%phylipDist, \%dist, "phylip: distances follow input-order");
  is_deeply(\%matrixDist, \%dist, "matrix: distances follow input-order");

  # Names that truncate to the same name keep their exact order
  my $collisionDb = Mashtree::Db->new("$tempdir/collision.db.tsv", {significant_figures=>10});
  my @collision = qw(sample.fastq zeta sample.fasta alpha);
  my %collisionDist;
  for my $i(0..$#collision){
    for my $j(0..$#collision){
      $collisionDist{$collision[$i]}{$collision[$j]} = abs($i - $j) / 10;
    }
  }
  $collisionDb->addDistancesFromHash(\%collisionDist);
  is_deeply(phylipNames($collisionDb->toString(\@collision, "phylip", "input-order")), \@collision, "phylip: input-order with truncation collisions");

  # Random order with a controlled RNG
  srand(95);
  my @expected = shuffle(@input);
  srand(95);
  is_deeply(phylipNames($db->toString(\@input, "phylip", "random")), \@expected, "phylip: random with a fixed seed");
  srand(95);
  is_deeply(phylipNames($db->toString(\@input, "phylip", "rand")), \@expected, "phylip: rand with a fixed seed");

  # The random order should not always be alphabetical
  my $numShuffled = grep {"@{ phylipNames($db->toString(\@input, 'phylip', 'random')) }" ne "@abc"} 1..20;
  cmp_ok($numShuffled, '>', 0, "phylip: random is not always alphabetical ($numShuffled/20 shuffled)");

  eval{ $db->toString(\@input, "phylip", "nonsense") };
  like($@, qr/sort-order/, "phylip: unknown sort order dies");
};

# Run mashtree with sample names in nonalphabetical order
my @sample = qw(sample4 sample2 sample3 sample1);
my @reads  = map{"$RealBin/lambda/$_.fastq.gz"} @sample;
my %runs;
for my $sortOrder(qw(ABC input-order random)){
  my $wd = "$tempdir/$sortOrder";
  mkdir $wd;
  # Use more than one cpu so that sketching reorders the samples
  my $cmd = "mashtree --sort-order $sortOrder --tempdir $wd --outmatrix $wd/matrix.tsv --genomesize 40000 --numcpus 2 @reads > $wd/tree.dnd 2> $wd/mashtree.log";
  note $cmd;
  system($cmd);
  if($?){
    note `cat $wd/mashtree.log`;
    BAIL_OUT("ERROR running mashtree with --sort-order $sortOrder");
  }
  local $/ = undef;
  open(my $phylipFh, "<", "$wd/distances.phylip") or die "ERROR: could not read $wd/distances.phylip: $!";
  open(my $matrixFh, "<", "$wd/matrix.tsv") or die "ERROR: could not read $wd/matrix.tsv: $!";
  $runs{$sortOrder} = {
    phylip => phylipNames(<$phylipFh>),
    matrix => matrixNames(<$matrixFh>),
  };
}

subtest 'mashtree --sort-order ABC' => sub {
  plan tests => 2;
  is_deeply($runs{ABC}{phylip}, [sort @sample], "phylip is alphabetical");
  is_deeply($runs{ABC}{matrix}, [sort @sample], "matrix is alphabetical");
};

subtest 'mashtree --sort-order input-order' => sub {
  plan tests => 2;
  is_deeply($runs{'input-order'}{phylip}, \@sample, "phylip is in input order");
  is_deeply($runs{'input-order'}{matrix}, \@sample, "matrix is in input order");
};

subtest 'mashtree --sort-order random' => sub {
  plan tests => 2;
  is_deeply([sort @{ $runs{random}{phylip} }], [sort @sample], "phylip has all samples");
  is_deeply($runs{random}{matrix}, $runs{random}{phylip}, "matrix is in the same order as the phylip file");
};

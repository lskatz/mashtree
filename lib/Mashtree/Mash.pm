#!/usr/bin/env perl
package Mashtree::Mash;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;

use lib dirname($INC{"Mashtree/Mash.pm"});
use lib dirname($INC{"Mashtree/Mash.pm"})."/..";

use Mashtree qw/logmsg/;
use JSON qw/from_json/;
use Bio::Seq;
use Bio::Tree::Tree;
use Bio::Matrix::Generic;
use Bio::SimpleAlign;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;

our @EXPORT_OK = qw(
         );

local $0=basename $0;

# If this is used in a scalar context, $self->toString() is called
use overload '""' => 'toString';

# Properties of this object:
sub new{
  my($class,$file,$settings)=@_;

  if(ref($file) ne 'ARRAY'){
    die "ERROR: the first parameter must be a list of mash file(s)";
  }

  my $self={
    file      => $file,
    info      => {},
    hashes    => {},
    distance  => {},
    cluster   => [],
    aln       => [],
    refinedTree=>"",  # set to bool false but will be set to Bio::Tree::Tree
  };
  bless($self,$class);

  # Gather info from each file. $self->{info} and
  # $self->{hashes} gets updated.
  for my $f(@$file){
    $self->info($f);
  }

  # TODO are all mash metadata compatible?  Else: ERROR

  return $self;
}

sub info{
  my($self,$msh)=@_;
  
  my %info = %{ $self->{info} };

  my $mashInfo=from_json(`mash info -d $msh`);

  for my $sketch(@{ $$mashInfo{sketches} }){
    #delete($$sketch{hashes});
    $info{$$sketch{name}}=$sketch;

    my %sketchHash;
    for my $pos(@{ $$sketch{hashes} }){
      $$self{hashes}{$pos}++;
      $sketchHash{$pos}=1;
    }
    $info{$$sketch{name}}{hashes}=\%sketchHash;
  }

  $self->{info}=\%info;

  return \%info;
}

sub mashDistances{
  my($self)=@_;
  # Don't recalculate
  return $self->{distance} if(keys(%{ $self->{distance} }) > 0);

  my @file=@{$self->{file}};
  my $numFiles=@file;
  
  my %distance;
  for(my $i=0;$i<$numFiles;$i++){
    for(my $j=0;$j<$numFiles;$j++){
      open(my $fh, "mash dist '$file[$i]' '$file[$j]' | ") or die "ERROR: running mash dist on $file[$i] and $file[$j]: $!";
      while(<$fh>){
        chomp;
        my($genome1,$genome2,$dist,$p,$sharedFraction)=split(/\t/,$_);
        $distance{$genome1}{$genome2}=$dist;
      }
      close $fh;
    }
  }
  $self->{distance}=\%distance;
  return \%distance;
}

# Returns clusters
sub clusters{
  my($self,$cutoff)=@_;
  $cutoff||=0.10;

  if(scalar(@{ $self->{cluster} }) > 0){
    return $self->{cluster};
  }

  my $distance=$self->mashDistances;

  my @genome=sort {$a cmp $b} keys(%$distance);
  my $numGenomes=@genome;
  # Initialize the set of clusters with the first genome.
  # The first element of each cluster is the 'seed'
  my @cluster=();
  $cluster[0]=[$genome[0]];
  for(my $i=1;$i<$numGenomes;$i++){
    # Compare each genome against each seed
    my $was_clustered=0;
    for(my $j=0;$j<@cluster;$j++){
      my $seed=$cluster[$j][0];
      # See if the genome is close to any clusters.
      if($$distance{$seed}{$genome[$i]} < $cutoff){
        push(@{ $cluster[$j] }, $genome[$i]);
        $was_clustered=1;
        # avoid having multiple genomes going to multiple clusters
        last;
      }
    }
    # If the genome didn't get clustered with anything,
    # make a new cluster.
    if(!$was_clustered){
      push(@cluster, [$genome[$i]]);
    }
  }
  $self->{cluster}=\@cluster;
  return \@cluster;
}

# Make sets of alignments
sub alignments{
  my($self)=@_;

  if(scalar(@{ $self->{aln} }) > 0){
    return $self->{aln};
  }

  # print an A for present, a T for not present
  my %nt=(present=>"A",absent=>"T");
  
  # Get clusters
  # TODO make dist cutoff value accessible here via $self->new()
  my $cluster=$self->clusters;

  # loop through each cluster
  # return set of alns
  my @aln;
  for my $genomeArr(@$cluster){
    my $aln=Bio::SimpleAlign->new();
    for my $name(@$genomeArr){
      my $sketch=$$self{info}{$name};
      my $sequence="";
      for my $hash(sort {$a<=>$b} keys(%{ $self->{hashes} })){
        if($$sketch{hashes}{$hash}){
          $sequence.= $nt{present};
        } else {
          $sequence.= $nt{absent};
        }
      } 
      my $seq=Bio::LocatableSeq->new(-id=>$name, -seq=>$sequence);
      $aln->add_seq($seq);
    }
    push(@aln, $aln);
  }
  
  $self->{aln}=\@aln;
  return \@aln;
}

# Make a tree for each cluster, then glue them together 
# using mash distances for branch length.
sub refinedTree{
  my($self)=@_;

  if($self->{refinedTree}){
    return $self->{refinedTree};
  }

  my $cluster=$self->clusters;
  my @seed = map{$$_[0]} @$cluster;
  my $alnArr=$self->alignments;
  my $distance=$self->mashDistances;

  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $matrix   = Bio::Matrix::Generic->new(-rownames=>\@seed,-colnames=>\@seed);

  die "TODO figure out to do with clusters < 3" if(@$cluster < 3);

  # Set the distances in the distance matrix
  for(my $i=0;$i<@seed;$i++){
    my $genome1=$seed[$i];
    for(my $j=0;$j<@seed;$j++){
      my $genome2=$seed[$j];
      $matrix->entry($genome1,$genome2,$$distance{$genome1}{$genome2});
    }
  }
  my $refinedTree = $dfactory->make_tree($matrix);

  my $node_with_longest_branch = _reroot_at_midpoint($refinedTree);
  my $refinedTree_root = $refinedTree->reroot_at_midpoint($node_with_longest_branch,"root_refined");
  $refinedTree_root->branch_length(0.01);

  # Make the individual trees
  for(my $i=0;$i<@$alnArr;$i++){
    # Don't run a new tree on a cluster of one.
    # Skip because it is already represented as a seed.
    if(@{ $$cluster[$i] } < 2){
      next;
    }

    my $aln=$$alnArr[$i];
    my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
    my $stats    = Bio::Align::DNAStatistics->new;
    my $matrix   = $stats->distance(-method=>'Uncorrected', -align => $aln);
    my $tree     = $dfactory->make_tree($matrix);

    # Root on the seed in the subtree.
    # Next, add the seed node from the subtree to the 
    # main tree.
    my $subtree_root = _reroot_at_midpoint($tree);
    my $refinedTree_seedNode=$refinedTree->find_node(-id=>$seed[$i]);
    my $ancestorNode=$refinedTree_seedNode->ancestor;
    for my $node($subtree_root->each_Descendent()){
      if($node->branch_length==0){
        $node->branch_length(0.001);
      }
      $ancestorNode->add_Descendent($node);
    }
    # Remove the seed node itself: it is duplicated in the
    # subtree.
    $ancestorNode->remove_Descendent($refinedTree_seedNode);
  }
  $refinedTree->contract_linear_paths;

  $self->{refinedTree}=$refinedTree;

  return $refinedTree;
}

sub _reroot_at_midpoint{
  my($tree)=@_;
  my $node_with_longest_branch = (sort{
    my $A=$a->branch_length || 0;
    my $B=$b->branch_length || 0;
    $B <=> $A
  } $tree->get_nodes)[0];
  my $rootNode=$tree->reroot_at_midpoint($node_with_longest_branch);
  if(!$rootNode->branch_length || $rootNode->branch_length < 0.01){
    $rootNode->branch_length(0.01);
  }
  return $rootNode;
}

sub toString{
  my($self)=@_;
  my $fileArr=$self->{file};
  my $return="Mashtree::Mash object with " .scalar(@$fileArr)." files:\n\n";
  
  for(@$fileArr){
    $return.="-------${_}------\n";
    my $info=`mash info '$_'`;
    chomp($info);
    $return.=$info;
    $return.="^^^^^^^${_}^^^^^^\n\n";
  }
  
  return $return;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


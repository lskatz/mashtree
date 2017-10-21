#!/usr/bin/env perl
package Mashtree::Mash;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;

use lib dirname($INC{"Mashtree/Mash.pm"});
use lib dirname($INC{"Mashtree/Mash.pm"})."/..";

use Mashtree qw/_truncateFilename logmsg/;
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


# Properties of this object:
sub new{
  my($class,$file,$settings)=@_;

  if(ref($file) ne 'ARRAY'){
    die "ERROR: parameter \$file must be a list of file(s)";
  }

  my $self={
    file      => $file,
    info      => {},
    hashes    => {},
    distance  => {},
    cluster   => [],
    aln       => [],
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
  $cutoff||=0.05;

  if(scalar(@{ $self->{cluster} }) > 0){
    return $self->{cluster};
  }

  my $distance=$self->mashDistances;

  my @genome=keys(%$distance);
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

  my $cluster=$self->clusters;
  my @seed = map{$$_[0]} @$cluster;
  my $alnArr=$self->alignments;
  my $distance=$self->mashDistances;

  my $comprehensiveTree=Bio::Tree::Tree->new;
  # Make a tree with all the seeds if there are at least
  # three clusters.
  # TODO if there are only one or two clusters, make a blank tree with a root node.
  if(@$cluster < 3){
    ...;
  } elsif(@$cluster >= 3){
    my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
    my $matrix   = Bio::Matrix::Generic->new(-rownames=>\@seed,-colnames=>\@seed);

    # Set the distances in the distance matrix
    for(my $i=0;$i<@seed;$i++){
      my $genome1=$seed[$i];
      for(my $j=0;$j<@seed;$j++){
        my $genome2=$seed[$j];
        $matrix->entry($genome1,$genome2,$$distance{$genome1}{$genome2});
      }
    }
        
    $comprehensiveTree = $dfactory->make_tree($matrix);
  }

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
    my $rootNode = $tree->get_root_node;

    # TODO HOW TO MERGE TREES????
    foreach my $node($rootNode->get_all_Descendents()){
      next if(!$node->is_Leaf);
      $comprehensiveTree->merge_lineage($node);
      last;
    }
    die $comprehensiveTree->as_text;
      
    my $seedNode = $tree->findnode_by_id($seed[$i]);
    die Dumper $seedNode;
    print $tree->as_text;
    $comprehensiveTree->merge_lineage($seedNode);
    print $tree->as_text()."\n";
  }
  die $comprehensiveTree->as_text;

  # glue together the trees

  return $comprehensiveTree;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


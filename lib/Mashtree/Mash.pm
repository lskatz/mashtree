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

=pod

=head1 NAME

Mashtree::Mash

=head1 SYNOPSIS

A module to read `mash info` output and transform it

  use strict;
  use warnings;
  use Mashtree::Mash

  # Quick example

  # Sketch all fastq files into one mash file.
  system("mash sketch *.fastq.gz > all.msh");
  die if $?;
  # Read the mash file.
  my $msh = Mashtree::Mash->new(["all.msh"]);
  # Get a Bio::Tree::Tree object
  my $tree = $msh->refinedTree;

  # Do something with the tree, e.g.,
  print $tree->as_text("newick");

=head1 DESCRIPTION

This is a module to read mash files produced by the Mash executable. For more information on Mash, see L<mash.readthedocs.org>.  This module is capable of reading mash files, estimating which entries are clustered, producing a pseudo-multiple sequence alignment (MSA), and creating a tree based on either mash distances or the pseudo-MSA.  The MSA is generated with A being present and T being absent.

=head1 METHODS

=over

=item Mashtree::Mash->new(\@listOfFiles,\%options);

Create a new instance of Mashtree::Mash.  One object per set of files.

  Applicable arguments:
  None so far (TODO)

=back

=cut

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

  # Test that all metadata for mash are the same as the
  # first genome in the set.  I.e., they are all the same.
  my $info=$self->{info};
  my @genomeName=keys(%$info);
  my $ref=shift(@genomeName); 
  for my $g(@genomeName){
    # TODO don't worry about sketch size b/c you can take
    # first X sketches from each where X is the smallest
    # number of sketches in the set. For now though to make
    # things simple, just take exact numbers.
    for my $key(qw(kmer alphabet preserveCase canonical sketchSize hashType hashBits hashSeed)){
      if($$info{$g}{$key} ne $$info{$ref}{$key}){
        die "ERROR: genomes $ref and $g are incompatible under property $key";
      }
    }
  }

  return $self;
}

=pod

=over

=item $msh->info()

Returns a hash ref that describes a single mash file. Updates the Mashtree::Mash object with this info. This method is ordinarily not used externally to the object.

  Arguments: One mash file
  Returns:   Reference to a hash

=back

=cut

sub info{
  my($self,$msh)=@_;
  
  my %info = %{ $self->{info} };

  if(!   $msh){
    logmsg "WARNING: no file was given to \$self->info.  Did you mean to call the subroutine or the hash element 'info'?";
    return {};
  }
  if(!-e $msh){
    die "ERROR: could not find file $msh";
  }

  my $mashInfo=from_json(`mash info -d $msh`);

  for my $sketch(@{ $$mashInfo{sketches} }){
    #delete($$sketch{hashes}); logmsg "DEBUG: removing hashes element";
    $info{$$sketch{name}}=$sketch;

    my %sketchHash;
    for my $pos(@{ $$sketch{hashes} }){
      $$self{hashes}{$pos}++;
      $sketchHash{$pos}=1;
    }
    $info{$$sketch{name}}{hashes}=\%sketchHash;

    # Also take on the general properties of the mash file
    for my $key(qw(kmer alphabet preserveCase canonical sketchSize hashType hashBits hashSeed)){
      $info{$$sketch{name}}{$key}=$$mashInfo{$key};
    }
  }

  $self->{info}=\%info;

  return \%info;
}

=pod

=over

=item $msh->mashDistances()

Returns a reference to a hash of mash distances

=back

=cut

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

=pod

=over

=item $msh->clusters

Estimates which entries belong to which cluster.  The default cutoff is a mash distance of 0.10, but the user can supply a different cutoff.  Supplying a cutoff of 0 indicates that each new mash profile is a separate cluster; 1 means that everything belongs to the same cluster.

This affects downstream methods that rely on knowing which genomes should be compared with each other in a more refined way.

  Arguments: cutoff (decimal)
  Returns:   List of lists.  Each sub-list is a cluster.

=back

=cut

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

=pod

=over

=item $msh->alignments()

Return a list of Bio::Align objects, one per cluster. If $msh->cluster is not called beforehand, it will be called internally with default values.

  Arguments: none
  Returns:   List of Bio::Align objects

=back

=cut

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

=pod

=over

=item $msh->refinedTree

Make a tree for each cluster, then glue them together using mash distances for branch length.  $msh->alignments is internally called if it has not already been called.

  Arguments: None
  Returns:   Bio::Tree::Tree object

=back

=cut

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
  #$refinedTree_root->branch_length(0.01);

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
        #$node->branch_length(0.001);
      }
      $ancestorNode->add_Descendent($node);
    }
    # Remove the seed node itself: it is duplicated in the
    # subtree.
    $ancestorNode->remove_Descendent($refinedTree_seedNode);
  }
  #$refinedTree->contract_linear_paths;

  $self->{refinedTree}=$refinedTree;

  return $refinedTree;
}

# TODO subroutine to make just the mash distances tree
# TODO subroutine to make individual trees based on kmer presence/absence


##### Utility methods

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

=pod

=head1 COPYRIGHT AND LICENSE

MIT license.

=head1 AUTHOR

Author:  Lee Katz <lkatz@cdc.gov>

For additional help, go to https://github.com/lskatz/Mashtree

CPAN module at http://search.cpan.org/~lskatz/Mashtree

=for html <a href="https://travis-ci.org/lskatz/mashtree"><img src="https://travis-ci.org/lskatz/mashtree.svg?branch=master"></a>

=cut

1; # gotta love how we we return 1 in modules. TRUTH!!!


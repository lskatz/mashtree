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
  my $tree = $msh->tree;

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
    names     => [],
    hashes    => {},
    distance  => {},
    aln       => [],
    tree=>"",  # set to bool false but will be set to Bio::Tree::Tree
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

  # Set a sorted list of names
  $self->{names}=[sort {$a cmp $b} keys(%{ $self->{info} })];

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

=item $msh->alignments()

Return a Bio::SimpleAlign object based on presence/absence of kmer hashes.

  Arguments: none
  Returns:   Bio::SimpleAlign

=back

=cut

sub alignment{
  my($self)=@_;

  if(scalar(@{ $self->{aln} }) > 0){
    return $self->{aln};
  }

  # print an A for present, a T for not present
  my %nt=(present=>"A",absent=>"T");

  my $aln=Bio::SimpleAlign->new();
  for my $name(keys(%{ $self->{info} })){
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
  
  $self->{aln}=$aln;
  return $aln;
}

=pod

=over

=item $msh->tree

Make a tree based on $msh->alignment. $msh->alignment is internally called if it has not already been called.

  Arguments: None
  Returns:   Bio::Tree::Tree object

=back

=cut

sub tree{
  my($self)=@_;

  if($self->{tree}){
    return $self->{tree};
  }

  my $aln=$self->alignment;

  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $stats    = Bio::Align::DNAStatistics->new;
  my $matrix   = $stats->distance(-method=>'Uncorrected', -align => $aln);
  my $tree     = $dfactory->make_tree($matrix);

  my $subtree_root = _reroot_at_midpoint($tree);

  $self->{tree}=$tree;

  return $tree;
}

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


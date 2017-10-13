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
    file=>$file,
    info=>{},
    hashes=>{},
    distance=>{},
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

# TODO make sets of alignments and also report the 
# mash distance between them all.
sub alignment{
  my($self)=@_;

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
  return $aln;
}

sub refinedTree{
  my($self,$aln)=@_;
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $stats    = Bio::Align::DNAStatistics->new;
  my $matrix   = $stats->distance(-method=>'Kimura', -align => $aln);
  my $tree     = $dfactory->make_tree($matrix);
  #print $tree->as_text()."\n";
  return $tree;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


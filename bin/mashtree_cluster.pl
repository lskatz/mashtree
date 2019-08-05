#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Read a mashtree database to cluster genomes

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename dirname fileparse/;
use File::Copy qw/copy/;
use POSIX qw/floor/;
use List::Util qw/min max/;
use Scalar::Util qw/looks_like_number/;

use DBI;

use threads;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg/;
#use Mashtree qw/logmsg @fastqExt @fastaExt @mshExt @richseqExt _truncateFilename createTreeFromPhylip $MASHTREE_VERSION/;
#use Mashtree::Db;
#use Graph::Dijkstra;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help threshold|cutoff=f numcpus=i nonzero=f)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{nonzero}||=1e-99;
  $$settings{threshold}||=0.1;

  my($dbFile) = @ARGV;

  die usage() if($$settings{help} || !$dbFile);

  if(!-e $dbFile){
    die "ERROR: could not find database file $dbFile";
  }

  my $clusters = makeClusters($dbFile,$settings);
  for my $clusterArr(values(%$clusters)){
    print join("\t",@$clusterArr)."\n";
  }

  return 0;
}

sub makeClusters{
  my($db,$settings)=@_;

  my %C; #clustering hash, where the key is the seed

  #my $graph = Graph::Dijkstra->new(edgedefault=>"undirected");
  my %graph=();# hash of pairwise nodes

  my $dbh = DBI->connect("dbi:SQLite:dbname=$db","","",{
      RaiseError => 1
  });

  my $sth = $dbh->prepare(qq(SELECT GENOME1,GENOME2,DISTANCE
    FROM DISTANCE
    WHERE GENOME1 != GENOME2
    ORDER BY DISTANCE ASC, 
             GENOME1  ASC, 
             GENOME2  ASC
    ));
  my $rv = $sth->execute() or die $DBI::errstr;
  if($rv < 0){
    die "ERROR: no distances were found in the database $db";
  }

  my $rowCounter=0;
  my %genomeIndex;
  while(my @row=$sth->fetchrow_array()){
    my($g1,$g2,$dist)=@row;
    $genomeIndex{$g1}=1;
    $genomeIndex{$g2}=1;

    #$dist=$$settings{nonzero} if(!$dist);
    my($G1,$G2)=sort{$a cmp $b} ($g1,$g2);
    $graph{$G1}{$G2}=$dist;
    $graph{$G2}{$G1}=$dist;
    $rowCounter++;

    if($rowCounter % 1000000 == 0){
      logmsg "Read the distances for $rowCounter pairs";
    }
  }
  logmsg "Done reading the distances for $rowCounter pairs";
  # TODO? close the database fh

  logmsg "Calculating clusters";
  #my @genome=sort {$a cmp $b} map {$$_{id}} $graph->nodeList;
  my @genome = sort{$a cmp $b} keys(%genomeIndex);
  my $numGenomes=@genome;
  for(my $i=0;$i<$numGenomes;$i++){
    # Is this genome $i close to anything?
    my $query = $genome[$i];
    logmsg "Querying with $query";
    my $hit = undef;

    # See if the query is close to any of the current seeds first
    my @seed           = sort{$a cmp $b} keys(%C);
    my $numSeeds       = scalar(@seed);
    my $possibleHit;
    my $closestSeedDist = ~0;
    for(my $j=0; $j<$numSeeds; $j++){
      next if(!defined($graph{$seed[$j]}{$query}));
      my $dist = $graph{$seed[$j]}{$query};
      next if($dist > $$settings{threshold});

      if($dist < $closestSeedDist){
        $possibleHit = $seed[$j];
        $closestSeedDist = $dist;
      }
    }

    # If it's close to a seed, then record it.
    if($closestSeedDist <= $$settings{threshold}){
      logmsg "  => $possibleHit";
      push(@{ $C{$possibleHit} }, $query);
    }
    # If it's not close to a seed, start a new seeded cluster.
    else {
      logmsg "  => new cluster seeded with $query";
      $C{$query} = [$query];
    }
    
  }
  return \%C;
}

sub usage{
  "Output clusters of genomes, based on a mashtree sqlite database.
  The stdout will be a tab-delimited list of genomes, one
  line per cluster.
  Usage: $0 file.sql > clusters.tsv

  --threshold  0.1       Maximum any two genomes can be
                         from the seed in a given cluster
  --nonzero    1e-99     Zero distance is not tolerated
                         in this script. Give a nonzero
                         value in case a zero distance
                         is found.
  --numcpus    1         Not currently used.
  "
}


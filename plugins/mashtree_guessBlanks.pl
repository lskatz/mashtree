#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Read a mashtree database to cluster genomes

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;
use File::Basename qw/basename dirname fileparse/;
use File::Copy qw/mv/;
use POSIX qw/floor/;
use List::Util qw/min max/;
use Scalar::Util qw/looks_like_number/;

use DBI;

use threads;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg/;
use Graph::Dijkstra;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s method=s threshold|cutoff=f numcpus=i nonzero=f)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{nonzero}||=1e-99;
  $$settings{threshold}||=0.1;

  my($dbFile) = @ARGV;
  die usage() if($$settings{help} || !$dbFile);

  if(!$$settings{method}){
    $$settings{method}="dijkstra";
    logmsg "Setting method to dijkstra.";
  }
  $$settings{method}=lc($$settings{method});
  $$settings{tempdir}||=tempdir("MASHTREE_GUESS_$$settings{method}.XXXXXX",CLEANUP=>1,TMPDIR=>1);

  if(!-e $dbFile){
    die "ERROR: could not find database file $dbFile";
  }

  optimizeDb($dbFile,$$settings{method},$settings);
  saveOptimizations($dbFile,$settings);

  return 0;
}

# Optimize the database using the method supplied. This
# subroutine will redirect the workload to whichever
# specialized subroutine.
sub optimizeDb{
  my($dbFile,$method,$settings)=@_;

  no strict 'refs';
  my $subroutine="optimizeDb_$method";
  
  if(!defined(&$subroutine)){
    die "ERROR: optimization method $method is not defined in this script";
  }

  return &$subroutine($dbFile,$settings);
}

sub optimizeDb_dijkstra{
  my($db,$settings)=@_;

  my $graph = Graph::Dijkstra->new(edgedefault=>"undirected");

  my $dbh = DBI->connect("dbi:SQLite:dbname=$db","","",{
      RaiseError => 1,
      AutoCommit => 0,
  });

  logmsg "Reading the database";
  my $sth = $dbh->prepare(qq(
    SELECT GENOME1,GENOME2,DISTANCE
    FROM DISTANCE
    WHERE GENOME1 != GENOME2
  ));
  my $rv = $sth->execute() or die $DBI::errstr;
  if($rv < 0){
    die "ERROR: no distances were found in the database $db";
  }

  my $rowCounter=0;
  my %seenNode; # a fast hash to track whether we need to define a node
  while(my @row=$sth->fetchrow_array()){
    my($g1,$g2,$dist)=@row;
    # next if($g1 eq $g2); # take care of this with SQL
    $dist=$$settings{nonzero} if(!$dist);

    if(!$seenNode{$g1}){
      $graph->node({id=>$g1});
      $seenNode{$g1}=1;
    }
    if(!$seenNode{$g2}){
      $graph->node({id=>$g2});
      $seenNode{$g2}=1;
    }
    $graph->edge({sourceID=>$g1,targetID=>$g2,weight=>$dist});
    $rowCounter++;

    if($rowCounter % 1000000 == 0){
      logmsg "$rowCounter distances read";
    }
  }
  undef(%seenNode);

  # Create the new table
  $sth = $dbh->prepare(qq(
    DROP TABLE IF EXISTS OPTIMIZED_DISTANCE
    ));
  $sth->execute();
  $sth = $dbh->prepare(qq(
    CREATE TABLE OPTIMIZED_DISTANCE(
      GENOME1     CHAR(255)    NOT NULL,
      GENOME2     CHAR(255)    NOT NULL,
      DISTANCE    REAL         NOT NULL,
      PRIMARY KEY(GENOME1,GENOME2)
    )) 
  );
  $sth->execute();

  my $insertSth = $dbh->prepare("INSERT INTO OPTIMIZED_DISTANCE(GENOME1,GENOME2,DISTANCE)
    VALUES(?,?,?);");
  $dbh->{AutoCommit}=0;

  my %distBuffer;       # buffer for db transactions
  my $bufferSize=10000; # how many transactions at a time for the db
  my @genome=sort {$a cmp $b} map {$$_{id}} $graph->nodeList;
  my $numGenomes=@genome;
  my $distancesCounter=0;
  logmsg "Optimizing $numGenomes genomes";
  for(my $i=0;$i<$numGenomes;$i++){
    logmsg "Optimizing $genome[$i] (".($i+1)."/$numGenomes)";
    for(my $j=$i;$j<$numGenomes;$j++){
      # Sort the genome names so that we only store the pair once
      # in a predictable fashion.
      my($G1,$G2)=sort{$a cmp $b} ($genome[$i],$genome[$j]);
      my %solution=(originID=>$G1, destinationID=>$G2);
      if($G1 eq $G2){
        $solution{weight}=0;
      } else {
        $graph->shortestPath(\%solution);
      }
      #if(!$solution{edges}){
      #  if($G1 eq $G2){
      #    $solution{weight} = 0;
      #  } else {
      #    logmsg "WARNING: could not find distance between $G1 and $G2. Setting to large INT.";
      #    $solution{weight}= ~0; # largest int
      #  }
      #  $solution{edges}=[];
      #}

      $distBuffer{$G1}{$G2}=$solution{weight};
      $distancesCounter++;

      # Insert the distances if the buffer is full
      if($distancesCounter % $bufferSize == 0){
        my $insertValues="";
        while(my($G1,$distances)=each(%distBuffer)){
          while(my($G2,$dist) = each(%$distances)){
            $insertSth->execute($G1,$G2,$dist);
          }
        }
        $dbh->commit();
        %distBuffer=(); # purge the buffer
      }
    }
  }

  # One last set of inserts
  while(my($G1,$distances)=each(%distBuffer)){
    while(my($G2,$dist) = each(%$distances)){
      $insertSth->execute($G1,$G2,$dist);
    }
  }
  $dbh->commit();
  #system("sqlite3 $db .dump");
}

sub saveOptimizations{
  my($dbFile, $settings)=@_;

  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbFile","","",{
      RaiseError => 1,
      AutoCommit => 0,
  });

  my $dropSth = $dbh->prepare("
    ALTER TABLE DISTANCE
    RENAME TO DISTANCE_BACKUP
  ");
  $dropSth->execute() or die "ERROR dropping table DISTANCE: ".$DBI::errstr;

  my $alterSth = $dbh->prepare("
    ALTER TABLE OPTIMIZED_DISTANCE
    RENAME TO DISTANCE
  ");
  $alterSth->execute() or die "ERROR renaming table OPTIMIZED_DISTANCE to DISTANCE: ".$DBI::errstr;

  $dbh->commit();
  #system("sqlite3 $dbFile .dump");

  return 1;  
}

sub usage{
  "Estimates distances for unknown values in a mashtree database using different methods
  Usage: $0 [options] mashtree.sqlite

  Note: only one optimization is available at this time

  --method   Dijkstra   Estimate distances using Dijkstra's algorithm
  "
}


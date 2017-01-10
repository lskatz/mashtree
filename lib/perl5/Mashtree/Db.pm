#!/usr/bin/env perl
package Mashtree::Db;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use DBI;

our @EXPORT_OK = qw(
         );

local $0=basename $0;

# Properties of this object:
#   dbFile
#   dbh
sub new{
  my($class,$dbFile,$settings)=@_;

  my $self={};
  bless($self,$class);

  $self->selectDb($dbFile);
  return $self;
}

# Create an SQLite database for genome distances.
sub selectDb{
  my($self, $dbFile)=@_;

  $self->{dbFile}=$dbFile;

  $self->connect();

  if(-e $dbFile && -s $dbFile > 0){
    return 0;
  }

  my $dbh=$self->{dbh};
  $dbh->do(qq(
    CREATE TABLE DISTANCE(
      GENOME1     CHAR(255)    NOT NULL,
      GENOME2     CHAR(255)    NOT NULL,
      DISTANCE    INT          NOT NULL,
      PRIMARY KEY(GENOME1,GENOME2)
    )) 
  );

  return 1;
}

sub connect{
  my($self)=@_;

  my $dbFile=$self->{dbFile};
  my $dbh=DBI->connect("dbi:SQLite:dbname=$dbFile","","",{
      RaiseError => 1
  });
  
  $self->{dbh}=$dbh;
  
  return $dbh;
}

sub addDistances{
  my($self,$distancesFile)=@_;

  my $dbh=$self->{dbh};

  open(my $fh, "<", $distancesFile) or die "ERROR: could not read $distancesFile: $!";
  my $query="";
  while(<$fh>){
    chomp;
    if(/^#\s*query\s+(.+)/){
      $query=$1;
      next;
    }
    my($subject,$distance)=split(/\t/,$_);
    
    next if(defined($self->findDistance($query,$subject)));

    $dbh->do(qq(
      INSERT INTO DISTANCE VALUES("$query", "$subject", $distance);
    ));
    
    if($dbh->err()){
      die "ERROR: could not insert $query,$subject,$distance into the database:\n  ".$dbh->err();
    }
  }
}

sub findDistance{
  my($self,$genome1,$genome2)=@_;

  my $dbh=$self->{dbh};
  
  my $sth=$dbh->prepare(qq(SELECT DISTANCE FROM DISTANCE WHERE GENOME1="$genome1" AND GENOME2="$genome2"));
  my $rv = $sth->execute() or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  # Distance will be undefined unless there is a result
  # on the SQL select statement.
  my $distance;
  while(my @row=$sth->fetchrow_array()){
    ($distance)=@row;
  }
  return $distance;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


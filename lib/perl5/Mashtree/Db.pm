#!/usr/bin/env perl
package Mashtree::Db;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use List::Util qw/shuffle/;
use Data::Dumper;
use DBI;

use lib dirname($INC{"Mashtree/Db.pm"});
use lib dirname($INC{"Mashtree/Db.pm"})."/..";

use Mashtree qw/_truncateFilename logmsg sortNames/;


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

sub clone{
  my($self)=@_;

  my $clone={};
  $clone={
    dbFile    => $self->{dbFile},
    dbh       => $self->{dbh}->clone,
  };
  bless($clone,"Mashtree::Db");
  return $clone;
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

# No need to call this usually because the object will 
# disconnect from the database when it is destroyed.
sub disconnect{
  my($self)=@_;
  $self->{dbh}->disconnect();
}

sub addDistances{
  my($self,$distancesFile)=@_;

  my $dbh=$self->{dbh};
  my $numInserted=0;   # how many are going to be inserted?

  my $insertSQL="INSERT INTO DISTANCE VALUES";
  open(my $fh, "<", $distancesFile) or die "ERROR: could not read $distancesFile: $!";
  my $query="";
  while(<$fh>){
    chomp;
    if(/^#\s*query\s+(.+)/){
      $query=$1;
      $query=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
      $query=_truncateFilename($query);
      next;
    }
    my($subject,$distance)=split(/\t/,$_);
    $subject=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
    $subject=_truncateFilename($subject);
    
    next if(defined($self->findDistance($query,$subject)));

    $numInserted++;
    $insertSQL.=qq( ("$query", "$subject", $distance), );

  }

  $insertSQL=~s/\s*,\s*$//; # remove whitespace and comma from the end

  if($numInserted == 0){
    return $numInserted;
  }
  
  $dbh->do($insertSQL);
  if($dbh->err()){
    die "ERROR: could not insert $distancesFile into the database with query\n  $insertSQL\n  ".$dbh->err();
  }

  return $numInserted;
}

# Find all distances with a single genome
sub findDistances{
  my($self,$genome1)=@_;

  my $dbh=$self->{dbh};
  
  my $sth=$dbh->prepare(qq(SELECT GENOME2,DISTANCE 
    FROM DISTANCE 
    WHERE GENOME1="$genome1"
    ORDER BY GENOME2
  ));
  my $rv = $sth->execute() or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  # Distance will be undefined unless there is a result
  # on the SQL select statement.
  my %distance;
  while(my @row=$sth->fetchrow_array()){
    $distance{$row[0]}=$row[1];
  }
  return \%distance;
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

# Format can be:
#   tsv    3-column format
#   phylip Phylip matrix format
# sortBy can be:
#   abc
#   rand
sub toString{
  my($self,$format,$sortBy)=@_;
  $format//="tsv";
  $format=lc($format);
  $sortBy//="abc";
  $sortBy=lc($sortBy);
  
  if($format eq "tsv"){
    return $self->toString_tsv($sortBy);
  } elsif($format eq "phylip"){
    return $self->toString_phylip($sortBy);
  }

  die "ERROR: could not format ".ref($self)." as $format.";
}

sub toString_tsv{
  my($self,$sortBy)=@_;
  my $dbh=$self->{dbh};

  my $str="";

  my $sql=qq(
    SELECT GENOME1,GENOME2,DISTANCE
    FROM DISTANCE
  );
  if($sortBy eq 'abc'){
    $sql.="ORDER BY GENOME1,GENOME2 ASC";
  } elsif($sortBy eq 'rand'){
    $sql.="ORDER BY NEWID()";
  }
  
  my $sth=$dbh->prepare($sql);
  my $rv=$sth->execute or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  while(my @row=$sth->fetchrow_array()){
    $str.=join("\t",@row)."\n";
  }
  return $str;
}

sub toString_phylip{
  my($self,$sortBy)=@_;
  my $dbh=$self->{dbh};

  my $str="";

  # The way phylip is, I need to know the genome names
  # a priori
  my @name;
  my $sth=$dbh->prepare(qq(
    SELECT DISTINCT(GENOME1) 
    FROM DISTANCE 
    ORDER BY GENOME1 ASC
  ));
  my $rv=$sth->execute or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  my $maxGenomeLength=0;
  while(my @row=$sth->fetchrow_array()){
    push(@name,$row[0]);
    $maxGenomeLength=length($row[0]) if(length($row[0]) > $maxGenomeLength);
  }

  # We are already sorted alphabetically, so just worry
  # about whether or not we sort by random
  if($sortBy eq 'rand'){
    @name = shuffle(@name);
  }

  my $numGenomes=@name;

  $str.=(" " x 4) . "$numGenomes\n";
  for(my $i=0;$i<$numGenomes;$i++){
    $str.=$name[$i];
    $str.=" " x ($maxGenomeLength - length($name[$i]) + 2);
    my $distanceHash=$self->findDistances($name[$i]);

    for(my $j=0;$j<$numGenomes;$j++){
      $str.=sprintf("%0.4f  ",$$distanceHash{$name[$j]});
    }
    $str=~s/ +$/\n/; # replace that trailing whitespace with a newline
  }
  return $str;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


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
  my $sth = $dbh->prepare(qq(
    CREATE TABLE DISTANCE(
      GENOME1     CHAR(255)    NOT NULL,
      GENOME2     CHAR(255)    NOT NULL,
      DISTANCE    INT          NOT NULL,
      PRIMARY KEY(GENOME1,GENOME2)
    )) 
  );
  $sth->execute();

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

sub addDistancesFromHash{
  my($self,$distHash)=@_;

  my $dbh=$self->{dbh};
  my $numInserted=0;   # how many are going to be inserted?

  my $insert = $dbh->prepare( "INSERT INTO DISTANCE VALUES ( ?, ?, ? )" );
  my $autocommit = $dbh->{AutoCommit};
  $dbh->{AutoCommit} = 0; # begin a new transaction
  my $query="";
  #my $genome = map {s/^\s+|\s+$//g; _truncateFilename($_);} keys(%$distHash);
  my @genome = keys(%$distHash);
  my $numGenomes = @genome;
  for(my $i=0;$i<$numGenomes;$i++){
    my $genomeI = $genome[$i];
    for(my $j=$i;$j<$numGenomes;$j++){
      next if(defined($self->findDistance($genomeI,$genome[$j])));

      $insert->execute($genomeI, $genome[$j], $$distHash{$genomeI}{$genome[$j]});
      if ( $dbh->err() ) {
        die "Error: could not insert distance between $genomeI and $genome[$j] into the database: ".$dbh->err."\n";
      }
      $numInserted++;
    }
  }
  $dbh->commit;
  $dbh->{AutoCommit} = $autocommit;

  return $numInserted;
}

sub addDistances{
  my($self,$distancesFile)=@_;

  my $dbh=$self->{dbh};
  my $numInserted=0;   # how many are going to be inserted?

  my $insert = $dbh->prepare( "INSERT INTO DISTANCE VALUES ( ?, ?, ? )" );
  my $autocommit = $dbh->{AutoCommit};
  $dbh->{AutoCommit} = 0; # begin a new transaction
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

    $insert->execute( $query, $subject, $distance );
    if ( $dbh->err() ) {
        die "Error: could not insert $distancesFile into the database.\n";
    }
    $numInserted++;
  }
  $dbh->commit;
  $dbh->{AutoCommit} = $autocommit;

  return $numInserted;
}

# Find all distances with a single genome
sub findDistances{
  my($self,$genome1)=@_;

  my $dbh=$self->{dbh};
  
  my $sth=$dbh->prepare(qq(SELECT GENOME2,DISTANCE 
    FROM DISTANCE 
    WHERE GENOME1=?
    ORDER BY GENOME2
  ));
  my $rv = $sth->execute( $genome1 ) or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  # Distance will be undefined unless there is a result
  # on the SQL select statement.
  my %distance;
  while(my @row=$sth->fetchrow_array()){
    $distance{$row[0]}=$row[1];
  }

  # Find it in GENOME2 also
  my $sth2=$dbh->prepare(qq(SELECT GENOME1,DISTANCE 
    FROM DISTANCE 
    WHERE GENOME2=?
    ORDER BY GENOME1
  ));
  my $rv2 = $sth2->execute( $genome1 ) or die $DBI::errstr;
  if($rv2 < 0){
    die $DBI::errstr;
  }
  while(my @row=$sth2->fetchrow_array()){
    $distance{$row[0]}=$row[1];
  }

  return \%distance;
}

sub findDistance{
  my($self,$genome1,$genome2)=@_;

  my $dbh=$self->{dbh};
  
  my $sth=$dbh->prepare(qq(SELECT DISTANCE FROM DISTANCE WHERE GENOME1=? AND GENOME2=?));
  my $rv = $sth->execute( $genome1, $genome2 ) or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  # Distance will be undefined unless there is a result
  # on the SQL select statement.
  my $distance;
  while(my @row=$sth->fetchrow_array()){
    ($distance)=@row;
  }

  # Look in reverse order too if we haven't found a value
  if(!defined $distance){
    my $sth2=$dbh->prepare(qq(SELECT DISTANCE FROM DISTANCE WHERE GENOME1=? AND GENOME2=?));
    my $rv2 = $sth2->execute( $genome2, $genome1 ) or die $DBI::errstr;
    if($rv2 < 0){
      die $DBI::errstr;
    }
    while(my @row=$sth2->fetchrow_array()){
      ($distance)=@row;
    }
  }
  return $distance;
}

# Format can be:
#   tsv    3-column format
#   matrix all-vs all tsv format
#   phylip Phylip matrix format
# sortBy can be:
#   abc
#   rand
sub toString{
  my($self,$genome,$format,$sortBy)=@_;
  $format//="tsv";
  $format=lc($format);
  $sortBy//="abc";
  $sortBy=lc($sortBy);
  
  if($format eq "tsv"){
    return $self->toString_tsv($genome,$sortBy);
  } elsif($format eq "matrix"){
    return $self->toString_matrix($genome,$sortBy);
  } elsif($format eq "phylip"){
    return $self->toString_phylip($genome,$sortBy);
  }

  die "ERROR: could not format ".ref($self)." as $format.";
}

sub toString_matrix{
  my($self,$genome,$sortBy)=@_;

  my %distance=$self->toString_tsv($genome,$sortBy);
  my @name=keys(%distance);
  my $numNames=@name;

  my $str=join("\t",".",@name)."\n";
  for(my $i=0;$i<$numNames;$i++){
    $str.=$name[$i];
    for(my $j=0;$j<$numNames;$j++){
      my $dist = $distance{$name[$i]}{$name[$j]} || $distance{$name[$j]}{$name[$i]};
      $str.="\t$dist";
    }
    $str.="\n";
  }

  return $str;
}

sub toString_tsv{
  my($self,$genome,$sortBy)=@_;
  $sortBy||="abc";
  $genome||=[];

  my $dbh=$self->{dbh};

  # Index the genome array
  my %genome;
  $genome{_truncateFilename($_)}=1 for(@$genome);

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

  my %distance;
  while(my @row=$sth->fetchrow_array()){
    # If the parameter was given to filter genome names,
    # do it here.
    my $truncatedName1 = _truncateFilename($row[0]);
    my $truncatedName2 = _truncateFilename($row[1]);
    next if(@$genome && (!$genome{$truncatedName1} || !$genome{$truncatedName2}));

    $_=~s/^\s+|\s+$//g for(@row); # whitespace trim
    $str.=join("\t",@row)."\n";
    $distance{$row[0]}{$row[1]}=$row[2];
  }
  return %distance if(wantarray);
  return $str;
}

sub toString_phylip{
  my($self,$genome,$sortBy)=@_;
  $sortBy||="abc";
  $genome||=[];
  my $dbh=$self->{dbh};

  # Index the genome array
  my %genome;
  $genome{_truncateFilename($_)}=1 for(@$genome);

  my $str="";

  # The way phylip is, I need to know the genome names
  # a priori. Get the genome names from the db.
  my @name;
  my $sql=qq(
    SELECT DISTINCT(GENOME1) 
    FROM DISTANCE 
    ORDER BY GENOME1 ASC\n
  );
  my $sth=$dbh->prepare($sql);
  my $rv=$sth->execute or die $DBI::errstr;
  if($rv < 0){
    die $DBI::errstr;
  }

  my $maxGenomeLength=0;
  while(my @row=$sth->fetchrow_array()){
    # If the parameter was given to filter genome names,
    # do it here.
    my $truncatedName = _truncateFilename($row[0]);
    next if(@$genome && !$genome{$truncatedName});

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
      if(!defined($$distanceHash{$name[$j]})){
        logmsg "WARNING: could not find distance for $name[$i] and $name[$j]";
      }
      $str.=sprintf("%0.10f  ",$$distanceHash{$name[$j]});
    }
    $str=~s/ +$/\n/; # replace that trailing whitespace with a newline
  }
  return $str;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


#!/usr/bin/env perl
package Mashtree::Db;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use List::Util qw/shuffle/;
use Data::Dumper;

use lib dirname($INC{"Mashtree/Db.pm"});
use lib dirname($INC{"Mashtree/Db.pm"})."/..";

use Mashtree qw/_truncateFilename logmsg sortNames/;


our @EXPORT_OK = qw(
         );

local $0=basename $0;

=pod

=head1 NAME

Mashtree::Db - functions for Mashtree databasing 

=head1 SYNOPSIS

  use strict;
  use warnings
  use Mashtree::Db;

  my $dbFile = "mashtree.tsv";
  my $db=Mashtree::Db->new($dbFile);

  # Add 10 distances from genome "test" to other genomes
  my %distHash;
  for(my $dist=0;$dist<10;$dist++){
    my $otherGenome = "genome" . $dist;
    $distHash{"test"}{$otherGenome} = $dist;
  }
  $db->addDistancesFromHash(\$distHash);

  my $firstDistance = $db->findDistance("test", "genome0");
  # => 0

=head1 DESCRIPTION

This is a helper module, usually not used directly.
This is how Mashtree reads and writes to the internal database.

=cut

=pod

=head1 METHODS

=over

=item Mashtree::Db->new($dbFile, \%settings)

Create a new Mashtree::Db object.

The database file is a tab-separated file and will be created if it doesn't exist.
If it does exist, then it will be read into memory.

Arguments:

  * $dbFile - a file path
  * $settings - a hash of key/values (currently unused)

=back

=cut

# Properties of this object:
#   dbFile
#   cache
#   settings, a hashref with keys:
#     significant_figures, how many sigfigs in mash distances default:10
sub new{
  my($class,$dbFile,$settings)=@_;

  # How many significant digits to go into the mash dists 
  $$settings{significant_figures} ||= 10;

  my $self={
    _significant_figures => $$settings{significant_figures},
    cache                => {},
  };
  bless($self,$class);

  $self->selectDb($dbFile);
  $self->readDatabase();

  return $self;
}

=pod

=over 

=item $db->selectDb

Selects a database. If it doesn't exist, then it will be created.
Then, it sets the object property `dbFile` to the file path.

=back

=cut

# Create a database from a TSV
# Returns 1 if created a new database
sub selectDb{
  my($self, $dbFile)=@_;

  $self->{dbFile}=$dbFile;

  if(!-e $dbFile){
    open(my $fh, '>', $dbFile) or die "ERROR: could not write to $dbFile: $!";
    print $fh join("\t", qw(genome1 genome2 distance))."\n";
    close $fh;
  }

  return 1;
}

=pod

=over

=item $db->readDatabase

Reads the database from the dbFile set by `selectDb`.
Returns a hash of distances, e.g., genome1 => {genome2=>dist}

Then, this hash of distances is set in the object property `cache`.

=back

=cut

sub readDatabase{
  my($self) = @_;

  my %dist;

  open(my $fh, "<", $self->{dbFile}) or die "ERROR: could not read from ".$self->{dbFile}.": $!";
  my $header = <$fh>;
  chomp($header);
  my @header = split(/\t/, $header);
  my $numHeaders = @header;
  while(my $line = <$fh>){
    chomp($line);
    my @F = split(/\t/, $line);
    my %F;
    for(my $i=0;$i<$numHeaders;$i++){
      $F{$header[$i]} = $F[$i];
    }
    
    $dist{$F{genome1}}{$F{genome2}} = $F{distance};
  }

  $self->{cache} = \%dist;
  return \%dist;
}

=pod

=over

=item addDistancesFromHash

Add distances from a perl hash, $distHash
$distHash is { genome1 => {$genome2 => $dist} }

=back

=cut

sub addDistancesFromHash{
  my($self,$distHash)=@_;

  my $numInserted=0;   # how many are going to be inserted?

  open(my $fh, ">>", $self->{dbFile}) or die "ERROR: could not append to ".$self->{dbFile}.": $!";
  for my $genome1(sort keys(%$distHash)){
    for my $genome2(sort keys(%{ $$distHash{$genome1} })){
      my $dist = $$distHash{$genome1}{$genome2};
      if(!defined($dist)){
        logmsg "WARNING: distance was not defined for $genome1 <=> $genome2";
        logmsg "  Setting distance to 1.";
        $dist = 1;
      }
      print $fh join("\t", $genome1, $genome2, $dist)."\n";
      $numInserted++;
    }
  }
  close $fh;

  return $numInserted;
}

=pod 

=over

=item $db->addDistances

Add distances from a TSV file.
TSV file should be a mash distances tsv file and is in the format of, e.g.,
   # query t/lambda/sample1.fastq.gz
   t/lambda/sample2.fastq.gz  0.059
   t/lambda/sample3.fastq.gz  0.061

=back

=cut

sub addDistances{
  my($self,$distancesFile)=@_;

  # update the cache just in case
  # it helps us skip some redundancies.
  $self->readDatabase();

  my $numInserted=0;   # how many are going to be inserted?

  open(my $inFh, "<", $distancesFile) or die "ERROR: could not read $distancesFile: $!";
  open(my $outFh, ">>", $self->{dbFile}) or die "ERROR: could not append to ".$self->{dbFile}.": $!";
  my $query = "";
  while(<$inFh>){
    chomp;
    if(/^#\s*query\s+(.+)/){
      $query=$1;
      $query=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
      $query=_truncateFilename($query);
      next;
    }
    die "ERROR: query was not stated in the mash dist -t output: $distancesFile" if(!$query);
    my($subject,$distance)=split(/\t/,$_);
    $subject=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
    $subject=_truncateFilename($subject);

    next if(defined($self->findDistance($query,$subject)));
    print $outFh join("\t", $query, $subject, $distance)."\n";
    $numInserted++;
  }
  close $outFh;
  close $inFh;

  return $numInserted;
}

=pod

=over

=item $db->findDistance

Find the distance between any two genomes.
Return undef if not found.

=back

=cut

sub findDistance{
  my($self, $query, $subject) = @_;
  
  my $dists = $self->{cache};
  if(defined($$dists{$query}{$subject})){
    return $$dists{$query}{$subject};
  }
  if(defined($$dists{$subject}{$query})){
    return $$dists{$subject}{$query};
  }
  
  return undef;
}

=pod

=over

=item $db->findDistances

Find all distances from one genome to all others
Return undef if not found.

=back

=cut

sub findDistances{
  my($self, $query) = @_;

  my $dists = $self->{cache};
  return $$dists{$query};
}

=pod

=over

=item $db->toString

Turn the database into a string representation.

Arguments:

  * genomeArray - list of genomes to include, or undef for all genomes
  * format - can be a string of one of these values:
    * tsv    3-column format (default)
    * matrix all-vs all tsv format
    * phylip Phylip matrix format
  * sortBy - can be:
    * abc (default)
    * rand

=back

=cut

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

# Return the database in a TSV formatted file along with a header,
# or if you ask for an array/hash, then it will return a hash of genome1=>{genome2=>dist}
sub toString_tsv{
  my($self,$genome,$sortBy)=@_;
  $sortBy||="abc";
  $genome||=[];

  my $dbh=$self->{dbh};

  # Index the genome array
  my %genome;
  $genome{_truncateFilename($_)}=1 for(@$genome);

  my %dist;

  open(my $fh, "<", $self->{dbFile}) or die "ERROR: could not read ".$self->{dbFile}.": $!";
  my $str = <$fh>;
  my @header = split(/\t/, $str);
  chomp(@header);
  while(my $line = <$fh>){
    $str .= $line;

    my %d;
    chomp($line);
    my @F = split(/\t/, $line);
    @d{@header} = @F;

    $dist{$d{genome1}}{$d{genome2}}=$d{distance};
  }
  return %dist if(wantarray);
  return $str;
}

sub toString_phylip{
  my($self,$genome,$sortBy)=@_;

  my %dist = $self->toString_tsv();
  $sortBy||="abc";
  $genome||=[];

  my $sigfigs = $$self{_significant_figures} || die "INTERNAL ERROR: could not set sigfigs";

  # Index the genome array
  my %genome;
  $genome{_truncateFilename($_)}=1 for(@$genome);

  my $str="";

  # The way phylip is, I need to know the genome names
  # a priori. Get the genome names from the db.
  my @rawName = sort{$a cmp $b} keys(%dist);

  my $maxGenomeLength=0;
  my @name;
  for my $name(@rawName){
    # If the parameter was given to filter genome names,
    # do it here.
    my $truncatedName = _truncateFilename($name);
    next if(@$genome && !$genome{$truncatedName});

    push(@name,$name);
    $maxGenomeLength=length($name) if(length($name) > $maxGenomeLength);
  }

  # We are already sorted alphabetically, so just worry
  # about whether or not we sort by random
  if($sortBy eq 'rand'){
    @name = shuffle(@name);
  }

  my $numGenomes=@name;

  $str.=(" " x 4) . "$numGenomes\n";
  for(my $i=0;$i<$numGenomes;$i++){

    # sanitize sample names by removing spaces
    my $sanitizedName = $name[$i];
    $sanitizedName =~ s/\s/_/g;

    $str.=$sanitizedName;
    $str.=" " x ($maxGenomeLength - length($name[$i]) + 2);
    my $distanceHash = $dist{$name[$i]};

    for(my $j=0;$j<$numGenomes;$j++){
      if(!defined($$distanceHash{$name[$j]})){
        logmsg "WARNING: could not find distance for $name[$i] and $name[$j]";
      }
      $str.=sprintf("%0.${sigfigs}f  ",$$distanceHash{$name[$j]});
    }
    $str=~s/ +$/\n/; # replace that trailing whitespace with a newline
  }
  return $str;
}

1; # gotta love how we we return 1 in modules. TRUTH!!!


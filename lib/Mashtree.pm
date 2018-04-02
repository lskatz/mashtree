#!/usr/bin/env perl
package Mashtree;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use List::Util qw/shuffle/;

use threads;
use threads::shared;

use lib dirname($INC{"Mashtree.pm"});
use Bio::Matrix::IO;
use Bio::TreeIO;

our @EXPORT_OK = qw(
           logmsg openFastq _truncateFilename distancesToPhylip createTreeFromPhylip sortNames treeDist
           @fastqExt @fastaExt @bamExt @vcfExt @richseqExt @mshExt
           $MASHTREE_VERSION
         );

local $0=basename $0;

######
# CONSTANTS

our $VERSION = "0.31";
our $MASHTREE_VERSION=$VERSION;
our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fas .fsa .fa);
our @bamExt=qw(.sorted.bam .bam);
our @vcfExt=qw(.vcf.gz .vcf);
our @mshExt=qw(.msh);
# Richseq extensions were obtained mostly from bioperl under
# the genbank, embl, and swissprot entries, under
# the source for Bio::SeqIO
our @richseqExt=qw(.gb .gbank .genbank .gbk .gbs .gbf .embl .ebl .emb .dat .swiss .sp);

# Helpful things
my $fhStick :shared;  # A thread can only open a fastq file if it has the talking stick.

#################################################
### COMMON SUBS/TOOLS (not object subroutines) ##
#################################################
# Redefine how scripts die
$SIG{'__DIE__'} = sub {
  local $0=basename($0);
  my $e = $_[0] || ""; 
  my $callerSub=(caller(1))[3] || (caller(0))[3] || "UnknownSub";

  $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; 
  die("$0: $callerSub: $e"); 
};
# Centralized logmsg
#sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}
sub logmsg {
  local $0=basename $0;
  my $parentSub=(caller(1))[3] || (caller(0))[3];
  $parentSub=~s/^main:://;

  # Find the thread ID and stringify it
  my $tid=threads->tid;
  $tid=($tid) ? "(TID$tid)" : "";

  my $msg="$0: $parentSub$tid: @_\n";

  print STDERR $msg;
}

# Opens a fastq file in a thread-safe way.
sub openFastq{
  my($fastq,$settings)=@_;

  my $fh;

  lock($fhStick);

  my @fastqExt=qw(.fastq.gz .fastq .fq.gz .fq);
  my($name,$dir,$ext)=fileparse($fastq,@fastqExt);
  if($ext =~/\.gz$/){
    open($fh,"zcat $fastq | ") or die "ERROR: could not open $fastq for reading!: $!";
  } else {
    open($fh,"<",$fastq) or die "ERROR: could not open $fastq for reading!: $!";
  }
  return $fh;
}

# Removes fastq extension, removes directory name,
# truncates to a length, and adds right-padding.
sub _truncateFilename{
  my($file,$settings)=@_;
  $$settings{truncLength}||=255;
  # strip off msh and get the basename
  my $name=basename($file,'.msh');
  # One more extension
  $name=basename($name,@fastqExt,@richseqExt,@fastaExt);
  # Truncate
  $name=substr($name,0,$$settings{truncLength}); 
  # Add in padding
  $name.=" " x ($$settings{truncLength}-length($name)); 
  return $name;
}


# 1. Read the mash distances
# 2. Create a phylip file
sub distancesToPhylip{
  my($distances,$outdir,$settings)=@_;

  my $phylip = "$outdir/distances.phylip"; 
  # NOTE: need to regenerate the combined distances each time
  # because I need to allow variation in the input samples.
  #return $phylip if(-e $phylip);

  # The way phylip is, I need to know the genome names
  # a priori
  my %name;
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";
  while(<MASHDIST>){
    next if(/^#/);
    my($name)=split(/\t/,$_);
    $name=~s/^\s+|\s+$//g;  # whitespace trim before right-padding is added
    $name{_truncateFilename($name,$settings)}=1;
  }
  close MASHDIST;
  my @name=sortNames([keys(%name)],$settings);
  # Index the array
  my $columnIndex=0;
  for(@name){
    $name{$_}=$columnIndex++;
  }

  # Load up the matrix object
  logmsg "Reading the distances file at $distances";
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";
  my $query="UNKNOWN"; # Default ID in case anything goes wrong
  my @m;
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      $query=_truncateFilename($1,$settings);
    } else {
      my ($reference,$distance)=split(/\t/,$_);
      $reference=_truncateFilename($reference,$settings);
      $distance=sprintf("%0.10f",$distance);
      $m[$name{$query}][$name{$reference}]=$distance;
      $m[$name{$reference}][$name{$query}]=$distance;
    }
  }
  close MASHDIST;
  #my $matrixObj=Bio::Matrix::Generic->new(-rownames=>\@name,-colnames=>\@name,-values=>\@m);

  # taking this method from write_matrix in http://cpansearch.perl.org/src/CJFIELDS/BioPerl-1.6.924/Bio/Matrix/IO/phylip.pm
  my $str;
  $str.=(" " x 4) . scalar(@name)."\n";
  for(my $i=0;$i<@name;$i++){
    $str.=$name[$i];
    my $count=0;
    for(my $j=0;$j<@name;$j++){
      if($count < $#name){
        $str.=$m[$i][$j]. "  ";
      } else {
        $str.=$m[$i][$j];
      }
      $count++;
    }
    $str.="\n";
  }
  open(PHYLIP,">",$phylip) or die "ERROR: could not write to $phylip: $!";
  print PHYLIP $str;
  close PHYLIP;
  return $phylip;
}

sub sortNames{
  my($name,$settings)=@_;
  my @sorted;
  if($$settings{'sort-order'} =~ /^(abc|alphabet)$/){
    @sorted=sort { $a cmp $b } @$name;
  } elsif($$settings{'sort-order'}=~/^rand(om)?/){
    @sorted=shuffle(@$name);
  } elsif($$settings{'sort-order'} eq 'input-order'){
    @sorted=@$name;
  } else {
    die "ERROR: I don't understand sort-order $$settings{'sort-order'}";
  }
  return @sorted;
}

# Create tree file with Quicktree but bioperl 
# as a backup.
sub createTreeFromPhylip{
  my($phylip,$outdir,$settings)=@_;

  my $treeObj;

  my $quicktreePath=`which quicktree 2>/dev/null`;
  # bioperl if there was an error with which quicktree
  if($?){
    logmsg "DEPRECATION WARNING: CANNOT FIND QUICKTREE IN YOUR PATH. I will use BioPerl to make the tree this time, but it will be removed in the next version.";
    logmsg "Creating tree with BioPerl";
    my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
    my $matrix   = Bio::Matrix::IO->new(-format=>"phylip", -file=>$phylip)->next_matrix;
    $treeObj = $dfactory->make_tree($matrix);
    open(TREE,">","$outdir/tree.dnd") or die "ERROR: could not open $outdir/tree.dnd: $!";
    print TREE $treeObj->as_text("newick");
    print TREE "\n";
    close TREE;
  }
  # quicktree
  else {
    logmsg "Creating tree with QuickTree";
    system("quicktree -in m $phylip > $outdir/tree.dnd.tmp");
    die "ERROR with quicktree" if $?;
    $treeObj=Bio::TreeIO->new(-file=>"$outdir/tree.dnd.tmp")->next_tree;
    my $outtree=Bio::TreeIO->new(-file=>">$outdir/tree.dnd", -format=>"newick");
    $outtree->write_tree($treeObj);

    unlink("$outdir/tree.dnd.tmp");
  }

  return $treeObj;

}

# Lee's implementation of a tree distance. The objective
# is to return zero if two trees are the same.
sub treeDist{
  my($treeObj1,$treeObj2)=@_;

  # If the tree objects are really strings, then make Bio::Tree::Tree objects
  if(!ref($treeObj1)){
    if(-e $treeObj1){ # if this is a file, get the contents
      $treeObj1=`cat $treeObj1`;
    }
    $treeObj1=Bio::TreeIO->new(-string=>$treeObj1)->next_tree;
  }
  if(!ref($treeObj2)){
    if(-e $treeObj2){ # if this is a file, get the contents
      $treeObj2=`cat $treeObj2`;
    }
    $treeObj2=Bio::TreeIO->new(-string=>$treeObj2)->next_tree;
  }
  for($treeObj1,$treeObj2){
    #$_->force_binary;
  }
  
  # Get all leaf nodes so that they can be compared
  my @nodes1=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $treeObj1->get_nodes;
  my @nodes2=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $treeObj2->get_nodes;
  my $numNodes=@nodes1;

  # Test 1: are these the same nodes?
  my $nodeString1=join(" ",map{$_->id} @nodes1);
  my $nodeString2=join(" ",map{$_->id} @nodes2);
  if($nodeString1 ne $nodeString2){
    # TODO print out the differing nodes?
    logmsg "ERROR: nodes are not the same in both trees!\n  $nodeString1\n  $nodeString2";
    return ~0; #largest int
  }

  # Find the number of branches it takes to get to each node.
  # Turn it into a Euclidean distance
  my $euclideanDistance=0;
  for(my $i=0;$i<$numNodes;$i++){
    for(my $j=$i+1;$j<$numNodes;$j++){
      my ($numBranches1,$numBranches2);

      my $lca1=$treeObj1->get_lca($nodes1[$i],$nodes1[$j]);
      my $lca2=$treeObj2->get_lca($nodes2[$i],$nodes2[$j]);
      
      # Distance in tree1
      my $distance1=0;
      my @ancestory1=reverse $treeObj1->get_lineage_nodes($nodes1[$i]);
      my @ancestory2=reverse $treeObj1->get_lineage_nodes($nodes1[$j]);
      for my $currentNode(@ancestory1){
        $distance1++;
        last if($currentNode eq $lca1);
      }
      for my $currentNode(@ancestory2){
        $distance1++;
        last if($currentNode eq $lca1);
      }
      
      # Distance in tree2
      my $distance2=0;
      my @ancestory3=reverse $treeObj2->get_lineage_nodes($nodes2[$i]);
      my @ancestory4=reverse $treeObj2->get_lineage_nodes($nodes2[$j]);
      for my $currentNode(@ancestory3){
        $distance2++;
        last if($currentNode eq $lca2);
      }
      for my $currentNode(@ancestory4){
        $distance2++;
        last if($currentNode eq $lca2);
      }

      if($distance1 != $distance2){
        logmsg "These two nodes do not have the same distance between trees: ".$nodes1[$i]->id." and ".$nodes1[$j]->id;
      }

      # Add up the Euclidean distance
      $euclideanDistance+=($distance1 - $distance2) ** 2;
    }
  }
  $euclideanDistance=sqrt($euclideanDistance);
  return $euclideanDistance;
}

1;


__END__

=head1 NAME

Mashtree - Create a tree using Mash distances.

=cut


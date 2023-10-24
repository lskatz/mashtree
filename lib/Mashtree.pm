#!/usr/bin/env perl
package Mashtree;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use List::Util qw/shuffle/;
use Scalar::Util qw/looks_like_number/;

use threads;
use threads::shared;

use lib dirname($INC{"Mashtree.pm"});
use Bio::Matrix::IO;
use Bio::TreeIO;

our @EXPORT_OK = qw(
           logmsg openFastq _truncateFilename distancesToPhylip createTreeFromPhylip sortNames treeDist mashDist mashHashes raw_mash_distance raw_mash_distance_unequal_sizes
           @fastqExt @fastaExt @bamExt @vcfExt @richseqExt @mshExt
           $MASHTREE_VERSION
         );

local $0=basename $0;

=pod

=head1 NAME Mashtree

=head1 SYNOPSIS

Helps run a mashtree analysis to make rapid trees for genomes.
Please see github.com/lskatz/Mashtree for more information.

=over

=item mashtree executables

This document covers the Mashtree library, but the highlight
the mashtree package is the executable `mashtree`.
See github.com/lskatz/Mashtree for more information.

Fast method:

    mashtree --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

More accurate method:

    mashtree --mindepth 0 --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

Bootstrapping and jackknifing

    mashtree_bootstrap.pl --reps 100 --numcpus 12 *.fastq.gz -- --min-depth 0 > mashtree.jackknife.dnd
    mashtree_jackknife.pl --reps 100 --numcpus 12 *.fastq.gz -- --min-depth 0 > mashtree.jackknife.dnd

=back

=head1 VARIABLES

=over

=item $VERSION

=item $MASHTREE_VERSION (same value as $VERSION)

=item @fastqExt = qw(.fastq.gz .fastq .fq .fq.gz)

=item @fastaExt = qw(.fasta .fna .faa .mfa .fas .fsa .fa)

=item @bamExt = qw(.sorted.bam .bam)

=item @vcfExt = qw(.vcf.gz .vcf)

=item @mshExt = qw(.msh)

=item @richseqExt = qw(.gb .gbank .genbank .gbk .gbs .gbf .embl .ebl .emb .dat .swiss .sp)

=item $fhStick :shared 

Used to mark whether a file is being read, so that Mashtree limits disk I/O

=back

=cut

######
# CONSTANTS

our $VERSION = "1.4.6";
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

=head1 METHODS

=over

=item $SIG{'__DIE__'}

Remakes how `die` works, so that it references the caller

=cut

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

=pod

=item logmsg

Prints a message to STDERR with the thread number and the program name, with a trailing newline.

=cut

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

=pod

=item openFastq

 Opens a fastq file in a thread-safe way.

=cut

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

=pod

=item _truncateFilename

 Removes fastq extension, removes directory name,

=cut

sub _truncateFilename{
  my($file,$settings)=@_;
  # strip off msh and any other known extentions
  my $name=$file;
  my $oldName="";
  # Strip until we get convergence
  while($name ne $oldName){
    $oldName = $name;
    $name = basename($name,@mshExt,@fastqExt,@richseqExt,@fastaExt);
  }
  return $name;
}

=pod

=item distancesToPhylip

1. Read the mash distances
2. Create a phylip file

Arguments: hash of distances, output directory, settings hash

=cut

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

=pod

=item sortNames

Sorts names.

Arguments: 

1. $name - array of names
2. $settings - options
  * $$settings{'sort-order'} is either "abc", "random", "input-order"

=cut

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

=pod

=item createTreeFromPhylip($phylip, $outdir, $settings)

 Create tree file with Quicktree but bioperl 
 as a backup.

=cut

sub createTreeFromPhylip{
  my($phylip,$outdir,$settings)=@_;

  my $treeObj;

  my $quicktreePath=`which quicktree 2>/dev/null`;
  # bioperl if there was an error with which quicktree
  if($?){
    logmsg "DEPRECATION WARNING: CANNOT FIND QUICKTREE IN YOUR PATH. I will use BioPerl to make the tree this time, but it will be removed in the next version.";
    #logmsg "Creating tree with BioPerl";
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
    #logmsg "Creating tree with QuickTree";
    system("quicktree -in m $phylip > $outdir/tree.dnd.tmp");
    die "ERROR with quicktree" if $?;
    $treeObj=Bio::TreeIO->new(-file=>"$outdir/tree.dnd.tmp")->next_tree;
    open(my $treeFh, ">", "$outdir/tree.dnd") or die "ERROR: could not write to $outdir/tree.dnd: $!";
    print $treeFh $treeObj->as_text("newick")."\n";
    #my $outtree=Bio::TreeIO->new(-file=>">$outdir/tree.dnd", -format=>"newick");
    #$outtree->write_tree($treeObj);

    unlink("$outdir/tree.dnd.tmp");
  }

  return $treeObj;

}

=pod

=item treeDist($treeObj1, $treeObj2)

 Lee's implementation of a tree distance. The objective
 is to return zero if two trees are the same.

=cut

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

=pod

=item mashDist($file1, $file2, $k, $settings)

Find the distance between two mash sketch files
Alternatively: two hash lists.

=cut

sub mashDist{
  my($file1, $file2, $k, $settings)=@_;

  my($hashes1, $hashes2, $kmer1, $kmer2);
  if(ref($file1) eq 'ARRAY'){
    $hashes1 = $file1;
    $kmer1 = -1;
  } else {
    ($hashes1, $kmer1) = mashHashes($file1);
  }
  if(ref($file2) eq 'ARRAY'){
    $hashes2 = $file2;
    $kmer2 = -1;
  } else {
    ($hashes2, $kmer2) = mashHashes($file2);
  }

  if($kmer1 ne $kmer2){
    die "ERROR: kmer lengths do not match($kmer1 vs $kmer2)";
  }

  # Set the default kmer length and perform sanity check
  $k ||= $kmer1;
  if(!looks_like_number($k)){
    die "ERROR: k was not set to an integer";
  }
  if($k < 1){
    die "ERROR: k was undefined or set to less than 1";
  }

  my($common, $total) = (0,0);
  if(scalar(@$hashes1) != scalar(@$hashes2)){
    ($common, $total) = raw_mash_distance_unequal_sizes($hashes1, $hashes2);
  } else {
    ($common, $total) = raw_mash_distance($hashes1, $hashes2);
  }
  my $jaccard = $common/$total;
  my $mash_distance = -1/$k * log(2*$jaccard / (1+$jaccard));
  #logmsg "========== $mash_distance = -1/$k * log(2*$jaccard / (1+$jaccard)) ==============";

  return $mash_distance;
}

=pod

=item mashHashes($sketch)

Return an array of hashes, the kmer length, and the genome estimated length

=cut

sub mashHashes{
  my($sketch)=@_;
  my @hash;
  my $length = 0;
  my $kmer   = 0;

  if(!-e $sketch){
    die "ERROR: file not found: $sketch";
  }

  my $fh;
  if($sketch =~ /\.msh$/){
    open($fh, "mash info -d $sketch | ") or die "ERROR: could not run mash info -d on $sketch: $!";
  } elsif($sketch =~ /\.json$/){
    open($fh, $sketch) or die "ERROR: could not read $sketch: $!";
  }
  while(<$fh>){
    if(/kmer\D+(\d+)/){
      $kmer = $1;
    }
    elsif(/length\D+(\d+)/){
      $length = $1;
    }
    elsif(/hashes/){
      while(<$fh>){
        last if(/\]/);
        next if(!/\d/);
        s/\D+//g;
        s/^\s+|\s+$//g;
        push(@hash, $_);
      }
    }
  }
  if(!@hash){
    die "ERROR: no hashes found in $sketch";
  }
  return (\@hash, $kmer, $length);
}

=pod

=item raw_mash_distance_unequal_sizes($hashes1, $hashes2)

Compare unequal sized hashes. Treat the first
set of hashes as the reference (denominator)
set.

=cut

sub raw_mash_distance_unequal_sizes{
  my($hashes1, $hashes2) = @_;

  my (%sketch1,%sketch2);
  @sketch1{@$hashes1} = (1) x scalar(@$hashes1);
  @sketch2{@$hashes2} = (1) x scalar(@$hashes2);

  my %union;
  for my $h(@$hashes1){
    if($sketch2{$h}){
      $union{$h}++;
    }
  }

  my $common = scalar(keys(%union));
  my $total  = scalar(@$hashes1);

  return($common,$total);
}

=pod

=item raw_mash_distance($hashes1, $hashes2)

Return the number of kmers in common and the number compared total. inspiration from
https://github.com/onecodex/finch-rs/blob/master/src/distance.rs#L34

=cut

sub raw_mash_distance{
  my($hashes1, $hashes2) = @_;

  my @sketch1 = sort {$a <=> $b} @$hashes1;
  my @sketch2 = sort {$a <=> $b} @$hashes2;

  my $i      = 0;
  my $j      = 0;
  my $common = 0;
  my $total  = 0;

  my $sketch_size = @sketch1;
  while($total < $sketch_size && $i < @sketch1 && $j < @sketch2){
    my $ltgt = ($sketch1[$i] <=> $sketch2[$j]); # -1 if sketch1 is less than, +1 if sketch1 is greater than

    if($ltgt == -1){
      $i += 1;
    } elsif($ltgt == 1){
      $j += 1;
    } elsif($ltgt==0) {
      $i += 1;
      $j += 1;
      $common += 1;
    } else {
      die "Internal error";
    }

    $total += 1;
  }

  if($total < $sketch_size){
    if($i < @sketch1){
      $total += @sketch1 - 1;
    }

    if($j < @sketch2){
      $total += @sketch2 - 1;
    }

    if($total > $sketch_size){
      $total = $sketch_size;
    }
  }

  return ($common, $total);
}

# Calculates the Transfer Bootstrap Expectation (TBE) for internal nodes based on
# the methods outlined in Lemoine et al, Nature, 2018.
# Currently experimental.
# I entered this sub into bioperl > v1.7.5 but it is worth
# having locally here to help maintain compatibility.
# The only difference is that it isn't an object method
# and that it is called without an OO implementation.

=item transfer_bootstrap_expectation

 Title   : transfer_bootstrap_expectation
 Usage   : my $tree_with_bs = transfer_bootstrap_expectation(\@bs_trees,$guide_tree);
 Function: Calculates the Transfer Bootstrap Expectation (TBE) for internal nodes based on 
           the methods outlined in Lemoine et al, Nature, 2018.
           Currently experimental.
 Returns : L<Bio::Tree::TreeI>
 Args    : Arrayref of L<Bio::Tree::TreeI>s
           Guide tree, L<Bio::Tree::TreeI>s

=back

=cut

sub transfer_bootstrap_expectation{
  my ($bs_trees,$guide_tree) = @_;

  if(!defined($bs_trees) || ref($bs_trees) ne 'ARRAY'){
    die "ERROR: second parameter in assess_bootstrap() must be a list";
  }
  my $num_bs_trees = scalar(@$bs_trees);
  if($num_bs_trees < 1){
    die "ERROR: no bootstrap trees were passed to ".(caller(0))[3];
  }

  # internal nodes are defined by their children
  my %internal = ();
  my %leafNameId = ();
  my @idLookup = ();
  my @internalLookup = ();
  my @tree = ($guide_tree, @$bs_trees);
  my $numTrees = scalar(@tree);
  for(my $i = 0; $i < $numTrees; $i++){ # guide tree's index is $i==0
    # Do this as a top down approach, can probably be
    # improved by caching internal node states, but not going
    # to worry about it right now.

    my @allnodes = $tree[$i]->get_nodes;
    my @internalnodes = grep { ! $_->is_Leaf } @allnodes;
    for my $node ( @internalnodes ) {
      my @tips = sort map { $_->id } 
                      grep { $_->is_Leaf() } $node->get_all_Descendents;
      my $id = join(",", @tips);
      # Map the concatenated-leaf ID to the internal ID on the guide tree
      if( $i == 0 ) {
        $internal{$id} = $node->internal_id;
        $leafNameId{$node->internal_id} = $id;
      }

      # Record the tips for each tree's internal node
      # ID lookup (concatenated string of leaf names)
      $idLookup[$i]{$id} = \@tips;
      # Internal ID lookup
      $internalLookup[$i]{$internal{$id}} = \@tips;
    }
  }

  # Find the average distance from branch b to all
  # bootstrap trees' branches b*
  my @id = sort keys %internal;
  my $numIds = @id;
  # Loop through all internal nodes of the guide tree
  for(my $j=0; $j<$numIds; $j++){
    my $refNode = $guide_tree->find_node(-internal_id => $internal{$id[$j]});
    my $refNodeId = $refNode->internal_id;
    my $refJoinId = $leafNameId{$refNodeId};
    my $refLeaves = $idLookup[0]{$refJoinId};
    my %refLeafIndex = map{$_=>1} @$refLeaves;
    #next if(!defined($refLeaves));

    # For each internal node, start calculating for
    # an average TBE distance.
    my $nodeTbeTotal = 0;
  
    # Loop through all bootstrap trees, skipping the 0th
    # tree which is the guide tree.
    for(my $i=1;$i<$numTrees;$i++){

      # Find the right branch to bootstrap with. The right
      # branch will be the one that has the smallest
      # TBE distance.
      my @bsNode   = grep {!$_->is_Leaf} $tree[$i]->get_nodes;
      my $numBsIds = scalar(@bsNode);
      my $minDistance = ~0; # large int
      for(my $k=0;$k<$numBsIds;$k++){
        my @queryLeaves = sort map { $_->id }
                      grep { $_->is_Leaf() } $bsNode[$k]->get_all_Descendents;

        my %queryLeafIndex = map{$_=>1} @queryLeaves;

        # How many moves does it take to go from query to ref?
        my $dist=0;
        for my $queryLeaf(@queryLeaves){
          if(!$refLeafIndex{$queryLeaf}){
            $dist++;
          }
        }
        for my $refLeaf(@$refLeaves){
          if(!$queryLeafIndex{$refLeaf}){
            $dist++;
          }
        }

        if($dist < $minDistance){
          $minDistance = $dist;
        }
      }
      $nodeTbeTotal += $minDistance;
    }
    my $avgTbe = $nodeTbeTotal / $numTrees;
    
    # Calculate the average of all b to b* distances
    # But it is also 1 - average.
    my $numRefLeaves = scalar(@$refLeaves);
    my $nodeTbe = 1 - $avgTbe/$numRefLeaves;
    # Round to an integer
    $refNode->bootstrap(sprintf("%0.0f",100 * $nodeTbe));
  }
   
  return $guide_tree;
}

1;


__END__


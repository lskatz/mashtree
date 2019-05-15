#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Uses Mash and BioPerl to create a NJ tree based on distances.
# Run this script with -h for help and usage.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/tempdir tempfile/;
use File::Basename qw/basename dirname fileparse/;
use File::Copy qw/cp mv/;
use List::Util qw/shuffle/;
use List::MoreUtils qw/part/;
use POSIX qw/floor/;
use IO::Handle; # allows me to autoflush file handles

#use Fcntl qw/:flock LOCK_EX LOCK_UN/;

use threads;
use Thread::Queue;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Mashtree qw/logmsg @fastqExt @fastaExt createTreeFromPhylip mashDist/;
use Mashtree::Db;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Tree::DistanceFactory;
use Bio::Tree::Statistics;

my $writeStick :shared;
my $readStick  :shared; # limit disk IO by reading some files one at a time

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  my @wrapperOptions=qw(help outmatrix=s tempdir=s reps=i numcpus=i);
  GetOptions($settings,@wrapperOptions) or die $!;
  $$settings{reps}//=1;
  $$settings{numcpus}||=1;
  die usage() if($$settings{help});
  die usage() if(@ARGV < 1);

  $$settings{tempdir}||=tempdir("MASHTREE_BOOTSTRAP.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  logmsg "Temporary directory will be $$settings{tempdir}";

  if($$settings{reps} < 1){
    die "ERROR: no reps were given!";
  }
  if($$settings{reps} < 100){
    logmsg "WARNING: You have very few reps planned on this mashtree run. Recommended reps are at least 100.";
  }
  
  ## Catch some options that are not allowed to be passed
  # Tempdir: All mashtree temporary directories will be under the
  # wrapper's tempdir.
  if(grep(/^\-+tempdir$/,@ARGV) || grep(/^\-+t$/,@ARGV)){
    die "ERROR: tempdir was specified for mashtree but should be an option for $0";
  }
  # Numcpus: this needs to be specified in the wrapper and will
  # appropriately be transferred to the mashtree script
  if(grep(/^\-+numcpus$/,@ARGV) || grep(/^\-+n$/,@ARGV)){
    die "ERROR: numcpus was specified for mashtree but should be an option for $0";
  }
  # Outmatrix: the wrapper script needs to control where
  # the matrix goes because it can only have the outmatrix
  # for the observed run and not the replicates for speed's
  # sake.
  if(grep(/^\-+outmatrix$/,@ARGV) || grep(/^\-+o$/,@ARGV)){
    die "ERROR: outmatrix was specified for mashtree but should be an option for $0";
  }
  
  # Separate flagged options from reads in the mashtree options
  my @reads = ();
  my @mashOptions = ();
  for(my $i=0;$i<@ARGV;$i++){
    if(-e $ARGV[$i]){
      push(@reads, $ARGV[$i]);
    } else {
      push(@mashOptions, $ARGV[$i]);
    }
  }
  my $mashOptions=join(" ",@mashOptions);
  my $reads = join(" ", @reads);
  
  # Some filenames we'll expect
  my $observeddir="$$settings{tempdir}/observed";
  my $obsDistances="$observeddir/distances.phylip";
  my $observedTree="$$settings{tempdir}/observed.dnd";
  my $outmatrix="$$settings{tempdir}/observeddistances.tsv";
  my $mshList = "$$settings{tempdir}/mash.list";
  my $mergedMash = "$$settings{tempdir}/merged.msh";

  # Make the observed directory and run Mash
  logmsg "Running mashtree on full data (".scalar(@reads)." targets)";
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree --outmatrix $outmatrix.tmp --tempdir $observeddir --numcpus $$settings{numcpus} $mashOptions $reads > $observedTree.tmp");
  die if $?;
  mv("$observedTree.tmp",$observedTree) or die $?;
  mv("$outmatrix.tmp",$outmatrix) or die $?;

  # Round-robin sample the replicate number, so that it's
  # well-balanced.
  my @repNumber = (1..$$settings{reps});
  my @reps;
  for(my $i=0;$i<$$settings{reps};$i++){
    my $thrIndex = $i % $$settings{numcpus};
    push(@{$reps[$thrIndex]}, $i);
  }

  # Sample the mash sketches with replacement for rapid bootstrapping
  my @bsThread;
  for my $i(0..$$settings{numcpus}-1){
    my %settingsCopy = %$settings;
    $bsThread[$i] = threads->new(\&mashtreeRepWorker, $reps[$i], $reads, $mashOptions, \%settingsCopy);
  }

  my @bsTree;
  for my $t(@bsThread){
    my $treeArr = $t->join;
    for my $file(@$treeArr){
      my $treein = Bio::TreeIO->new(-file=>$file);
      while(my $tree=$treein->next_tree){
        push(@bsTree, $tree);
      }
    }
  }
  
  # Combine trees into a bootstrapped tree and write it 
  # to an output file. Then print it to stdout.
  logmsg "Adding bootstraps to tree";
  my $guideTree=Bio::TreeIO->new(-file=>"$observeddir/tree.dnd")->next_tree;
  my $bsTree=assess_bootstrap($guideTree,\@bsTree,$guideTree);
  for my $node($bsTree->get_nodes){
    next if($node->is_Leaf);
    my $id=$node->bootstrap||$node->id||0;
    $node->id($id);
  }
  open(my $treeFh,">","$$settings{tempdir}/bstree.dnd") or die "ERROR: could not write to $$settings{tempdir}/bstree.dnd: $!";
  print $treeFh $bsTree->as_text('newick');
  print $treeFh "\n";
  close $treeFh;

  system("cat $$settings{tempdir}/bstree.dnd"); die if $?;

  if($$settings{'outmatrix'}){
    cp($outmatrix,$$settings{'outmatrix'});
  }
  
  return 0;
}

sub mashtreeRepWorker{
  my($reps, $reads, $mashtreeOptions, $settings) = @_;

  my @tree;
  my $largeInt = ~0;

  for my $rep(@$reps){
    logmsg "Rep $rep";
    my $tempdir = $$settings{tempdir}."/$rep";
    my $seed = int(rand(4.29497e+09)); # apparently 4.29497e+09 is the largest int for mash sketch
    my $repMashtreeOptions .= " --seed $seed";
    my $stdout = `$FindBin::RealBin/mashtree --tempdir $tempdir --numcpus 1 $repMashtreeOptions $reads 2>&1`;
    die "ERROR on rep $rep\n$stdout" if $?;

    push(@tree, "$$settings{tempdir}/$rep/tree.dnd");
  }

  return \@tree;
}


# Fixed bootstrapping function from bioperl
# https://github.com/bioperl/bioperl-live/pull/304
sub assess_bootstrap{
   my ($self,$bs_trees,$guide_tree) = @_;
   my @consensus;

   if(!defined($bs_trees) || ref($bs_trees) ne 'ARRAY'){
     die "ERROR: second parameter in assess_bootstrap() must be a list";
   }
   my $num_bs_trees = scalar(@$bs_trees);
   if($num_bs_trees < 1){
     die "ERROR: no bootstrap trees were passed to assess_bootstrap()";
   }

   # internal nodes are defined by their children

   my (%lookup,%internal);
   my $i = 0;
   for my $tree ( $guide_tree, @$bs_trees ) {
       # Do this as a top down approach, can probably be
       # improved by caching internal node states, but not going
       # to worry about it right now.

       my @allnodes = $tree->get_nodes;
       my @internalnodes = grep { ! $_->is_Leaf } @allnodes;
       for my $node ( @internalnodes ) {
           my @tips = sort map { $_->id } 
                      grep { $_->is_Leaf() } $node->get_all_Descendents;
           my $id = "(".join(",", @tips).")";
           if( $i == 0 ) {
               $internal{$id} = $node->internal_id;
           } else { 
               $lookup{$id}++;
           }
       }
       $i++;
   }
   #my @save; # unsure why this variable is needed
   for my $l ( keys %lookup ) {
       if( defined $internal{$l} ) {#&& $lookup{$l} > $min_seen ) {
           my $intnode = $guide_tree->find_node(-internal_id => $internal{$l});
           $intnode->bootstrap(sprintf("%d",100 * $lookup{$l} / $num_bs_trees));
       }
   }
   return $guide_tree;
}

#######
# Utils
#######

sub usage{
  my $usage="$0: a wrapper around mashtree for bootstrapping.
  Usage: $0 [options] [-- mashtree options] *.fastq.gz *.fasta > tree.dnd
  --outmatrix          ''   Output file for distance matrix
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping.
                            Bootstrapping will only work on compressed fastq
                            files.
  --numcpus            1    This will be passed to mashtree and will
                            be used to multithread reps.
  
  --                        Used to separate options for $0 and mashtree
  MASHTREE OPTIONS:\n".
  # Print the mashtree options starting with numcpus,
  # skipping the tempdir option.
  `mashtree --help 2>&1 | grep -A 999 "TREE OPTIONS" | grep -v ^Stopped`;

  return $usage;
}


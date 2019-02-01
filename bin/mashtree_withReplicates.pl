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
use File::Slurp qw/read_file write_file/;
use List::Util qw/shuffle/;
use List::MoreUtils qw/part/;
use POSIX qw/floor/;
use JSON;

use Fcntl qw/:flock LOCK_EX/;

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
use Bio::Matrix::Generic;

my $writeStick :shared;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  my @wrapperOptions=qw(help sketch-size=i outmatrix=s tempdir=s reps=i numcpus=i);
  GetOptions($settings,@wrapperOptions) or die $!;
  $$settings{reps}||=0;
  $$settings{numcpus}||=1;
  $$settings{'sketch-size'}||=10000;
  die usage() if($$settings{help});
  die usage() if(@ARGV < 1);

  $$settings{tempdir}||=tempdir("MASHTREE_WRAPPER.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  logmsg "Temporary directory will be $$settings{tempdir}";

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
  if(grep(/^\-+sketch\-size$/,@ARGV)){
    die "ERROR: --sketch-size was specified for mashtree but should be an option for $0";
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

  # Make the observed directory and run Mash
  logmsg "Running mashtree on full data";
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree --outmatrix $outmatrix.tmp --tempdir $observeddir --numcpus $$settings{numcpus} $mashOptions $reads > $observedTree.tmp");
  die if $?;
  mv("$observedTree.tmp",$observedTree) or die $?;
  mv("$outmatrix.tmp",$outmatrix) or die $?;

  my @msh = glob("$observeddir/*.msh");

  # Sample the mash sketches with replacement for rapid bootstrapping
  my @bsThread;
  my @repNumber = (1..$$settings{reps});
  my $bs_per_thread = int($$settings{reps} / $$settings{numcpus}) + 1;
  for my $i(0..$$settings{numcpus}-1){
    my @reps = splice(@repNumber, 0, $bs_per_thread);
    $bsThread[$i] = threads->new(sub{
      my($reps) = @_;
      my @bootstrapTree;
      for my $rep(@$reps){
        my $subsampleDir = "$$settings{tempdir}/rep$rep";
        mkdir $subsampleDir;
        logmsg "Initializing rep $rep => $subsampleDir";
        my $treeFile = subsampleMashSketches(\@msh, $subsampleDir, $settings);
        push(@bootstrapTree, $treeFile);
      }
      return \@bootstrapTree;
    }, \@reps);
  }

  my @bsTree;
  for(@bsThread){
    for my $file(@{ $_->join }){
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

sub subsampleMashSketches{
  my($mshArr, $subsampleDir, $settings) = @_;

  my %dist;

  # Make JSON files that represent sampling with replacement
  # of mash hashes
  my $json = JSON->new;
  $json->utf8;           # If we only expect characters 0..255. Makes it fast.
  $json->allow_nonref;   # can convert a non-reference into its corresponding string
  $json->allow_blessed;  # encode method will not barf when it encounters a blessed reference 
  $json->pretty;         # enables indent, space_before and space_after 
  for my $msh(@$mshArr){
    # Make a lock on this msh file while we read it and
    # convert it to json.
    open(my $lockFh, ">", "$msh.lock") or die "ERROR: could not make lockfile $msh.lock: $!";
    flock($lockFh, LOCK_EX) or die "ERROR locking file $msh.lock: $!";
    my $jsonText = `mash info -d $msh`;
    close $lockFh;

    # Manipulate the json data
    my $jsonHash = $json->decode($jsonText);
    my $numHashes = scalar(@{ $$jsonHash{sketches}[0]{hashes} });
    my @randHashesWithReplacement;
    for(1..$numHashes){
      my $randIndex = int(rand($numHashes));
      push(@randHashesWithReplacement,
        $$jsonHash{sketches}[0]{hashes}[$randIndex]
      );
    }
    $$jsonHash{sketches}[0]{hashes} = \@randHashesWithReplacement;
    
    # Write the json data to file
    write_file("$subsampleDir/".basename($msh).'.json', $json->encode($jsonHash));
  }

  # Make distances between all .msh.json files.
  my @jsonFile = glob("$subsampleDir/*.json");
  my @name = map {basename($_,qw(.fastq.gz.msh.json))} @jsonFile;
  my $matrix = Bio::Matrix::Generic->new(
    -rownames => \@name,
    -colnames => \@name
  );
  for(my $i=0;$i<@jsonFile;$i++){
    my $nameI = $name[$i];
    $matrix->entry($nameI,$nameI,0.0);
    for(my $j=$i+1; $j<@jsonFile;$j++){
      my $nameJ = $name[$j];
      my $distance = mashDist($jsonFile[$i], $jsonFile[$j]);
      $matrix->entry($nameI,$nameJ,$distance);
      $matrix->entry($nameJ,$nameI,$distance);
    }
  }

  my $numNames = @name;
  my $phylip = '    ' . $numNames ."\n";
  my $phylipFile = "$subsampleDir/distances.phylip";
  for(my $i=0;$i<$numNames;$i++){
    $phylip .= $name[$i];
    for(my $j=0;$j<$numNames;$j++){
      my $distance = sprintf("%0.10f", $matrix->entry($name[$i], $name[$j]));
      $phylip .= '  ' . $distance;
    }
    $phylip .= "\n";
  }

  open(my $phylipFh, ">", $phylipFile) or die "ERROR: could not write to $phylipFile: $!";
  print $phylipFh $phylip;
  close $phylipFh;

  my $treeObj = createTreeFromPhylip($phylipFile, $subsampleDir, $settings);

  return "$subsampleDir/tree.dnd";
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
  my $usage="$0: a wrapper around mashtree.
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


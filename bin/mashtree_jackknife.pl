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
# if JSON::XS is installed, it will be automatically used.
# Specify () so that no functions are exported. We will use OO.
use JSON (); 
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

  $$settings{tempdir}||=tempdir("MASHTREE_WRAPPER.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  logmsg "Temporary directory will be $$settings{tempdir}";

  if($$settings{reps} < 1){
    die "ERROR: no reps were given!";
  }
  if($$settings{reps} < 100){
    logmsg "WARNING: You have very few reps planned on this mashtree run. Recommended reps are at least 100.";
  }
  # Give a warning if JSON is going to be slow later
  my $jsonTmp = JSON->new();
  if(! $jsonTmp->is_xs){
    logmsg "WARNING: the currently installed JSON module will make this script very slow when jack knifing. To avoid this error, install the JSON::XS module like so: `cpanm ~l ~ JSON::XS`.";
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
  my $mergedJSON = "$$settings{tempdir}/merged.msh.json.gz";

  # Make the observed directory and run Mash
  logmsg "Running mashtree on full data (".scalar(@reads)." targets)";
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree --outmatrix $outmatrix.tmp --tempdir $observeddir --numcpus $$settings{numcpus} $mashOptions $reads > $observedTree.tmp");
  die if $?;
  mv("$observedTree.tmp",$observedTree) or die $?;
  mv("$outmatrix.tmp",$outmatrix) or die $?;

  # Merge the mash files to make the threads go faster later
  my @msh = glob("$observeddir/*.msh");
  open(my $mshListFh, ">", $mshList) or die "ERROR writing to mash list $mshList: $!";
  for(@msh){
    print $mshListFh $_."\n";
  }
  close $mshListFh;
  unlink($mergedMash) if(-e $mergedMash); # mash complains about overwriting an existing file
  system("mash paste -l $mergedMash $mshList >&2"); die "ERROR merging mash files" if $?; 
  # max compression on the json file so that it can be a
  # smaller footprint and so that it can be read faster
  # in each thread.
  system("mash info -d $mergedMash | gzip -c9 > $mergedJSON");
  die "ERROR with mash info | gzip -c9" if $?;
  unlink($_) for(@msh); # remove redundant files to the merged msh
  
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
    $bsThread[$i] = threads->new(\&subsampleMashSketchesWorker, $mergedJSON, $reps[$i], \%settingsCopy);
  }

  my @bsTree;
  for my $thr(@bsThread){
    my $fileArr = $thr->join;
    if(ref($fileArr) ne 'ARRAY'){
      die "ERROR: one or more threads did not return an array of jack knife tree files as expected.";
    }
    for my $file(@$fileArr){
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

sub subsampleMashSketchesWorker{
  my($mergedJSON, $reps, $settings) = @_;
  $reps //= [];

  my @treeFile;
  return \@treeFile if(@$reps <1);

  # Initialize our JSON reader
  my $json = JSON->new;
  $json->utf8;           # If we only expect characters 0..255. Makes it fast.
  $json->allow_nonref;   # can convert a non-reference into its corresponding string
  $json->allow_blessed;  # encode method will not barf when it encounters a blessed reference 
  $json->pretty;         # enables indent, space_before and space_after 
  # Only let one thread at a time read the large JSON file
  my $mashInfoStr=""; # must be declared outside of the block because we are parsing it after the block
  {
    lock($readStick);
    logmsg "Reading huge JSON file describing all mash distances, $mergedJSON";
    $mashInfoStr = `gzip -cd $mergedJSON`; 
    die "ERROR running gzip -cd $mergedJSON: $!" if $?;
  }
  my $mashInfoHash = $json->decode($mashInfoStr);
  $mashInfoStr=""; # clear some ram
  my $kmerlength = $$mashInfoHash{kmer}; # find the kmer length right away

  my $numReps = @$reps;
  my $numSketches = scalar(@{ $$mashInfoHash{sketches} });
  my $msg="Initializing thread with $numReps replicates: ";
  for my $i(0..2){
    my $repID = $$reps[$i] // "";
    $msg.="$repID, ";
  }
  $msg=~s/[ ,]+$//; # right trim
  if($numReps > 3){
    $msg.="...";
  }
  logmsg $msg;
  for(my $repI=0;$repI<@$reps;$repI++){
    my $rep = $$reps[$repI];

    my $subsampleDir = "$$settings{tempdir}/rep$rep";
    mkdir $subsampleDir;
    my $log = "$subsampleDir/jackknife.log";
    open(my $logFh, ">", $log) or die "ERROR: cannot write to log $log: $!";
    $logFh->autoflush();
    logmsg "rep$rep: $subsampleDir (".($repI+1)." out of $numReps in this thread)";

    my @name;
    # Start off a distances string for printing to file later
    my %dist = ();
    for(my $sketchCounter=0; $sketchCounter<$numSketches; $sketchCounter++){
      my $nameI = basename($$mashInfoHash{sketches}[$sketchCounter]{name},(@fastqExt,@fastaExt));
      print $logFh "Subsampling sketches from entry $sketchCounter, $nameI\n";

      # subsample the hashes of one genome at a time as compared
      # to the other genomes, to get a jack knife distance.
      my $numHashes = scalar(@{ $$mashInfoHash{sketches}[$sketchCounter]{hashes} });
      my $keepHashes = int($numHashes / 2); # based on half the hashes
      my @subsampleHash = @{ $$mashInfoHash{sketches}[$sketchCounter]{hashes} };
      @subsampleHash = (shuffle(@subsampleHash))[0..$keepHashes-1];
      @subsampleHash = sort{$a<=>$b} @subsampleHash;

      print $logFh "Distances between $nameI and other genomes";
      # Initialize the distances string with the query name.
      # Use this string as a buffer for the distances file
      # to help avoid too much disk IO.
      push(@name, $nameI);
      for(my $j=0; $j<$numSketches; $j++){
        my $nameJ = basename($$mashInfoHash{sketches}[$j]{name},(@fastqExt,@fastaExt));
        my $distance = mashDist(\@subsampleHash, $$mashInfoHash{sketches}[$j]{hashes}, $kmerlength);
        $dist{$nameI}{$nameJ} = $distance;
        print $logFh ".";
      }
      print $logFh "\n";
    }

    # Add distances to database
    print $logFh "Creating database, $subsampleDir/distances.sqlite\n";
    my $mashtreeDb = Mashtree::Db->new("$subsampleDir/distances.sqlite");
    $mashtreeDb->addDistancesFromHash(\%dist);
    # Convert to Phylip
    my $phylipFile = "$subsampleDir/distances.phylip";
    print $logFh "Creating phylip file, $phylipFile\n";
    open(my $phylipFh, ">", $phylipFile) or die "ERROR: could not write to $phylipFile: $!";
    print $phylipFh $mashtreeDb->toString(\@name, "phylip");
    close $phylipFh;

    print $logFh "Creating tree file, $subsampleDir/tree.dnd\n";
    my $treeObj = createTreeFromPhylip($phylipFile, $subsampleDir, $settings);
    push(@treeFile, "$subsampleDir/tree.dnd");

    print $logFh "DONE!\n";
    close $logFh;
  }

  return \@treeFile;
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


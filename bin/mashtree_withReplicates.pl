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
use Mashtree qw/logmsg @fastqExt @fastaExt createTreeFromPhylip/;
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

  if($$settings{reps} < 10){
    logmsg "WARNING: You have very few reps planned on this mashtree run. Recommended reps are at least 10 or 100.";
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

  # New idea!  Mess around with a large sketch file
  my $sketchPoolDir = makeSketchPool(\@reads, $settings);
  my @bsTree;
  for my $rep(0 .. $$settings{reps}-1){
    my $subsampleDir = "$$settings{tempdir}/rep$rep";
    mkdir $subsampleDir;
    my $jackknifeTree = subsampleMashSketches($sketchPoolDir, $subsampleDir, $settings);
    my $treeObj = Bio::TreeIO->new(-file=>$jackknifeTree)->next_tree;
    push(@bsTree, $treeObj);
  }
  
  # Make the observed directory and run Mash
  logmsg "Running mashtree on full data";
  mkdir($observeddir);
  system("$FindBin::RealBin/mashtree --outmatrix $outmatrix.tmp --tempdir $observeddir --numcpus $$settings{numcpus} $mashOptions $reads > $observedTree.tmp");
  die if $?;
  mv("$observedTree.tmp",$observedTree) or die $?;
  mv("$outmatrix.tmp",$outmatrix) or die $?;

  # Combine trees into a bootstrapped tree and write it 
  # to an output file. Then print it to stdout.
  logmsg "Adding bootstraps to tree";
  my $biostat=Bio::Tree::Statistics->new;
  my $guideTree=Bio::TreeIO->new(-file=>"$observeddir/tree.dnd")->next_tree;
  my $rootLeaf = (sort map {$_->id} grep {$_->is_Leaf}  $guideTree->get_nodes)[0];
  my $sample1Node = (grep{defined($_->id) && $_->id eq $rootLeaf} $guideTree->get_nodes)[0];
  $guideTree->reroot($sample1Node);
  my $bsTree=$biostat->assess_bootstrap(\@bsTree,$guideTree);
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

sub makeSketchPool{
  my($reads, $settings)=@_;

  my @readsCopy = @$reads;

  my @reads_per_thread = ([@readsCopy]);
  my $tmp_i = 0;
  if($$settings{numcpus} > 1){
    @reads_per_thread = part { $tmp_i++ % ($$settings{numcpus}-1) } @readsCopy;
  }

  my $largeSketchDir = "$$settings{tempdir}/largeSketches";
  mkdir $largeSketchDir;
  # Make the sketch size slightly more than double because
  # some will be truncated.
  my $sketchSize = int($$settings{'sketch-size'} * 2.0);

  my @thr;
  for my $i(0..$$settings{numcpus}-1){
    $thr[$i] = threads->new(sub{
      my($reads) = @_;
      for my $r(@$reads){
        my $prefix = "$largeSketchDir/".basename($r);
        system("mash sketch -s $sketchSize -o $prefix $r");
        die if $?;
      }
    }, $reads_per_thread[$i]);
  }

  for(@thr){
    $_->join;
  }

  # Truncate all sketches to the largest in-common int.
  # This will remove any specificity issues for integers
  # that would have fallen outside of the in-common range
  # of integers.
  # First find that largest integer.
  my $json = JSON->new;
  $json->allow_nonref;
  $json->allow_blessed;
  my $largestInt = ~0; # largest int to start
  for my $sketch(glob("$largeSketchDir/*.msh")){
    my $jsonText = `mash info -d $sketch`; die if $?;
    my $jsonHash = $json->decode($jsonText);
    #my $jsonHash = from_json($json);
    my @ints = sort{$b<=>$a} @{ $$jsonHash{sketches}[0]{hashes} };
    if($ints[0] < $largestInt){
      $largestInt = $ints[0];
    }
  }
  # Truncate all sketches to largest in-common int.
  for my $sketch(glob("$largeSketchDir/*.msh")){
    my $jsonText = `mash info -d $sketch`; die if $?;
    my $jsonHash = $json->decode($jsonText);
    my $ints = $$jsonHash{sketches}[0]{hashes};
    my @newInts;
    for(my $i=0;$i<@$ints;$i++){
      if($$ints[$i] <= $largestInt){
        push(@newInts, $$ints[$i]);
      }
    }
    $$jsonHash{sketches}[0]{hashes} = \@newInts;
    #logmsg scalar(@$ints)." => ".scalar(@newInts);

    open(my $fh, ">", "$sketch.json") or die "ERROR writing to $sketch.json: $!";
    print $fh $json->encode($jsonHash);
    close $fh;
  }

  return $largeSketchDir;
}

sub subsampleMashSketches{
  my($sketchPoolDir, $subsampleDir, $settings) = @_;
  
  my %dist;

  my $json = JSON->new;
  $json->allow_nonref;
  $json->allow_blessed;

  # Make distances:
  # Random integers 10k hashes from reference file.
  # Compare to see if they exist in other files.
  my @jsonFile = glob("$sketchPoolDir/*.msh.json");
  my @name = map {basename($_,qw(.fastq.gz.msh.json))} @jsonFile;
  my $matrix = Bio::Matrix::Generic->new(
    -rownames => \@name,
    -colnames => \@name
  );
  for(my $i=0;$i<@jsonFile;$i++){
    my $jsonText = read_file($jsonFile[$i]);
    my $jsonHash = $json->decode($jsonText);
    my $intI = $$jsonHash{sketches}[0]{hashes};
    my @randInts = shuffle(@$intI);
    splice(@randInts,$$settings{'sketch-size'});
    my %intI;
    @intI{@randInts} = (1) x scalar(@randInts);
    $matrix->entry($name[$i],$name[$i],0); # diagonal
    for(my $j=$i+1; $j<@jsonFile;$j++){
      $jsonText = read_file($jsonFile[$j]);
      $jsonHash = $json->decode($jsonText);
      my $distance = 0;
      for my $int(@{ $$jsonHash{sketches}[0]{hashes} }){
        if($intI{$int}){
          $distance++;
        }
      }
      $distance = $distance / scalar(@randInts);
      logmsg $distance, scalar(@randInts);
      $matrix->entry($name[$i],$name[$j],$distance);
      $matrix->entry($name[$j],$name[$i],$distance);
    }
  }

  if(0){
    print "$subsampleDir\n";
    for(my $i=0;$i<@name;$i++){
      print $name[$i];
      for(my $j=0; $j<@name;$j++){
        print "\t".$matrix->get_entry($name[$i],$name[$j]);
      }
      print "\n";
    }
    print "\n";
  }

  # Make tree from distances
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $treeObj = $dfactory->make_tree($matrix);

  # Reroot to ensure bootstrapping (not sure if this is necessary but just to be sure)
  my $rootLeaf = (sort @name)[0];
  my $sample1Node = (grep{defined($_->id) && $_->id eq $rootLeaf} $treeObj->get_nodes)[0];
  $treeObj->reroot($sample1Node);
  print $treeObj->as_text("newick")."\n";
  # Print tree to file
  open(my $treeFh, ">", "$subsampleDir/tree.dnd") or die "ERROR writing to $subsampleDir/tree.dnd: $!";
  print $treeFh $treeObj->as_text("newick")."\n";
  close $treeFh;
  
  return "$subsampleDir/tree.dnd";
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


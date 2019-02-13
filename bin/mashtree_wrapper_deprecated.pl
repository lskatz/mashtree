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
use Bio::Matrix::IO;

my $writeStick :shared;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  my @wrapperOptions=qw(help outmatrix=s tempdir=s reps=i numcpus=i);
  GetOptions($settings,@wrapperOptions) or die $!;
  $$settings{reps}||=0;
  $$settings{numcpus}||=1;
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

  # Copy reads over to the temp storage where I assume it is faster
  # and where we can write .lock files.
  my $inputdir = "$$settings{tempdir}/input";
  mkdir $inputdir;
  my $tmp_i= 0;
  my @reads_per_thread = ([@reads]);
  if($$settings{numcpus} > 1){
    @reads_per_thread = part { $tmp_i++ % ($$settings{numcpus}-1) } @reads;
  }
  my @cpThread;
  for(0..$$settings{numcpus}-1){
    $cpThread[$_] = threads->new(sub{
      my($fileArr) = @_;
      my @copiedReads;
      for my $file(@$fileArr){
        logmsg "Copying $file to temp space - $inputdir";
        my $copiedFile = "$inputdir/".basename($file);
        cp($file, $copiedFile);
        push(@copiedReads, $copiedFile);
      }
      return \@copiedReads;
    }, $reads_per_thread[$_]);
  }
  @reads = ();
  for(@cpThread){
    my $tmp = $_->join;
    push(@reads, @$tmp);
  }
  logmsg "Finished copying input data to $inputdir";

  my $mashOptions=join(" ",@mashOptions);
  my $reads = join(" ", @reads);
  
  # Some filenames we'll expect
  my $observeddir="$$settings{tempdir}/observed";
  my $obsDistances="$observeddir/distances.phylip";
  my $observedTree="$$settings{tempdir}/observed.dnd";
  my $outmatrix="$$settings{tempdir}/observeddistances.tsv";

  # Multithreaded reps
  my @rep_id = (1..$$settings{reps});
  my $repsPerThread = int($$settings{reps} / $$settings{numcpus}) + 1;
  my @thr;
  for(0..$$settings{numcpus}-1){
    my @theseReps = splice(@rep_id, 0, $repsPerThread);
    $thr[$_]=threads->new(\&repWorker, \@mashOptions, \@reads, \@theseReps, $settings);
  }

  my @bsTree;
  for(@thr){
    my $treeArr=$_->join;
    for(@$treeArr){
      push(@bsTree,Bio::TreeIO->new(-file=>$_)->next_tree);
    }
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

sub repWorker{
  my($mashOptions, $reads, $reps,$settings)=@_;
  my @bsTree;
  if(!defined($reps) || ref($reps) ne 'ARRAY' || !@$reps){
    return \@bsTree;
  }

  my @mashOptions = @$mashOptions;
  my @reads       = @$reads;

  my $numcpus = floor($$settings{numcpus}/$$settings{reps});
  $numcpus = 1 if($numcpus < 1);

  for my $rep(@$reps){
    my $repTempdir="$$settings{tempdir}/rep$rep";
    mkdir $repTempdir;
    logmsg "Starting mashtree replicate $rep - $repTempdir";
    
    #logmsg "Downsampling reads (replicate $rep).";
    # Downsample the reads
    my @downsampledReads=();
    for my $r(@reads){
      my $newReads = "$repTempdir/".basename($r);
      #logmsg "DEBUG";push(@downsampledReads, $newReads);next;

      my @buffer = ();
      open(my $lockFh, ">", "$r.lock") or die "ERROR: could not make lockfile $r.lock: $!";
      flock($lockFh, LOCK_EX) or die "ERROR locking file $r.lock: $!";

      open(my $inFh, "zcat $r | ") or die "ERROR reading $r for downsampling: $!";
      open(my $outFh," | gzip -c > $newReads") or die "ERROR gzipping to $newReads: $!";
      close $outFh;
      while(my $id=<$inFh>){
        my $seq =<$inFh>;
        my $plus=<$inFh>;
        my $qual=<$inFh>;
        if(rand(1) < 0.5){
          push(@buffer, $id.$seq.$plus.$qual);
          # The buffer size is something like 600 bytes per entry
          # times 100,000 entries times 12 threads = 720Mb.
          if(scalar(@buffer) > 100000){
            {
              # Only let one thread at a time write to the file but
              # this makes us open the file in append mode.
              lock($writeStick);
              open($outFh, " | gzip -c >> $newReads") or die "ERROR gzipping to $newReads: $!";
              print $outFh join("",@buffer);
              close $outFh;
            }
            @buffer = (); # flush the buffer
          }
        }
      } 
      # Finish the remaining entries in the buffer
      {
        lock($writeStick);
        open($outFh, " | gzip -c >> $newReads") or die "ERROR gzipping to $newReads: $!";
        print $outFh join("",@buffer);
        close $outFh;
      }
      close $inFh;
      close $lockFh;

      push(@downsampledReads, $newReads);
    }

    logmsg "Done downsampling for replicate $rep. Running mashtree on files in $repTempdir";
    my $log = `mashtree --numcpus $numcpus @mashOptions @downsampledReads 2>&1 > $repTempdir/tree.dnd`;
    if($?){
      die "ERROR with mashtree on rep $rep (exit code $?):\n$log";
    }

    logmsg "Finished with rep $rep";
    push(@bsTree,"$repTempdir/tree.dnd");
  }

  return \@bsTree;
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


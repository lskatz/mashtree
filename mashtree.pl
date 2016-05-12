#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Uses Mash and BioPerl to create a NJ tree based on distances.
# Run this script with -h for help and usage.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use File::Temp qw/tempdir tempfile/;
use File::Basename qw/basename dirname fileparse/;
use Bio::Tree::DistanceFactory;
use Bio::Matrix::IO;
use Bio::Tree::Statistics;

use threads;
use Thread::Queue;
use threads::shared;

local $0=basename $0;
my @fastqExt=qw(.fastq.gz .fastq .fq.gz .fq);
my $fhStick :shared;  # A thread can only open a fastq file if it has the talking stick.

# Smart log messaging with thread IDs
sub logmsg{ 
  my $tid=threads->tid || 0;
  if($tid > 0){
    $tid = " (TID$tid)";
  }
  print STDERR "$0$tid: @_\n";
}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i genomesize=i mindepth=i reps=i truncLength=i warn-on-duplicate validate-reads save-space)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{truncLength}||=250;  # how long a genome name is
  $$settings{reps}||=0;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  logmsg "Temporary directory will be $$settings{tempdir}";

  # Mash-specific options
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}||=2;

  die usage() if($$settings{help});

  my @reads=@ARGV;
  die usage() if(@reads < 2);

  # Check for prereq executables.
  for my $exe(qw(mash)){
    system("$exe -h > /dev/null 2>&1");
    die "ERROR: could not find $exe in your PATH" if $?;
  }

  logmsg "$0 on ".scalar(@reads)." files";

  validateFastq(\@reads,$settings) if($$settings{'validate-reads'});

  # This step will return an empty array list if there are no reps
  my $repsDirs=makeBootstrapReads(\@reads,$$settings{reps},$settings);

  my $primarySketches=sketchAll(\@reads,"$$settings{tempdir}/msh",$settings);

  my @bsSketches;
  for my $rep(@$repsDirs){
    my $bsSketches=sketchAll([glob("$rep/*.fastq")],$rep,$settings);
    push(@bsSketches,$bsSketches);
  }
  

  # Now that the sketches are all done, do the same steps on both
  # bootstrap samples and the real sample set.
  my @trees;
  my @sketches=($primarySketches, @bsSketches);
  for(my $i=0;$i<@sketches;$i++){

    my $subTempdir="$$settings{tempdir}/rep$i";
    mkdir $subTempdir;

    my $distances=mashDistance($sketches[$i],$subTempdir,$settings);

    my $phylip = distancesToPhylip($distances,$subTempdir,$settings);

    my $treeObj = createTree($phylip,$subTempdir,$settings);

    push(@trees,$treeObj);
  }


  # Make bootstraps but move the ID to the bootstrap field for
  # compatibility with Newick and tree drawing programs like MEGA.
  # TODO: move this to a subroutine to keep main() clean.
  my $guideTree=shift(@trees);
  my $stat=Bio::Tree::Statistics->new;
  my $bs_tree=$stat->assess_bootstrap(\@trees,$guideTree);
  for my $node(grep { ! $_->is_Leaf } $bs_tree->get_nodes){
    my $id=$node->bootstrap || 0;
    $node->id($id);
  }

  print $bs_tree->as_text('newick');
  
  return 0;
}

sub validateFastq{
  my($reads,$settings)=@_;

  logmsg "Validating your read sets";
  
  # Check whether the filenames are unique enough
  my %seen;
  for my $r(@$reads){
    my $trunc=_truncateFilename($r,$settings);
    if($seen{$trunc}){
      my $msg="I have already seen $r as $seen{$trunc} (truncated name: $trunc)";
      if($$settings{'warn-on-duplicate'}){
        logmsg "WARNING: $msg";
      } else {
        die "ERROR: $msg";
      }
    }
    $seen{$trunc}=$r;
  }

  # Check whether there could be enough reads in the read set.
  # Do this by going to the 10k-th line and seeing if we've
  # reached the end of the file yet.
  for my $r(@$reads){
    my $fastqFh=openFastq($r,$settings);
    my $numLines=0;
    while(<$fastqFh>){
      $numLines++;
      last if($numLines > 100000);
    }
    if(eof($fastqFh)){
      die "ERROR: $r has too few reads: ". ($numLines/4);
    }
    close $fastqFh;
  }
  logmsg "Done validating reads";

  # TODO: use validateFastq.pl?  Or is that outside of this 
  # script's scope?

}

sub makeBootstrapReads{
  my($reads,$reps,$settings)=@_;
  return [] if($reps < 1);
  
  # Enqueue the reads with a replicate ID
  my $readsQ=Thread::Queue->new();
  for(my $i=1;$i<=$reps;$i++){
    for my $r(@$reads){
      $readsQ->enqueue([$i,$r]);
    }
  }

  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&subsampleReads, $readsQ, $settings);
  }

  my @terminator=(undef) x scalar(@thr);
  $readsQ->enqueue(@terminator);
  my @reads;
  for(@thr){
    my $rList=$_->join;
    push(@reads,@$rList);
  }

  my @dir=uniq(map{dirname($_)} @reads);
  return \@dir;
}

sub subsampleReads{
  my($readsQ,$settings)=@_;
  my @fastqOut;
  while(defined(my $tmp=$readsQ->dequeue())){
    my($rep,$r)=@$tmp;
    my $readsFh=openFastq($r,$settings);

    # subsample each fastq file
    my $readCount=0;
    my $outdir="$$settings{tempdir}/subsampledReads/$rep";
    system("mkdir -p $outdir");
    my $outfile="$outdir/".basename($r,@fastqExt).".fastq";
    push(@fastqOut,$outfile);

    logmsg "(rep $rep) Getting random reads from $r";
    
    # Use CGP to get random entries via downsampling.
    # Get the first 88888 lines of the fastq file.
    # Any number ending in 888 is a multiple of 8, so it doesn't
    # split up any fastq entries.
    printRandomReadsToFile($r,$outfile,$settings);
  }

  # Return a unique list because the re-enqueuing nonsense
  # might make duplicate entries.
  # TODO: I probably don't need to uniq this because I'm not requeuing anymore.
  return [uniq(@fastqOut)];
}

# Run mash sketch on everything, multithreaded.
sub sketchAll{
  my($reads,$sketchDir,$settings)=@_;

  mkdir $sketchDir;

  my $readsQ=Thread::Queue->new(@$reads);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashSketch,$sketchDir,$readsQ,$settings);
  }
  
  $readsQ->enqueue(undef) for(@thr);

  my @mshList;
  for(@thr){
    my $mashfiles=$_->join;
    for my $file(@$mashfiles){
      push(@mshList,$file);
    }
  }

  return \@mshList;
}

# Individual mash sketch
sub mashSketch{
  my($sketchDir,$Q,$settings)=@_;

  my @msh;
  while(defined(my $fastq=$Q->dequeue)){
    logmsg "Sketching $fastq";
    my $outPrefix="$sketchDir/".basename($fastq);
    if(-e "$outPrefix.msh"){
      logmsg "WARNING: ".basename($fastq)." was already mashed. You need unique filenames for this script. This file will be skipped: $fastq";
    } elsif(-s $fastq < 1){
      logmsg "WARNING: $fastq is a zero byte file. Skipping.";
    } else {
      system("mash sketch -k 21 -s 10000 -m $$settings{mindepth} -c 10 -g $$settings{genomesize} -o $outPrefix $fastq > /dev/null 2>&1");
      die if $?;
    }

    push(@msh,"$outPrefix.msh");
  }

  return \@msh;
}

# Parallelized mash distance
sub mashDistance{
  my($mshList,$outdir,$settings)=@_;

  # Make a temporary file with one line per mash file.
  # Helps with not running into the max number of command line args.
  my $mshListFilename="$outdir/mshList.txt";
  open(my $mshListFh,">",$mshListFilename) or die "ERROR: could not write to $mshListFilename: $!";
  print $mshListFh $_."\n" for(@$mshList);
  close $mshListFh;

  my $mshQueue=Thread::Queue->new(@$mshList);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$mshQueue,$mshListFilename,$settings);
  }

  $mshQueue->enqueue(undef) for(@thr);

  my $distfile="$outdir/distances.tsv";
  open(DIST,">",$distfile) or die "ERROR: could not open $distfile for writing: $!";
  for(@thr){
    my $distfiles=$_->join;
    for my $file(@$distfiles){
      # Print the contents of each dist file to the
      # main dist file.
      open(ONEDISTFILE,"<",$file) or die "ERROR: could not open $file for reading: $!";
      while(<ONEDISTFILE>){
        print DIST $_;
      }
      close ONEDISTFILE;
    }
  }
  close DIST;

  return $distfile;
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshQueue,$mshList,$settings)=@_;
  my @dist;
  while(defined(my $msh=$mshQueue->dequeue)){
    my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh -l $mshList > $outfile");
    die if $?;

    push(@dist,$outfile);
  }

  return \@dist;
}

# 1. Read the mash distances
# 2. Create a phylip file
sub distancesToPhylip{
  my($distances,$outdir,$settings)=@_;

  my $phylip = "$outdir/distances.phylip"; 
  return $phylip if(-e $phylip);

  logmsg "Reading the distances file at $distances";
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";

  my $id="UNKNOWN"; # Default ID in case anything goes wrong
  my %m; #matrix for distances
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      $id=_truncateFilename($1,$settings);
    } else {
      my @F=split(/\t/,$_);
      $F[0]=_truncateFilename($F[0],$settings);
      $m{$id}{$F[0]}=sprintf("%0.6f",$F[1]);
    }
  }
  close MASHDIST;

  # Create the phylip file.
  # Make the text first so that we can edit it a bit.
  # TODO I should probably make the matrix the bioperl way.
  logmsg "Creating the distance matrix file for fneighbor.";
  my %seenTruncName;
  my $phylipText="";
  my @genome=sort{$a cmp $b} keys(%m);
  for(my $i=0;$i<@genome;$i++){ 
    my $name=_truncateFilename($genome[$i],$settings);
    $phylipText.="$name  "; 
    if($seenTruncName{$name}++){
      
    }
    for(my $j=0;$j<@genome;$j++){
      $phylipText.=$m{$genome[$i]}{$genome[$j]}."  ";
    }
    $phylipText.= "\n";
  }
  $phylipText=~s/  $//gm;

  # Make the phylip file.
  open(PHYLIP,">",$phylip) or die "ERROR: could not open $phylip for writing: $!";
  print PHYLIP "    ".scalar(@genome)."\n";
  print PHYLIP $phylipText;
  close PHYLIP;

  return $phylip;
}

# Create tree file with BioPerl
sub createTree{
  my($phylip,$outdir,$settings)=@_;

  logmsg "Creating a NJ tree with BioPerl";
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $matrix   = Bio::Matrix::IO->new(-format=>"phylip", -file=>$phylip)->next_matrix;
  my $treeObj = $dfactory->make_tree($matrix);
  open(TREE,">","$outdir/tree.dnd") or die "ERROR: could not open $outdir/tree.dnd: $!";
  print TREE $treeObj->as_text("newick");
  close TREE;

  return $treeObj;

}

#######
# Utils
#######

sub printRandomReadsToFile{
  my($infile,$outfile,$settings)=@_;

  my $numEntries=0;
  open(my $outFh,">",$outfile) or die "ERROR: could not open $outfile for writing: $!";
  my $fh=openFastq($infile,$settings);
  while(my $entry=<$fh> . <$fh> . <$fh> . <$fh>){
    next if(rand() < 0.5);

    print $outFh $entry;
    last if(++$numEntries > 10000);
  }
  close $fh;

  if($$settings{'save-space'}){
    system("gzip $outfile");
    die if $?;
  }
  
  return $numEntries;
}

# Removes fastq extension, removes directory name,
# truncates to a length, and adds right-padding.
sub _truncateFilename{
  my($file,$settings)=@_;
  my $name=basename($file,@fastqExt);
  $name=substr($name,0,$$settings{truncLength}); 
  $name.=" " x ($$settings{truncLength}-length($name)); 
  return $name;
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

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz > tree.dnd
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.
  --truncLength        100  How many characters to keep in a filename
  --warn-on-duplicate       Warn instead of die when a duplicate
                            genome name is found
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping. 
  --validate-reads          Do you want to see if your reads will work
                            with $0?
                            Currently checks number of reads and 
                            uniqueness of filename.
  --save-space              Save space in the temporary directory
                            where possible

  MASH SKETCH OPTIONS
  --genomesize   5000000
  --mindepth     2     
  "
}

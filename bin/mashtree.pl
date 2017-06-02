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
use File::Copy qw/copy/;
use POSIX qw/floor/;
use List::Util qw/min max/;

use threads;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Mashtree qw/logmsg @fastqExt @fastaExt @richseqExt _truncateFilename createTreeFromPhylip $MASHTREE_VERSION/;
use Mashtree::Db;
use Bio::Tree::DistanceFactory;
use Bio::Matrix::IO;
use Bio::Tree::Statistics;
use Bio::SeqIO;

my %delta :shared=(); # change in amplitude for peak detection, for each fastq
my $scriptDir=dirname $0;
my $dbhLock :shared;  # Use this as a lock so that only one thread writes to the db at a time
local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outmatrix=s tempdir=s numcpus=i genomesize=i mindepth|min-depth=i truncLength=i kmerlength=i sort-order=s sketch-size=i version)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{truncLength}||=250;  # how long a genome name is
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{'sort-order'}||="ABC";

  # Mash-specific options
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}//=5;
  $$settings{kmerlength}||=21;
  $$settings{'sketch-size'}||=10000;

  # Make some settings lowercase
  for(qw(sort-order)){
    $$settings{$_}=lc($$settings{$_});
  }

  die usage() if($$settings{help});
  if($$settings{version}){
    print "Mashtree version $MASHTREE_VERSION\n";
    return 0;
  }

  # "reads" are either fasta assemblies or fastq reads
  my @reads=@ARGV;
  die usage() if(@reads < 2);

  # Check for prereq executables.
  for my $exe(qw(mash)){
    system("$exe -h > /dev/null 2>&1");
    die "ERROR: could not find $exe in your PATH" if $?;
  }

  # Distributed cpus if we have few genomes but high numcpus
  $$settings{cpus_per_mash}=floor($$settings{numcpus}/@reads);
  $$settings{cpus_per_mash}=1 if($$settings{cpus_per_mash} < 1);
  $$settings{numthreads}=min(scalar(@reads), $$settings{numcpus});

  #die Dumper [$$settings{cpus_per_mash},$$settings{numthreads},\@reads];
  #$$settings{cpus_per_mash}=1;
  #$$settings{numthreads}=$$settings{numcpus};

  logmsg "Temporary directory will be $$settings{tempdir}";
  logmsg "$0 on ".scalar(@reads)." files";

  my %seen;
  my @tmp;
  for my $reads(@ARGV){
    my $basename=basename($reads);
    if($seen{$basename}++){
      logmsg "Skipping $reads: already seen $basename";
      next;
    }
    push(@tmp,$reads);
  }
  @reads=@tmp;

  my $sketches=sketchAll(\@reads,"$$settings{tempdir}",$settings);

  my $phylip = mashDistance($sketches,\@reads,$$settings{tempdir},$settings);

  logmsg "Creating a NJ tree with BioPerl";
  my $treeObj = createTreeFromPhylip($phylip,$$settings{tempdir},$settings);

  print $treeObj->as_text('newick');
  print "\n";

  return 0;
}

# Run mash sketch on everything, multithreaded.
sub sketchAll{
  my($reads,$sketchDir,$settings)=@_;

  mkdir $sketchDir if(!-d $sketchDir);

  # Make an array of genomes that would distribute well
  # across threads.  For example, don't put all raw-read
  # genomes into a single thread and all the assemblies
  # into another.
  my %filesize=();
  for(@$reads){
    $filesize{$_} = -s $_;
  }
  my @sortedReads=sort {$filesize{$a} <=> $filesize{$b}} @$reads;
  my @threadArr=();
  for(my $i=0; $i<@sortedReads; $i++){
    # Since each genome is sorted smallest to leargest,
    # they can be sent round-robin to each thread to 
    # ensure balance.
    my $threadIndex = $i % $$settings{numcpus};
    push(@{ $threadArr[$threadIndex] }, $sortedReads[$i]);
  }

  # Initiate the threads
  my @thr;
  for(0..$$settings{numthreads}-1){
    $thr[$_]=threads->new(\&mashSketch,$sketchDir,$threadArr[$_],$settings);
  }

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
  my($sketchDir,$genomeArr,$settings)=@_;

  # If any file needs to be converted, it will end up in
  # this directory.
  my $tempdir=tempdir("$$settings{tempdir}/convertFasta.XXXXXX", CLEANUP=>1);

  my @msh;
  # $fastq is a misnomer: it could be any kind of accepted sequence file
  for my $fastq(@$genomeArr){
    my($fileName,$filePath,$fileExt)=fileparse($fastq,@fastqExt,@fastaExt,@richseqExt);

    # Unzip the file. This temporary file will
    # only exist if the correct extensions are detected.
    my $unzipped="$tempdir/".basename($fastq);
    $unzipped=~s/\.(gz|bz2?|zip)$//i;
    if($fastq=~/\.gz$/i){
      system("gzip  -cd $fastq > $unzipped");
      die "ERROR with gzip  -cd $fastq" if $?;
    } elsif($fastq=~/\.bz2?$/i){
      system("bzip2 -cd $fastq > $unzipped");
      die "ERROR with bzip2 -cd $fastq" if $?;
    } elsif($fastq=~/\.zip$/i){
      system("unzip -p  $fastq > $unzipped");
      die "ERROR with unzip -p  $fastq" if $?;
    }

    # If the fastq filename is not the same as the unzipped
    # filename, then presumably it is unzipped and the 
    # filename needs to be updated.
    if($fastq ne $unzipped){
      $fastq=$unzipped;
      ($fileName,$filePath,$fileExt)=fileparse($fastq,@fastqExt,@fastaExt,@richseqExt);
    }

    # If we see a richseq (e.g., gbk or embl), then convert it to fasta
    # TODO If Mash itself accepts richseq, then consider
    # doing away with this section.
    if(grep {$_ eq $fileExt} @richseqExt){
      # Make a temporary fasta file, but it needs to have a
      # consistent name in case Mashtree is being run with
      # the wrapper for bootstrap values.
      # I can't exactly make a consistent filename in case
      # different mashtree invocations collide, so
      # I need to make a new temporary directory with a 
      # consistent filename.
      my $tmpfasta="$tempdir/$fileName$fileExt.fasta";
      my $in=Bio::SeqIO->new(-file=>$fastq);
      my $out=Bio::SeqIO->new(-file=>">$tmpfasta", -format=>"fasta");
      while(my $seq=$in->next_seq){
        $out->write_seq($seq);
      }
      logmsg "Wrote $tmpfasta";

      # Update our filename for downstream
      $fastq=$tmpfasta;
      ($fileName,$filePath,$fileExt)=fileparse($tmpfasta, @fastaExt);
    }

    # Do different things depending on fastq vs fasta
    my $sketchXopts="";
    if(grep {$_ eq $fileExt} @fastqExt){
      my $minDepth=determineMinimumDepth($fastq,$$settings{mindepth},$$settings{kmerlength},$settings);
      $sketchXopts.="-m $minDepth -g $$settings{genomesize} ";
    } elsif(grep {$_ eq $fileExt} @fastaExt) {
      $sketchXopts.=" ";
    } else {
      logmsg "WARNING: I could not understand what kind of file this is by its extension ($fileExt): $fastq";
    }
      
    my $outPrefix="$sketchDir/".basename($fastq);

    # See if the user already mashed this file locally
    if(-e "$fastq.msh"){
      logmsg "Found locally mashed file $fastq.msh. I will use it.";
      copy("$fastq.msh","$outPrefix.msh");
    }

    if(-e "$outPrefix.msh"){
      logmsg "WARNING: ".basename($fastq)." was already mashed.";
    } elsif(-s $fastq < 1){
      logmsg "WARNING: $fastq is a zero byte file. Skipping.";
    } else {
      logmsg "Sketching $fastq";
      system("mash sketch -p $$settings{cpus_per_mash} -k $$settings{kmerlength} -s $$settings{'sketch-size'} $sketchXopts -o $outPrefix $fastq  1>&2");
      die if $?;
    }

    push(@msh,"$outPrefix.msh");
  }

  return \@msh;
}

# Parallelized mash distance
sub mashDistance{
  my($mshList,$reads,$outdir,$settings)=@_;

  # Make a temporary file with one line per mash file.
  # Helps with not running into the max number of command line args.
  my $mshListFilename="$outdir/mshList.txt";
  open(my $mshListFh,">",$mshListFilename) or die "ERROR: could not write to $mshListFilename: $!";
  print $mshListFh $_."\n" for(@$mshList);
  close $mshListFh;

  # Instatiate the database and create the table before the threads get to it
  my $mashtreeDbFilename="$outdir/distances.sqlite";
  my $mashtreeDb=Mashtree::Db->new($mashtreeDbFilename);

  # Make an array of distance files for each thread.
  # Because distance files take about the same amount
  # of time to analyze, there is no need to sort.
  my @threadArr=();
  for(my $i=0; $i<@$mshList; $i++){
    my $threadIndex = $i % $$settings{numcpus};
    push(@{ $threadArr[$threadIndex] }, $$mshList[$i]);
  }

  # Initialize the threads
  my @thr;
  for(0..$$settings{numthreads}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$threadArr[$_],$mshListFilename,$mashtreeDbFilename,$settings);
  }

  logmsg "Joining $$settings{numthreads} threads";
  for(@thr){
    my $distfiles=$_->join;
  }

  my $phylip = "$outdir/distances.phylip";
  logmsg "Converting to phylip format into $phylip";
  open(my $phylipFh, ">", $phylip) or die "ERROR: could not write to $phylip: $!";
  print $phylipFh $mashtreeDb->toString($reads,"phylip");
  close $phylipFh;

  if($$settings{outmatrix}){
    logmsg "Writing a distance matrix to $$settings{outmatrix}";
    open(my $matrixFh, ">", $$settings{outmatrix}) or die "ERROR: could not write to $$settings{outmatrix}: $!";
    print $matrixFh $mashtreeDb->toString($reads,"matrix");
    close $matrixFh;
  }
  

  return $phylip;
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshArr,$mshList,$mashtreeDbFilename,$settings)=@_;

  # One distance file for all queries in this thread
  my($distFileFh,$distFile)=tempfile("mashdistXXXXXX", SUFFIX=>".tsv", DIR=>$outdir);

  my $numQueries=0;
  my $mashtreeDb=Mashtree::Db->new($mashtreeDbFilename);
  for my $msh(@$mshArr){
    #my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh -l $mshList >> $distFile");
    die "ERROR with 'mash dist -t $msh -l $mshList'" if $?;
    $numQueries++;
  }

  # If there is anything to add to the database, lock
  # the database.  Only lock it once per thread, optimally.
  if($numQueries > 0){
    lock($dbhLock);
    $mashtreeDb->addDistances($distFile);
  }

  # I think that the thread disconnects the db when
  # this sub ends but I wanted to do it directly and
  # in a readable fashion.
  $mashtreeDb->disconnect();
}

sub determineMinimumDepth{
  my($fastq,$mindepth,$kmerlength,$settings)=@_;

  $delta{$fastq}//=100;
  my $defaultDepth=1; # if no valley is detected

  return $mindepth if($mindepth > 0);

  my $basename=basename($fastq,@fastqExt);
  
  # Run the min abundance finder to find the valleys
  my $minAbundanceTempdir="$$settings{tempdir}/$basename.minAbundance.tmp";
  mkdir $minAbundanceTempdir;
  my $histFile="$$settings{tempdir}/$basename.hist.tsv";
  my @valley=`$scriptDir/min_abundance_finder.pl --numcpus $$settings{cpus_per_mash} $fastq --kmer $kmerlength --tempdir $minAbundanceTempdir --valleys --nopeaks --hist $histFile --delta $delta{$fastq}`;
  die "ERROR with min_abundance_finder.pl on $fastq: $!" if($?);
  chomp(@valley);
  # Some cleanup of large files
  unlink $_ for(glob("$minAbundanceTempdir/*"));
  rmdir $minAbundanceTempdir;

  # If there is no valley, return a default value
  if(!defined $valley[1]){
    $delta{$fastq}=int($delta{$fastq}/2);
    if($delta{$fastq} > 10){
      logmsg "Trying again to determine a min depth with delta==$delta{$fastq} on $fastq";
      return determineMinimumDepth(@_);
    }
    logmsg "WARNING: no valleys were found! Reporting minimum kmer coverage as $defaultDepth.";
    return $defaultDepth;
  }
  
  # Discard the header but keep the first line
  my($minKmerCount, $countOfCounts)=split(/\t/,$valley[1]);
  $minKmerCount=1 if(!defined($minKmerCount) || $minKmerCount < 1);

  logmsg "Setting the min depth as $minKmerCount for $fastq (delta==$delta{$fastq})";

  return $minKmerCount;
}

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 [options] *.fastq *.fasta *.gbk > tree.dnd
  NOTE: fastq files are read as raw reads;
        fasta, gbk, and embl files are read as assemblies;
        Input files can be gzipped.
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.
  --outmatrix          ''   If specified, will write a distance matrix
                            in tab-delimited format
  --version                 Display the version and exit

  TREE OPTIONS
  --truncLength        250  How many characters to keep in a filename
  --sort-order         ABC  For neighbor-joining, the sort order can
                            make a difference. Options include:
                            ABC (alphabetical), random, input-order

  MASH SKETCH OPTIONS
  --genomesize         5000000
  --mindepth           5    If mindepth is zero, then it will be
                            chosen in a smart but slower method,
                            to discard lower-abundance kmers.
  --kmerlength         21
  --sketch-size        10000
  "
}


#!/usr/bin/env perl
# Find the minimum abundance of kmers
# Original script was in Python at 
#   https://gist.github.com/alexjironkin/4ed43412878723491240814a0d5a6ed6/223dea45d70c9136703a4afaab0178cdbfbd2042
# Original author of python script: @alexjironkin, Public Health Englad
# I wanted to increase capatibility with Perl and have a more
# standalone script instead of relying on the khmer package.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper qw/Dumper/;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use List::Util qw/max/;
use IO::Uncompress::Gunzip qw/gunzip/;

use threads;
use Thread::Queue;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help kmerlength|kmer=i kmerCounter=s delta=i gt|greaterthan=i hist|histogram=s valleys! peaks! tempdir=s numcpus=i)) or die $!;
  $$settings{kmerlength} ||=21;
  $$settings{kmerCounter}||="";
  $$settings{delta}      ||=100;
  $$settings{gt}         ||=1;
  $$settings{tempdir}    ||=tempdir(TEMPLATE=>"$0.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{numcpus}    ||=1;

  $$settings{peaks}     //=0;
  $$settings{valleys}   //=1;

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  # Find peaks and valleys, but use the histogram file if it exists
  # and if the user specified it.
  my $histogram=[];
  if($$settings{hist} && -e $$settings{hist}){
    # Read the cached histogram file
    $histogram=readHistogram($$settings{hist},$settings);
  } else {
    # Count kmers and make a histogram.  kmerToHist()
    # will optionally write the histogram file.
    my $kmercount=countKmers($fastq,$$settings{kmerlength},$settings);
    $histogram=kmerToHist($kmercount,$settings);
  }
  # Find peaks and valleys using the simplified 
  # 'delta' algorithm.
  my $peaks=findThePeaksAndValleys($histogram,$$settings{delta},$settings);

  # Configure the output
  my @outputPos=();
  push(@outputPos, @{$$peaks{valleys}}) if($$settings{valleys});
  push(@outputPos, @{$$peaks{peaks}})   if($$settings{peaks});
  @outputPos=sort {$$a[0] <=> $$b[0]} @outputPos;

  if(!@outputPos){
    logmsg "WARNING: no peaks or valleys were reported";
  }

  # Header for output table
  print join("\t",qw(kmer count))."\n";
  # Finish off the table of output.
  for my $position (@outputPos){
    print join("\t",@$position)."\n";
  }

  return 0;
}

sub readHistogram{
  my($infile,$settings)=@_;
  my @hist=(0);
  open(HIST,$infile) or die "ERROR: could not read $infile: $!";
  logmsg "Reading histogram from $infile";
  while(<HIST>){
    chomp;
    my($count,$countOfCounts)=split /\t/;
    $hist[$count]=$countOfCounts;
  }
  close HIST;

  # ensure defined values
  for(my $i=0;$i<@hist;$i++){
    $hist[$i] //= 0;
  }

  return \@hist;
}

# TODO: count kmers with faster programs in this order of
# priority: jellyfish, KAnalyze, Khmer
# and lastly, pure perl.

sub countKmers{
  my($fastq,$kmerlength,$settings)=@_;

  my $kmerHash={};

  if($$settings{kmerCounter} =~ /^\s*$/){
    logmsg "Auto-detecting which kmer counter program to use.";
    if(which("jellyfish")){
      $$settings{kmerCounter}="jellyfish";
    } else {
      $$settings{kmerCounter}="pureperl";
    }
  }

  # try/catch kmer counting and if it fails, do the
  # pure perl method.  Don't redo pure perl if
  # it was already set up that way.
  eval{
    if($$settings{kmerCounter}=~ /(pure)?.*perl/i){
      logmsg "Counting with pure perl";
      $kmerHash=countKmersPurePerl($fastq,$kmerlength,$settings);
    } elsif($$settings{kmerCounter} =~ /jellyfish/i){
      logmsg "Counting with Jellyfish";
      $kmerHash=countKmersJellyfish($fastq,$kmerlength,$settings);
    } else {
      die "ERROR: I do not understand the kmer counter $$settings{kmerCounter}";
    }
  };

  if($@){
    if($$settings{kmerCounter}=~ /(pure)?.*perl/i){
      die "ERROR counting kmers with pure perl";
    }
    logmsg "Error detected.  Trying again by counting with pure perl";
    $kmerHash=countKmersPurePerl($fastq,$kmerlength,$settings);
  }

  return $kmerHash;
}

sub countKmersPurePerl{
  my($fastq,$kmerlength,$settings)=@_;

  # Multithreading
  my $seqQ=Thread::Queue->new;
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&countKmersPurePerlWorker,$kmerlength,$seqQ,$settings);
  }

  # Pure perl to make this standalone... the only reason
  # we are counting kmers in Perl instead of C.
  my $fastqFh=openFastq($fastq,$settings);
  my $i=0;
  my @buffer=();
  while(<$fastqFh>){ # burn the read ID line
    $i++;
    my $seq=<$fastqFh>;
    push(@buffer, $seq);

    if($i % 1000000 == 0){
      logmsg "Enqueuing ".scalar(@buffer)." reads for kmer counting";
      $seqQ->enqueue(@buffer);
      @buffer=();
    }
    
    # Burn the quality score lines
    <$fastqFh>;
    <$fastqFh>;
  }
  close $fastqFh;

  logmsg "Enqueuing ".scalar(@buffer)." reads for kmer counting";
  $seqQ->enqueue(@buffer);
  logmsg $seqQ->pending . " reads still pending...";

  # Send the termination signal
  $seqQ->enqueue(undef) for(@thr);

  my %kmer=();
  for(@thr){
    my $threadKmer=$_->join;
    for my $kmer(keys(%$threadKmer)){
      $kmer{$kmer}+=$$threadKmer{$kmer};
    }
  }

  return \%kmer;
}

sub countKmersPurePerlWorker{
  my($kmerlength,$seqQ,$settings)=@_; 

  my %kmer;
  while(defined(my $seq=$seqQ->dequeue)){

    my $numKmersInRead=length($seq)-$kmerlength+1;

    # Count kmers in a sliding window.
    # We must keep this loop optimized for speed.
    for(my $j=0;$j<$numKmersInRead;$j++){
      $kmer{substr($seq,$j,$kmerlength)}++;
    }

  }

  return \%kmer;
}


sub countKmersJellyfish{
  my($fastq,$kmerlength,$settings)=@_;
  my $basename=basename($fastq);
  my %kmer=();

  # Version checking
  my $jfVersion=`jellyfish --version`;
  # e.g., jellyfish 2.2.6
  if($jfVersion =~ /(jellyfish\s+)?(\d+)?/){
    my $majorVersion=$2;
    if($majorVersion < 2){
      die "ERROR: Jellyfish v2 or greater is required";
    }
  }
  
  my $outprefix="$$settings{tempdir}/$basename.mer_counts";
  my $jfDb="$$settings{tempdir}/$basename.merged.jf";
  my $kmerTsv="$$settings{tempdir}/$basename.jf.tsv";

  # Counting
  logmsg "Counting kmers in $fastq";
  my $jellyfishCountOptions="-s 10000000 -m $kmerlength -o $outprefix -t $$settings{numcpus}";
  my $uncompressedFastq="$$settings{tempdir}/$basename.fastq";
  if($fastq=~/\.gz$/i){
    logmsg "Decompressing fastq for jellyfish into $uncompressedFastq";
    system("zcat $fastq > $uncompressedFastq"); die if $?;
    system("jellyfish count $jellyfishCountOptions $uncompressedFastq");
  } else {
    system("jellyfish count $jellyfishCountOptions $fastq");
  }
  die "Error: problem with jellyfish" if $?;

  logmsg "Jellyfish dump $outprefix";
  my $lowerCount=$$settings{gt}+1;
  system("jellyfish dump --lower-count=$lowerCount --column --tab -o $kmerTsv $outprefix");
  die if $?;

  # Load kmers to memory
  logmsg "Reading jellyfish kmers to memory";
  open(TSV,$kmerTsv) or die "ERROR: Could not open $kmerTsv: $!";
  while(<TSV>){
    chomp;
    my @F=split /\t/;
    $kmer{$F[0]}=$F[1];
  }
  close TSV;

  # cleanup
  for($jfDb, $kmerTsv, $outprefix, $uncompressedFastq){
    unlink $_ if($_ && -e $_);
  }

  return \%kmer;
}

sub kmerToHist{
  my($kmercountHash,$settings)=@_;
  my %hist=();
  my @hist=(0);

  for my $kmercount(values(%$kmercountHash)){
    $hist{$kmercount}++;
  }

  # Turn this hash into an array
  for(1..max(keys(%hist))){
    $hist[$_] = $hist{$_} || 0;
  }

  if($$settings{hist}){
    logmsg "Writing histogram to $$settings{hist}";
    open(HIST, ">", $$settings{hist}) or die "ERROR: could not write histogram to $$settings{hist}: $!";
    for(my $i=0;$i<@hist;$i++){
      print HIST "$i\t$hist[$i]\n";
    }
    close HIST;
  }

  return \@hist;
}

sub findThePeaksAndValleys{
  my($hist, $delta, $settings)=@_;

  my($min,$max)=(MAXINT,MININT);
  my($minPos,$maxPos)=(0,0);
  my @maxTab=();
  my @minTab=();

  my $lookForMax=1;

  my $numZeros=0; # If we see too many counts of zero, then exit.
  
  for(my $kmerCount=$$settings{gt}+1;$kmerCount<@$hist;$kmerCount++){
    my $countOfCounts=$$hist[$kmerCount];
    if($countOfCounts == 0){ 
      $numZeros++;
    }
    if($countOfCounts > $max){
      $max=$countOfCounts;
      $maxPos=$kmerCount;
    }
    if($countOfCounts < $min){
      $min=$countOfCounts;
      $minPos=$kmerCount;
    }

    if($lookForMax){
      if($countOfCounts < $max - $delta){
        push(@maxTab,[$maxPos,$max]);
        $min=$countOfCounts;
        $minPos=$kmerCount;
        $lookForMax=0;
      }
    }
    else{
      if($countOfCounts > $min + $delta){
        push(@minTab,[$minPos,$min]);
        $max=$countOfCounts;
        $maxPos=$kmerCount;
        $lookForMax=1;
      }
    }

    last if($numZeros > 3);
  }

  return {peaks=>\@maxTab, valleys=>\@minTab};
}

sub findTheValley{
  my($peak1,$peak2,$hist,$settings)=@_;

  my $valley=$$peak1[1];
  my $kmerCount=$$peak1[0];
  for(my $i=$$peak1[0]+1;$i<=$$peak2[0];$i++){
    if($valley < $$hist[$i]){
      $valley=$$hist[$i];
      $kmerCount=$i;
    } else {
      last;
    }
  }
  
  return [$kmerCount,$valley];
}

# Opens a fastq file in a smart way
sub openFastq{
  my($fastq,$settings)=@_;

  my $fh;

  my @fastqExt=qw(.fastq.gz .fastq .fq.gz .fq);
  my($name,$dir,$ext)=fileparse($fastq,@fastqExt);

  # Open the file in different ways, depending on if it
  # is gzipped or if the user has gzip installed.
  if($ext =~/\.gz$/){
    # use binary gzip if we can... why not take advantage
    # of the compiled binary's speedup?
    if(-e "/usr/bin/gzip"){
      open($fh,"gzip -cd $fastq | ") or die "ERROR: could not open $fastq for reading!: $!";
    }else{
      $fh=new IO::Uncompress::Gunzip($fastq) or die "ERROR: could not read $fastq: $!";
    }
  } else {
    open($fh,"<",$fastq) or die "ERROR: could not open $fastq for reading!: $!";
  }
  return $fh;
}

# http://www.perlmonks.org/?node_id=761662
sub which{
  my($exe,$settings)=@_;
  
  my $tool_path="";
  for my $path ( split /:/, $ENV{PATH} ) {
      if ( -f "$path/$exe" && -x "$path/$exe" ) {
          $tool_path = "$path/$exe";
          last;
      }
  }
  
  return $tool_path;
}

sub usage{
  "
  $0: 
       Find the valley between two peaks on a set of kmers
       such that you can discard the kmers that are 
       likely representative of contamination.
       This script does not require any dependencies.

  Usage: $0 file.fastq[.gz]
  --gt     1   Look for the first peak at this kmer count
               and then the next valley.
  --kmer   21  kmer length
  --delta  100 How different the counts have to be to
               detect a valley or peak
  --numcpus  1

  OUTPUT
  --hist   ''  A file to write the histogram, or a file
               to read the histogram if it already exists.
               Useful if you want to rerun this script.
  --valleys,   To output valleys (default)
  --novalleys  To supppress valleys
  --peaks,     To output maximum peaks (default)
  --nopeaks    To suppress maximum peaks

  MISC
  --kmerCounter ''  The kmer counting program to use.
                    Default: (empty string) auto-choose
                    Options: perl, jellyfish
  "
}


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
use List::Util qw/max/;
use IO::Uncompress::Gunzip qw/gunzip/;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help kmerlength|kmer=i delta=i gt|greaterthan=i hist|histogram=s valleys! peaks!)) or die $!;
  $$settings{kmerlength}||=21;
  $$settings{delta}     ||=100;
  $$settings{gt}        ||=0;

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

sub countKmers{
  my($fastq,$kmerlength,$settings)=@_;
  my %kmer=();

  # Pure perl to make this standalone... the only reason
  # we are counting kmers in Perl instead of C.
  my $fastqFh=openFastq($fastq,$settings);
  my $i=0;
  while(<$fastqFh>){ # burn the read ID line
    my $seq=<$fastqFh>;
    chomp($seq);
    my $numKmersInRead=length($seq)-$kmerlength+1;

    # Count kmers in a sliding window.
    # We must keep this loop optimized for speed.
    for(my $j=0;$j<$numKmersInRead;$j++){
      $kmer{substr($seq,$j,$kmerlength)}++;
    }

    # Burn the quality score lines
    <$fastqFh>;
    <$fastqFh>;

  }
  close $fastqFh;

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


sub usage{
  "
  $0: 
       Find the valley between two peaks on a set of kmers
       such that you can discard the kmers that are 
       likely representative of contamination.
       This script does not require any dependencies.

  Usage: $0 file.fastq[.gz]
  --gt     0   Look for the first peak at this kmer count
               and then the next valley.
  --kmer   21  kmer length
  --delta  100 How different the counts have to be to
               detect a valley or peak

  OUTPUT
  --hist   ''  A file to write the histogram, or a file
               to read the histogram if it already exists.
               Useful if you want to rerun this script.
  --valleys,   Valleys will be in the output by default
  --novalleys
  --peaks,     Maximum peaks will not be in the output
  --nopeaks    by default
  "
}


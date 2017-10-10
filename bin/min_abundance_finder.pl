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

use Bio::Kmer;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help kmerlength|kmer=i kmerCounter=s delta=i gt|greaterthan=i valleys! peaks! tempdir=s numcpus=i)) or die $!;
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

  my $kmer=Bio::Kmer->new($fastq,{
      numcpus       => $$settings{numcpus},
      kmercounter   => "perl",
  });
  my $histogram=$kmer->histogram();

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


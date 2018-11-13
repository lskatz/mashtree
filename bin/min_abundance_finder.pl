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
#use Bio::Kmer;

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
  GetOptions($settings,qw(help kmerlength|kmer=i kmerCounter=s delta=i gt|greaterthan=i tempdir=s numcpus=i)) or die $!;
  $$settings{kmerlength} ||=21;
  $$settings{kmerCounter}||="";
  $$settings{delta}      ||=10;
  $$settings{gt}         ||=1;
  $$settings{tempdir}    ||=tempdir(TEMPLATE=>"$0.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{numcpus}    ||=1;

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  my %firstValleyVote;
  for(my $i=9; $i<=23; $i+=2){
    $$settings{kmerlength}=$i;
    my $histogram=mashHistogram($fastq,$settings);
    my $tmp=findThePeaksAndValleys($histogram,$$settings{delta},$settings);
    my $firstValley=$$tmp{valleys}[0][0] || next;
    $firstValleyVote{$firstValley}++;
  }


  logmsg "For various values of k, valleys were found:";
  my $firstValley=0;
  # Sort bins by their votes, highest to lowest
  for my $bin(sort{$firstValleyVote{$b}<=>$firstValleyVote{$a}} keys(%firstValleyVote)){
    my $value=$firstValleyVote{$bin};
    $firstValley||=$bin; # set the valley to the first bin we come to
    logmsg "  $bin: $value votes";
  }

  print join("\t",qw(kmer count))."\n";
  print join("\t", $firstValley, 1)."\n";
  return 0;

}

# Poor man's way of subsampling
# Thanks to Nick Greenfield for pointing this out.
sub mashHistogram{
  my($fastq,$settings)=@_;
  my $sketch="$$settings{tempdir}/sketch.msh";
  system("mash sketch -k $$settings{kmerlength} -m 2 -o $sketch $fastq > /dev/null 2>&1");
  die if $?;
  
  my @histogram;
  open(my $fh, "mash info -c $sketch | ") or die "ERROR: could not get mash info on sketch $sketch";
  while(my $line=<$fh>){
    $line=~s/^\s+|\s+$//g;  # whitespace trim
    next if($line=~/^#/);
    my($filename, $bin, $frequency)=split(/\t/, $line);
    $histogram[$bin]=$frequency;
  }
  close $fh;
  $histogram[$_]||=0 for(0..@histogram);
  return \@histogram;
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
  --gt     1   Look for the first peak at this kmer count
               and then the next valley.
  --kmer   21  kmer length
  --delta  100 How different the counts have to be to
               detect a valley or peak
  --numcpus  1 (not currently used)

  MISC
  --kmerCounter ''  The kmer counting program to use.
                    Default: (empty string) auto-choose
                    Options: perl, jellyfish
  "
}


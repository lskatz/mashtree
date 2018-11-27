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

  # Find valleys with multithreading
  my @thr;
  for(my $kmerlength=5; $kmerlength<=32; $kmerlength+=3){
    logmsg "Counting $kmerlength-kmers in $fastq";
    my $delta = $$settings{delta}; # make a copy for threads
    push(@thr,
      threads->new(sub{
        my($fastq, $kmerlength, $delta)=@_;
          #my $kmerCounter = Bio::Kmer->new($fastq,{numcpus=>$$settings{numcpus},kmerlength=>$kmerlength,sample=>0.01});
          #my $histogram = $kmerCounter->histogram();
          my $histogram   = mashHistogram($fastq,$kmerlength,$settings);
          my $firstValley = findFirstValley($histogram, $settings);
          return $firstValley;

        },$fastq, $kmerlength, $delta
      )
    );
  }

  my %firstValleyVote;
  for(@thr){
    my $firstValley = $_->join;
    next if(!$firstValley);
    $firstValleyVote{$firstValley}++;
  }

  # If no valleys were found, simply set the stage so that
  # the minimum depth will be set to 0.
  if(keys(%firstValleyVote) < 1){
    logmsg "NOTE: no valleys were found and so I am inserting an imaginary vote for a valley at cov=1";
    $firstValleyVote{0} = 1;
  }

  logmsg "For various values of k, valleys were found:";
  my $firstValley=0;
  # Average out the votes
  my $totalFirstValley = 0;
  my $totalVotes = 0;
  my @vote;
  # Sort bins by their votes, highest to lowest
  for my $bin(sort{$firstValleyVote{$b}<=>$firstValleyVote{$a} || $a<=>$b} keys(%firstValleyVote)){
    my $value=$firstValleyVote{$bin};
    $firstValley||=$bin; # set the valley to the first bin we come to
    for(1..$value){
      push(@vote, $bin);
    }
    logmsg "  $bin: $value votes";

    $totalFirstValley += $bin * $value;
    $totalVotes += $value;
  }
  @vote = sort {$a<=>$b} @vote;
  my $medianFirstValley = $vote[ int(scalar(@vote)/2) ];

  # Get the average first valley across many kmers
  my $avgFirstValley = $totalFirstValley/$totalVotes;
  logmsg "    Average first valley is $avgFirstValley";
  logmsg "    However, I will use the median valley: $medianFirstValley";
  printf("%0.0f\n", $medianFirstValley);

  #print join("\t",qw(kmer count))."\n";
  #print join("\t", $firstValley, 1)."\n";
  return 0;

}

# Poor man's way of subsampling
# Thanks to Nick Greenfield for pointing this out.
sub mashHistogram{
  my($fastq,$k,$settings)=@_;
  my $sketch="$$settings{tempdir}/sketch.msh";
  system("mash sketch -k $k -m 2 -o $sketch $fastq > /dev/null 2>&1");
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


# https://www.perlmonks.org/?node_id=629742
sub localMinimaMaxima{
  my($array, $settings)=@_;

  my @minima;
  my @maxima;
  my $prev_cmp = 0;

  my $num = @$array - 2;
  for my $i (0 .. $num){
    my $cmp = $$array[$i] <=> $$array[$i+1];
    if ($cmp != $prev_cmp) {
      if($cmp < 0){
        push @minima, $i;
      }
      elsif($cmp > 0){
        push @maxima, $i;
      }
      # when this and next elements are ==, defer checking for
      # minima/maxima till next loop iteration
      $prev_cmp = $cmp if $cmp;
    }
  }

  # Uncomment the following if we want to look at the very
  # last number in the array.
  #if (@$array) {
  #  push @minima, $num if $prev_cmp >= 0;
  #  push @maxima, $num if $prev_cmp <= 0;
  #}

  return(\@minima, \@maxima);
}

sub findFirstValley{
  my($array, $settings)=@_;
  my($minima, $maxima) = localMinimaMaxima($array, $settings);

  # Return the first minimum if it's not the first element
  # or if there are no other minima.
  if($$minima[0] > 0 || !$$minima[1]){
    return $$minima[0];
  } else {
    return $$minima[1];
  }
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
  --numcpus  1 (not currently used)

  MISC
  --kmerCounter ''  The kmer counting program to use.
                    Default: (empty string) auto-choose
                    Options: perl, jellyfish
  "
}


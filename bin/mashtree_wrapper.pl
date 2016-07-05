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

use threads;
use Thread::Queue;
use threads::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Mashtree qw/logmsg @fastqExt @fastaExt/;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Tree::DistanceFactory;
use Bio::Tree::Statistics;
use Bio::Matrix::IO;

local $0=basename $0;
my $fhStick :shared;  # A thread can only open a fastq file if it has the talking stick.
my $writeStick :shared;  # Only one thread can write at a time

exit main();

sub main{
  my $settings={};
  my @wrapperOptions=qw(help outmatrix=s distance-matrix tempdir=s numcpus=i save-space reps=i);
  my @mashOptions=qw(kmerlength=i truncLength=i genomesize=i mindepth=i warn-on-duplicate);
  GetOptions($settings,@wrapperOptions,@mashOptions) or die $!;
  $$settings{numcpus}||=1;
  $$settings{reps}||=0;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});

  logmsg "Temporary directory will be $$settings{tempdir}";

  die usage() if($$settings{help});

  # "reads" are either fasta assemblies or fastq reads
  my @reads=@ARGV;
  die usage() if(@reads < 2);
  
  # Figure out options for mashtree.pl
  my $mashtreexopts="--numcpus $$settings{numcpus} ";
  for my $getopt(@mashOptions){
    my $option=$getopt;
    $option=~s/=.+$//; # remove equals sign and after

    # Treat boolean and options with values differently.
    # If the option never had equals sign in it, it's boolean
    # and the value should be stripped.
    my $value=$$settings{$option};
    if($getopt eq $option){
      $value="";
    }
    $mashtreexopts.="--$option $value " if($$settings{$option});
  }

  system("$FindBin::RealBin/mashtree.pl --tempdir $$settings{tempdir} $mashtreexopts @reads > $$settings{tempdir}/tree.dnd");
  die if $?;
  my $guideTree=Bio::TreeIO->new(-file=>"$$settings{tempdir}/tree.dnd")->next_tree;

  my @bsTree;
  for(my $i=1;$i<=$$settings{reps}; $i++){
    my $tempdir="$$settings{tempdir}/rep$i";
    mkdir($tempdir);
    logmsg "Made replicate directory $tempdir";
    next if(-e "$tempdir/tree.dnd");

    # 1. Downsample/subsample
    my $reads=subsampleAll(\@reads,$tempdir,$settings);
    # 2. Run Mashtree on Subsampled
    system("$FindBin::RealBin/mashtree.pl --tempdir $tempdir $mashtreexopts @$reads > $tempdir/tree.dnd.tmp && mv $tempdir/tree.dnd.tmp $tempdir/tree.dnd");
    die if $?;

    push(@bsTree,Bio::TreeIO->new(-file=>"$tempdir/tree.dnd")->next_tree);
  }
  
  # Combine trees into a bootstrapped tree and write it 
  # to an output file. Then print it to stdout.
  my $biostat=Bio::Tree::Statistics->new;
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

  if($$settings{'distance-matrix'}){
    cp("$$settings{tempdir}/distances.tsv",$$settings{'distance-matrix'});
  }
  
  return 0;
}

#######
# Utils
#######

sub subsampleAll{
  my($reads,$tempdir,$settings)=@_;
  
  my $readsQ=Thread::Queue->new(@$reads);
  $readsQ->enqueue((undef) x $$settings{numcpus});
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&subsample, $tempdir, $readsQ, $settings);
  }

  my @outfile;
  for(@thr){
    my $outfile=$_->join;
    push(@outfile,@$outfile);
  }
  return \@outfile;
}

sub subsample{
  my($tempdir,$Q,$settings)=@_;
  # Subsample if it's an assembly.
  # Downsample if it's reads.
  my @outfile;
  while(defined(my $read=$Q->dequeue)){
    my($name,$dir,$ext)=fileparse($read,@fastqExt,@fastaExt);
    my $outfile;
    if(grep(/\Q$ext\E/,@fastqExt)){
      $outfile="$tempdir/$name.fastq";
      printRandomReadsToFile($read,$outfile,$settings);
    } elsif (grep(/\Q$ext\E/,@fastaExt)){
      $outfile="$tempdir/$name.fasta";
      subsampleAssembly($read,$outfile,$settings);
    } else {
      die "ERROR: I do not understand extension $ext";
    }
    push(@outfile,$outfile);
  }
  return \@outfile;
}

sub subsampleAssembly{
  my($infile,$outfile,$settings)=@_;

  my $seqin=Bio::SeqIO->new(-file=>$infile);
  my $seqout=Bio::SeqIO->new(-file=>">$outfile");
  while(my $seq=$seqin->next_seq){

    # Sample 50% of the contig
    my $randStart=int(rand($seq->length));
    my $sublength=$seq->length/2; 
    my $subsequence=substr($seq->seq,$randStart-1,$sublength);
    if($sublength > length($subsequence)){
      $subsequence.=substr($seq->seq,0,($sublength - length($subsequence)));
    }

    my $subseq=Bio::Seq->new(-id=>$seq->id,-seq=>$subsequence);
    $seqout->write_seq($subseq);
  }
}

sub printRandomReadsToFile{
  my($infile,$outfile,$settings)=@_;

  logmsg "Downsampling $infile";

  return 0 if(-e $outfile);

  my $numEntries=0;
  my $fastqCache="";
  open(my $outFh,">","$outfile.tmp") or die "ERROR: could not open $outfile.tmp for writing: $!";
  my $fh=openFastq($infile,$settings);
  while(my $entry=<$fh> . <$fh> . <$fh> . <$fh>){
    next if(rand() < 0.5);

    $fastqCache.=$entry;
    ++$numEntries;

    if($numEntries % 1000000 == 0){
      lock($writeStick);
      print $outFh $fastqCache;
      $fastqCache="";
      logmsg "Wrote $numEntries records to $outfile";
    }
  }
  print $outFh $fastqCache; # clear out the rest of the cache
  logmsg "Finished writing $numEntries records to $outfile";
  close $outFh;
  close $fh;

  if($$settings{'save-space'}){
    die "TODO";
    system("gzip $outfile.tmp && mv $outfile.tmp.gz $outfile.gz");
    die if $?;
  } else {
    mv("$outfile.tmp",$outfile);
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
  Usage: $0 *.fastq.gz *.fasta > tree.dnd
  NOTE: fasta files are read as assembly files; fastq files
        are read as raw reads. Fastq file can be gzipped.
  MASH TREE OPTIONS
  --tempdir                 If not specified, one will be made for you
                            and then deleted at the end of this script.
  --numcpus            1    This script uses Perl threads.
  --truncLength        250  How many characters to keep in a filename
  --warn-on-duplicate       Warn instead of die when a duplicate
                            genome name is found
  MISC OPTIONS
  --distance-matrix    ''   Output file for distance matrix
  --reps               0    How many bootstrap repetitions to run;
                            If zero, no bootstrapping.
  --validate-reads          Do you want to see if your reads will work
                            with $0?
                            Currently checks number of reads and
                            uniqueness of filename.
  --save-space              Save space in the temporary directory
                            where possible
  MASH SKETCH OPTIONS
  --genomesize         5000000
  --mindepth           5     
  --kmerlength         21
  "
}

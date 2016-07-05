package Mashtree;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use threads;

use lib dirname($INC{"Mashtree.pm"})."/lib/perl5";

our @EXPORT_OK = qw(
           logmsg 
           @fastqExt @fastaExt @bamExt @vcfExt @richseqExt
         );

local $0=basename $0;

######
# CONSTANTS

our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fas .fa);
our @bamExt=qw(.sorted.bam .bam);
our @vcfExt=qw(.vcf.gz .vcf);
our @richseqExt=qw(.gbk .gbf .gb .embl);

#################################################
### COMMON SUBS/TOOLS (not object subroutines) ##
#################################################
# Redefine how scripts die
$SIG{'__DIE__'} = sub {
  local $0=basename($0);
  my $e = $_[0] || ""; 
  my $callerSub=(caller(1))[3] || (caller(0))[3] || "UnknownSub";

  $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; 
  die("$0: $callerSub: $e"); 
};
# Centralized logmsg
#sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}
sub logmsg {
  local $0=basename $0;
  my $parentSub=(caller(1))[3] || (caller(0))[3];
  $parentSub=~s/^main:://;

  # Find the thread ID and stringify it
  my $tid=threads->tid;
  $tid=($tid) ? "(TID$tid)" : "";

  my $msg="$0: $parentSub$tid: @_\n";

  print STDERR $msg;
}


1; # gotta love how we we return 1 in modules. TRUTH!!!


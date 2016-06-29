# mashtree
Create a tree using Mash distances

## Examples

    mashtree.pl *.fastq.gz > tree.dnd
    mashtree.pl --reps 100 --numcpus 12 --validate-reads *.fastq.gz > tree.dnd

## Usage

    mashtree.pl: use distances from Mash (min-hash algorithm) to make a NJ tree.
    Input files can be fastq or fasta.  Fastq files are assumed to be reads
    while fasta will be assumed to be assemblies.
      Usage: mashtree.pl *.fastq.gz > tree.dnd
      --tempdir                 If not specified, one will be made for you
                                and then deleted at the end of this script.
      --numcpus            1    This script uses Perl threads.
      --truncLength        250  How many characters to keep in a filename
      --warn-on-duplicate       Warn instead of die when a duplicate
                                genome name is found
      --reps               0    How many bootstrap repetitions to run;
                                If zero, no bootstrapping.
      --validate-reads          Do you want to see if your reads will work
                                with mashtree.pl?
                                Currently checks number of reads and
                                uniqueness of filename.
      --save-space              Save space in the temporary directory
                                where possible

      MASH SKETCH OPTIONS
      --genomesize   5000000
      --mindepth     2

## Requirements

* Mash >= v1.1
* Perl with multithreading

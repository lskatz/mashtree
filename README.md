# mashtree
Create a tree using Mash distances.

For simple usage, see `mashtree.pl`.  For advanced options, look at `mashtree_wrapper.pl`.

## Examples

    mashtree.pl --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

**Note**: fastq files are interpreted as raw read files. Fasta,
GenBank, and EMBL files are interpreted as genome
assemblies.

**Note**: Compressed files are also accepted of any of the
above file types.  You can compress with gz, bz2, or zip.

### Advanced

    mashtree_wrapper.pl --reps 100 -- --numcpus 12 *.fastq.gz > mashtree.dnd


## Usage

    mashtree.pl: use distances from Mash (min-hash algorithm) to make a NJ tree
      Usage: mashtree.pl [options] *.fastq *.fasta *.gbk > tree.dnd
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


Also see `mashtree_wrapper.pl` for advanced usage. Run either script with
`--help` for additional information.

## Requirements

* Mash >= v1.1
* SQLite3
* Perl 
  * multithreading 
  * BioPerl library
  * `DBD::SQLite`

## Installation

After downloading the latest release, go into the directory and run `make`

    $ cd mashtree
    $ perl Makefile.PL 
    $ make test

Add `mashtree/bin` to `PATH` and you're good to go!

## References

*  Mash: http://mash.readthedocs.io
*  BioPerl: http://bioperl.org

## GitHub stickers for Mashtree

[![Build Status](https://travis-ci.org/lskatz/mashtree.svg?branch=master)](https://travis-ci.org/lskatz/mashtree)

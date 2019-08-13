# mashtree

Create a tree using Mash distances.

For simple usage, see `mashtree`.  For advanced options, look at `mashtree_wrapper.pl`.

## Two modes: fast or accurate

**Input files**: fastq files are interpreted as raw read files. Fasta,
GenBank, and EMBL files are interpreted as genome
assemblies. Compressed files are also accepted of any of the
above file types.  You can compress with gz, bz2, or zip.

**Output files**: Newick (.dnd).  If `--outmatrix` is supplied, then
a distance matrix too.

### Faster

    mashtree --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

### More accurate

You can get a more accurate tree with the minimum abundance finder. Simply
give `--mindepth 0`.  This step helps ignore very unique kmers that are 
more likely read errors.

    mashtree --mindepth 0 --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

### Adding confidence values

Mashtree can add confidence values using jack knifing. For each
jack knife tree, 50% of hashes are used. Confidence values are calculated from
the jack knife trees using BioPerl. When using this method, you can pass
flags to `mashtree` using the double-dash like in the example below.

Added in version 0.40.

    mashtree_jackknife.pl --reps 100 --numcpus 12 *.fastq.gz -- --min-depth 0 > mashtree.jackknife.dnd
    mashtree_jackknife.pl --help # additional usage help

Bootsrapping was added in version 0.55.  This runs mashtree itself multiple times, each
with a random seed.

    mashtree_bootstrap.pl --reps 100 --numcpus 12 *.fastq.gz -- --min-depth 0 > mashtree.bootstrap.dnd

## Usage

    Usage: mashtree [options] *.fastq *.fasta *.gbk *.msh > tree.dnd
    NOTE: fastq files are read as raw reads;
          fasta, gbk, and embl files are read as assemblies;
          Input files can be gzipped.
    --tempdir            ''   If specified, this directory will not be
                              removed at the end of the script and can
                              be used to cache results for future
                              analyses.
                              If not specified, a dir will be made for you
                              and then deleted at the end of this script.
    --numcpus            1    This script uses Perl threads.
    --outmatrix          ''   If specified, will write a distance matrix
                              in tab-delimited format
    --file-of-files           If specified, mashtree will try to read
                              filenames from each input file. The file of
                              files format is one filename per line. This
                              file of files cannot be compressed.
    --outtree                 If specified, the tree will be written to
                              this file and not to stdout. Log messages
                              will still go to stderr.
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

* Mash >= v2.0
* SQLite3
* Perl 
  * multithreading 
  * BioPerl library
  * `DBD::SQLite`
* Quicktree

## Installation

### Installation from Git

After downloading the latest release, go into the directory and run `make`

    $ cd mashtree
    $ perl Makefile.PL 
    $ make test

Add `mashtree/bin` to `PATH` and you're good to go!

### Installation from CPAN

Installing from CPAN installs the latest stable version of Mashtree.  This method _should_ add the Mashtree perl modules to the correct place in your home directory and _should_ add the executables to your local bin directory.  However, I am new to CPAN, so please give me feedback via the issues tab if this is not correct.

    $ cpanm -L ~ Mashtree
    $ export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
    $ mashtree --help # verify it shows usage and not an error

### Alternate method of installing from CPAN

    $ cpan # initiates the CPAN command line prompt
    cpan[1]> install Mashtree
    cpan[2]> exit
    $ export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
    $ mashtree --help # verify it shows usage and not an error

### Uninstallation from CPAN

I'm not sure _why_ you'd want to uninstall Mashtree but here is how you would clean it up.

    $ cpanm --uninstall Mashtree --local-lib=$HOME

## Further documentation

For more information and help please see the [docs folder](docs/)

For more information on plugins, see the [plugins folder](plugins). (in development)

## References

*  Mash: http://mash.readthedocs.io
*  BioPerl: http://bioperl.org

## Citation

The paper is in preparation but for now, this is a valid citation:

Katz, L. S., Griswold, T., & Carleton, H. A. (2017, October 8-11). [_Generating WGS Trees with Mashtree_](misc/mashtree%20ASM%20NGS.pptx). Poster presented at the American Society for Microbiology Conference on Rapid Applied Microbial Next-Generation Sequencing and Bioinformatic Pipelines, Washington, DC. 

## GitHub stickers for Mashtree

[![Build Status](https://travis-ci.org/lskatz/mashtree.svg?branch=master)](https://travis-ci.org/lskatz/mashtree)


# mashtree

Create a tree using Mash distances.

For simple usage, see `mashtree --help`. For confidence values, run either with `--help`: `mashtree_bootstrap.pl` or `mashtree_jackknife.pl`.

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


## Installation

**Please see [INSTALL.md](docs/INSTALL.md)**

## Further documentation

For more information and help please see the [docs folder](docs/)

For more information on plugins, see the [plugins folder](plugins). (in development)

For more information on contributions, please see [CONTRIBUTING.md](CONTRIBUTING.md).

## References

*  Mash: http://mash.readthedocs.io
*  BioPerl: http://bioperl.org

## Citation

The paper is in preparation but for now, this is a valid citation:

Katz, L. S., Griswold, T., & Carleton, H. A. (2017, October 8-11). [_Generating WGS Trees with Mashtree_](misc/mashtree%20ASM%20NGS.pptx). Poster presented at the American Society for Microbiology Conference on Rapid Applied Microbial Next-Generation Sequencing and Bioinformatic Pipelines, Washington, DC. 

## GitHub stickers for Mashtree

[![Build Status](https://travis-ci.org/lskatz/mashtree.svg?branch=master)](https://travis-ci.org/lskatz/mashtree)


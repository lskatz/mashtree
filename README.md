# mashtree
Create a tree using Mash distances. For an overview of Mash, please see http://mash.readthedocs.io

For simple usage, see `mashtree.pl`.  For advanced options, look at `mashtree_wrapper.pl`.

## Examples

    mashtree.pl --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

### Advanced

    mashtree_wrapper.pl --reps 100 --numcpus 12 --validate-reads *.fastq.gz > mashtree.dnd


## Usage

    mashtree.pl: use distances from Mash (min-hash algorithm) to make a NJ tree.
    Input files can be fastq or fasta.  Fastq files are assumed to be reads
    while fasta will be assumed to be assemblies.
      Usage: mashtree.pl *.fastq.gz > tree.dnd
      --tempdir                 If not specified, one will be made for you
                                and then deleted at the end of this script.
      --numcpus            1    This script uses Perl threads.
      --truncLength        250  How many characters to keep in a filename
      --save-space              Save space in the temporary directory
                                where possible

      MASH SKETCH OPTIONS
      --genomesize   5000000
      --mindepth     2

Also see `mashtree_wrapper.pl` for advanced usage. Run either script with
`--help` for additional information.

## Requirements

* Mash >= v1.1
* Perl with multithreading

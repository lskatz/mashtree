# mashtree
Create a tree using Mash distances.

For simple usage, see `mashtree.pl`.  For advanced options, look at `mashtree_wrapper.pl`.

## Examples

    mashtree.pl --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd

### Advanced

    mashtree_wrapper.pl --reps 100 -- --numcpus 12 *.fastq.gz > mashtree.dnd


## Usage

    mashtree.pl: use distances from Mash (min-hash algorithm) to make a NJ tree
      Usage: mashtree.pl *.fastq.gz *.fasta > tree.dnd
      NOTE: fasta files are read as assembly files; fastq files
            are read as raw reads. Fastq file can be gzipped.
      --tempdir                 If not specified, one will be made for you
                                and then deleted at the end of this script.
      --numcpus            1    This script uses Perl threads.

      TREE OPTIONS
      --truncLength        250  How many characters to keep in a filename
      --sort-order         ABC  For neighbor-joining, the sort order can
                                make a difference. Options include:
                                ABC (alphabetical), random, input-order

      MASH SKETCH OPTIONS
      --genomesize         5000000
      --mindepth           0    If mindepth is zero, then it will be
                                chosen in a smart but slower method,
                                to discard lower-abundance kmers.
      --kmerlength         21
      --sketch-size        10000


Also see `mashtree_wrapper.pl` for advanced usage. Run either script with
`--help` for additional information.

## Requirements

* Mash >= v1.1
* Perl with multithreading and with the BioPerl library

## References

*  Mash: http://mash.readthedocs.io
*  BioPerl: http://bioperl.org

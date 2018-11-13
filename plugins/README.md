# Mashtree plugins

Mashtree plugins are executed such that a mashtree database
is the first positional parameter. Flag options are allowed.
Additional tables within the database might get created, but
overwriting the DISTANCE table is not allowed.

The first script typically will be mashtree\_init, which creates
the database.

For more information on the standard database, please see
[the documentation](../docs/SQL.md).

## Example workflow

The most basic workflow is to initialize the database, run
Mash and find distances, and then print the tree.

    mashtree_init.pl lambda.sqlite
    mashtree_mash lambda.sqlite t/lambda/*.fastq.gz
    mashtree_tree lambda.sqlite

## Synopses for individual plugins

Run any plugin with `--help` for usage.

### mashtree\_init

initialize an empty mashtree database
    
    mashtree_init.pl db.sqlite

### mashtree\_mash

Run mash on genomes and insert their distances into the database

    mashtree_mash.pl [options] mash.sqlite *.fastq *.fasta *.gbk *.msh

### mashtree\_dump

Dump distances from a mashtree database.  Note: you can also run
`sqlite3 file.sqlite .dump` to view raw data.

    mashtree_dump.pl db.sqlite

### mashtree\_optimize

This plugin finds and estimates distances between all genomes. Currently
only uses Dijkstra's algorithm but potentially could be expanded. Edits
the database in-place.

    mashtree_optimize.pl mashtree.sqlite

### mashtree\_tree

Calculate and print the tree

    mashtree_tree.pl mash.sqlite

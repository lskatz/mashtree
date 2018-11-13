# Mashtree plugins

Mashtree plugins are executed such that a mashtree database
is the first positional parameter. Flag options are allowed.
Additional tables within the database might get created, but
overwriting the DISTANCE table is not allowed.

The first script typically will be mashtree\_init, which creates
the database.

For more information on the standard database, please see
[the documentation](../docs/SQL.md).

## Synopses for individual plugins

Run any plugin with `--help` for usage.

### mashtree\_init

initialize an empty mashtree database

### mashtree\_mash

Run mash on genomes and insert their distances into the database

### mashtree\_dump

Dump distances from a mashtree database.  Note: you can also run
`sqlite3 file.sqlite .dump` to view raw data.

### mashtree\_optimize

This plugin finds and estimates distances between all genomes. Currently
only uses Dijkstra's algorithm but potentially could be expanded.

    Usage: mashtree_optimize.pl [options] mashtree.sqlite

### mashtree\_tree

Calculate and print the tree


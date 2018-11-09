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

### mashtree\_optimize

This plugin optimizes distances between all genomes. Currently
only uses Dijkstra's algorithm but potentially could be expanded.

    Usage: mashtree_optimize.pl [options] mashtree.sqlite


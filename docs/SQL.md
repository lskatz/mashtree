# SQL database

Mashtree saves distances in an SQLite database in the temporary folder.
This folder is automatically deleted at the end of the run unless `--tempdir` is supplied.
This database can be read using SQLite if saved.

## SQLite3

SQLite is a local database management system. To install SQLite on most systems, one of these commands would work.

    yum install sqlite

    apt-get install sqlite

## Reading the database

The database file is `distances.sqlite`. The scheme is very simple with three fields: GENOME1, GENOME2, DISTANCE.

### Sample query

    sqlite3 ./mashtree.tmp/distances.sqlite "SELECT * FROM DISTANCE LIMIT 10"

## Writing to the database

Not recommended but what the heck

### Sample query

    sqlite3 ./mashtree.tmp/distances.sqlite "INSERT INTO DISTANCE VALUES('somegenome1','othergenome',0.01);"
    

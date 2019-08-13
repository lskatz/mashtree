# Troubleshooting

## General

Here are some general tips for when things go wrong.

### Read integrity

Your reads might be badly formatted or the number of reads might be too low.
For validating read format, the github version of Mashtree comes bundled 
with [https://github.com/lskatz/ross](ROSS) which has a read validator. You can optionally download
ROSS and run it separately.
Then, run `friends_ross.pl --help` and follow instructions on how to check
each read set.  


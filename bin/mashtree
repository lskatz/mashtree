#!/bin/bash
set -e
set -U

db=$1; shift
reads=$@
dirname=$(realpath $0);

perl $dirname/../plugins/mashtree_init.pl $db
perl $dirname/../plugins/mashtree_mash.pl $db $@
perl $dirname/../plugins/mashtree_tree.pl $db


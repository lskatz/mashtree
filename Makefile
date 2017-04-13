# Author: Lee Katz <lkatz@cdc.gov>
# Mashtree
 
PROFILE := $(HOME)/.bashrc
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

###################################

.DEFAULT: install

.PHONY: install install-perlModules

.DELETE_ON_ERROR:

install: perl-modules

perl-modules: check-perl check-bioperl
	cpanm -v -L . Bio::Tree::DistanceFactory
	cpanm -v -L . Bio::Matrix::IO
	cpanm -v -L . Bio::Tree::Statistics
	cpanm -v -L . DBD::SQLite

check: check-perl check-bioperl check-sqlite3 check-mash
	@echo "Mashtree seems to be working! Be sure to add $$(pwd)/bin to your path!"

check-perl:
	@perl -e 1 || (echo "ERROR: perl not found" && exit 1)
	@perl -Mthreads -e 1 || (echo "ERROR: this version of perl does not have threads" && exit 1)

check-bioperl: check-perl 
	@perl -MBio::Perl -e 1 || (echo "ERROR: BioPerl not found! For installation instructions: http://bioperl.org/INSTALL.html" && exit 1)
	
check-sqlite3:
	@sqlite3 --version >/dev/null || (echo "ERROR: sqlite3 not found! For installation: https://sqlite.org/download.html" && exit 1)

check-mash:
	@mash --version >/dev/null || (echo "ERROR: mash not found! For installation: http://mash.readthedocs.io" && exit 1)

test:
	# Test to see if the first line of the usage is correct, 
	# ie that the script runs at all
	if [[ $$(mashtree.pl --help 2>&1 | head -n 1  | cut -f 1 -d ' ') != "mashtree.pl:" ]]; then exit 1; fi;
	@echo "Passed all tests!";

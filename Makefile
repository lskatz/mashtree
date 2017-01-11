# Author: Lee Katz <lkatz@cdc.gov>
# Mashtree
 
PROFILE := $(HOME)/.bashrc
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

###################################

.DEFAULT: install

.PHONY: core extras install install-perlModules install-config core extras

.DELETE_ON_ERROR:

install: perl-modules

perl-modules:
	@perl -Mthreads -e 1 || (echo "ERROR: this version of perl does not have threads" && exit 1)
	cpanm -v -L . Bio::Tree::DistanceFactory
	cpanm -v -L . Bio::Matrix::IO
	cpanm -v -L . Bio::Tree::Statistics
	cpanm -v -L . DBD::SQLite

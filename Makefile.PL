use strict;
use warnings;
 
use 5.008;
 
use ExtUtils::MakeMaker;
WriteMakefile
(
  NAME         => 'Mashtree',
  VERSION_FROM => 'lib/Mashtree.pm',
  EXE_FILES    => [
    "bin/mashtree.pl",
    "bin/mashtree_wrapper.pl",
    "bin/min_abundance_finder.pl",
  ],
  PREREQ_PM    => {
    # Core modules
    'File::Basename' => 0,
    'Data::Dumper'   => 0,
    'List::Util'     => 0,
    'Exporter'       => 0,
    # Threads modules
    'threads'        => 0,
    'threads::shared'=> 0,
    'Thread::Queue'  => 0,
    # Not core (I think?)
    'DBI'            => 0,
    'DBD::SQLite'    => 0,
    'Bio::Matrix::IO'=> 0,
    'Bio::Tree::DistanceFactory'=> 0,
    'Bio::Tree::Statistics'=> 0,
    'Bio::Kmer'      => 0.19,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
);


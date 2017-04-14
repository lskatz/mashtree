use strict;
use warnings;
 
use 5.008;
 
use ExtUtils::MakeMaker;
WriteMakefile
(
  NAME         => 'Mashtree',
  VERSION_FROM => 'lib/perl5/Mashtree.pm',
  PREREQ_PM    => {
    # Core modules
    'File::Basename' => 1,
    'Data::Dumper'   => 1,
    'List::Util'     => 1,
    'Exporter'       => 1,
    # Threads modules
    'threads'        => 1,
    'threads::shared'=> 1,
    'Thread::Queue'  => 1,
    # Not core (I think?)
    'DBI'            => 1,
    'DBD::SQLite'    => 0,
    'Bio::Matrix::IO'=> 0,
    'Bio::Tree::DistanceFactory'=> 0,
    'Bio::Tree::Statistics'=> 0,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
);


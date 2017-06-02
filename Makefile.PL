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
    'Bio::Matrix::IO'=> 0,
    'Bio::Tree::DistanceFactory'=> 0,
    'Bio::Tree::Statistics'=> 0,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
);


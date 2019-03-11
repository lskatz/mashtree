use strict;
use warnings;
 
use 5.010;
 
use ExtUtils::MakeMaker;

# Some help from
# http://blogs.perl.org/users/michal_wojciechowski/2011/11/github-friendly-readme-files-with-extutils-makemaker-and-module-build.html
my $preop = 'true';
  #'sed -i.bak "s/.5f/.10f/g" bin/quicktree-2.3/src/tree.c bin/quicktree-2.3/src/distancemat.c;' .
  #'make -C bin/quicktree-2.3;' .
  #'perldoc -uT $(VERSION_FROM) | tee $(DISTVNAME)/README.pod > README.pod;' .
  #'pod2text README.pod | tee $(DISTVNAME)/README > README';

WriteMakefile1(
  NAME         => 'Mashtree',
  VERSION_FROM => 'lib/Mashtree.pm',
  ABSTRACT_FROM => 'lib/Mashtree.pm',
  AUTHOR => q{Lee S. Katz <gzu2@cdc.gov>},
  EXE_FILES    => [
    "bin/mashtree",
    "bin/mashtree_wrapper_deprecated.pl",
    "bin/mashtree_jackknife.pl",
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
    #'Bio::Kmer'      => 0.19,
    'Graph::Dijkstra'=> 0,
    'JSON'           => '2.90',
    'JSON::XS'       => '4.0',
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
  MIN_PERL_VERSION => '5.10.0',
  LICENSE          => "gpl_3",
  dist=>{
    COMPRESS => "gzip -9f",
    SUFFIX   => "gz",
    PREOP    => $preop,
  },
  META_MERGE => {
    'meta-spec' => { version => 2 },
    resources => {
      repository => {
        type => 'git',
        url => 'git://github.com/lskatz/mashtree.git',
        web => 'https://github.com/lskatz/mashtree',
      },
    },
  },
);

sub WriteMakefile1 {  #Compatibility code for old versions of EU::MM. Written by Alexandr Ciornii, version 2. Added by eumm-upgrade.
  my %params=@_;
  my $eumm_version=$ExtUtils::MakeMaker::VERSION;
  $eumm_version=eval $eumm_version;
  die "EXTRA_META is deprecated" if exists $params{EXTRA_META};
  die "License not specified" if not exists $params{LICENSE};
  if ($params{AUTHOR} and ref($params{AUTHOR}) eq 'ARRAY' and $eumm_version < 6.5705) {
    $params{META_ADD}->{author}=$params{AUTHOR};
    $params{AUTHOR}=join(', ',@{$params{AUTHOR}});
  }
  if ($params{TEST_REQUIRES} and $eumm_version < 6.64) {
    $params{BUILD_REQUIRES}={ %{$params{BUILD_REQUIRES} || {}} , %{$params{TEST_REQUIRES}} };
    delete $params{TEST_REQUIRES};
  }
  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
    #EUMM 6.5502 has problems with BUILD_REQUIRES
    $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
    delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;

  WriteMakefile(%params);
}

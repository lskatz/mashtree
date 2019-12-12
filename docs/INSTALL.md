# Installation

## Prerequisite software

* Mash >= v2.0
* SQLite3
* Perl
  * multithreading
  * BioPerl library
  * `DBD::SQLite`
* Quicktree

### Environment

#### Environment variables

For most Linux OSs, you will need to set up your environment like this:

    export PATH=$HOME/bin:$PATH
    export PERL5LIB=$PERL5LIB:$HOME/lib/perl5

#### System packages

Some system packages are needed.  On Ubuntu, this is how you might install these packages.

    apt update
    apt install -y build-essential cpanminus libexpat1-dev wget sqlite3 libsqlite3-dev

#### Perl packages

There are some perl packages too and so this is how you would install those on most Linux OSs:

    cpanm -l ~ --notest BioPerl Bio::Sketch::Mash DBD::SQLite DBI 

#### Quicktree

    mkdir -pv $HOME/bin/build
    cd $HOME/bin/build
    wget https://github.com/khowe/quicktree/archive/v2.5.tar.gz
    tar xvf v2.5.tar.gz 
    cd quicktree-2.5
    make
    mv quicktree $HOME/bin/

#### Mash

    mkdir -pv $HOME/bin/build
    cd $HOME/bin/build
    wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar
    tar xvf mash-Linux64-v2.2.tar
    mv -v mash-Linux64-v2.2/mash $HOME/bin/

## Mashtree Installation from CPAN

Installing from CPAN installs the latest stable version of Mashtree.  This method adds the Mashtree perl modules to the correct place in your home directory and adds the executables to your home bin directory.

    export PERL5LIB=$PERL5LIB:$HOME/lib/perl5
    cpanm -l ~ Mashtree
    mashtree --help # verify it shows usage and not an error

## Uninstallation from CPAN

    cpanm --uninstall Mashtree --local-lib=$HOME

## Other sources for Mashtree

Other instances of Mashtree can be found in the wild. Although I did not create these,
others have found them useful. I cannot provide support for these outside instances.

### Docker

* https://hub.docker.com/r/bioedge/mashtree
* https://hub.docker.com/r/staphb/mashtree

### Conda

* https://anaconda.org/bioconda/mashtree



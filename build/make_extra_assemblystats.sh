#!/bin/bash

# This script builds a portable assembly-stats binary,
# using the Holy build box environment to create the portable binaries.
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)

# The Holy build box environment requires docker, make sure it is installed
# e.g. on ubuntu and derivatives
# sudo apt install docker.io
# Also make sure you have permission to use docker
# sudo usermod -a -G docker $USER

# stop on error
set -e

# Prepare and start docker with Holy Build box
# ============================================

script="$(readlink -f "$0")"
dir="$(dirname "$script")"
source "${dir}/start_hbb.sh"

# settings
# ========

# release version does not work with musl, used master version (at given patch) instead (next)
assembly_statsversion=1.0.1

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# does not compile with these settings (X)
# only use HBB for glibc compat, not static libs
# source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall cmake

function download {
    cd /build
    url=$1
    if [ "$2" = "" ] ; then
        filename=$(basename $url)
    else
        filename="$2"
    fi
    if [ ! -f $filename ] ; then
        wget -c -O $filename $url
    fi
    ext=${filename##*.}
    if [ "$ext" = "bz2" ] ; then
        tar xvjf $filename
    elif [ "$ext" = "xz" ] ; then
        tar xvJf $filename
    else
        tar xvzf $filename
    fi
}


# assembly-stats
# ---
rm -rf /build/assembly_stats-$assembly_statsversion || true
mkdir /build/assembly_stats-$assembly_statsversion
cd /build/assembly_stats-$assembly_statsversion
wget https://github.com/sanger-pathogens/assembly-stats/archive/refs/tags/v1.0.1-docker1.tar.gz
tar xvzf v1.0.1-docker1.tar.gz
cd assembly-stats-$assembly_statsversion*
mkdir build
cd build
cmake ..
make
make test

cp -af ./assembly-stats /build/assembly-stats-$assembly_statsversion-$arch
cd /build
ln -s assembly-stats-$assembly_statsversion-$arch assembly-stats

echo "Finished building $builddir/assembly-stats-$assembly_statsversion-$arch"


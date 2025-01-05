#!/bin/bash

# This script downloads a build of pandoc, and puts it in a portable directory
# 32 bit option is ignored
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
pandocversion=3.1.6

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
yuminstall wget

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

# pandoc
# ------
if [ ! "$arch" = "linux-x86_64"] ; then
	echo "Only linux-x86_64 supported currently"
	exit 1
fi

cd /build

wget https://github.com/jgm/pandoc/releases/download/$pandocversion/pandoc-$pandocversion-linux-amd64.tar.gz
tar xvzf pandoc-$pandocversion-linux-amd64.tar.gz
mv pandoc-3.1.6 pandoc-3.1.6-$arch
ln -sf pandoc-3.1.6-$arch/bin/pandoc .
rm pandoc-3.1.6-linux-amd64.tar.gz

echo "Finished building pandoc-$pandocversion-$arch"


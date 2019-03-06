#!/bin/bash

# This script builds the genomecomb extension using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build in (default ~/build/tcl$arch)
# -v|-tclversion|--tclversion: tcl version (default 8.5.19)

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

# Parse arguments
# ===============

tclversion=8.5.19
clean=1
while [[ "$#" -gt 0 ]]; do case $1 in
	-v|-tclversion|--tclversion) tclversion="$2"; shift;;
	-c|-clean|--clean) clean="$2"; shift;;
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# Tk does not compile with these settings (X)
# only use HBB for glibc compat, not static libs
# source /hbb_shlib/activate

# print all executed commands to the terminal
set -x

# set up environment
# ------------------

# locations
tcldir=/build/tcl$tclversion
tkdir=/build/tk$tclversion
dirtcldir=/build/dirtcl$tclversion-$arch
destdir=$dirtcldir/exts

# put dirtcl tclsh in PATH
mkdir /build/bin || true
cd /build/bin
ln -sf $dirtcldir/tclsh8.5 .
PATH=/build/bin:$PATH

# Build
# -----

# Compile
cd /io
mkdir -p linux-x86_64
cd linux-x86_64
../build/version.tcl
if [ "$clean" = 1 ] ; then
	make distclean || true
fi
../configure --disable-threads --prefix="$dirtcldir"
make

echo "Finished building genomecomb extension"

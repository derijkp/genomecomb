#!/bin/bash

# This script builds the genomecomb extension using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build in (default ~/build/tcl$arch)
# -v|-tclversion|--tclversion: tcl version (default 8.5.19)
# -d|-debug|--debug: 1 to make binaries with debug info (-g) 0

# stop on error
set -e

script="$(readlink -f "$0")"
dir="$(dirname "$script")"

# Parse arguments
# ===============

tclversion=8.5.19
clean=1
debug=0
while [[ "$#" -gt 0 ]]; do case $1 in
	-v|-tclversion|--tclversion) tclversion="$2"; shift;;
	-c|-clean|--clean) clean="$2"; shift;;
	-d|-debug|--debug) debug="$2"; shift;;
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

# set up environment
# ------------------

# locations
tcldir=/build/tcl$tclversion
tkdir=/build/tk$tclversion
dirtcldir=/build/dirtcl$tclversion-$arch
destdir=$dirtcldir/exts

# Build
# -----

# Compile
cd $dir
mkdir -p linux-$arch
cd linux-$arch
../build/version.tcl
if [ "$clean" = 1 ] ; then
	make distclean || true
	if [ "$debug" = 1 ] ; then
		../configure --enable-symbols --disable-threads --prefix="$dirtcldir"
	else
		../configure --disable-threads --prefix="$dirtcldir"
	fi
fi

# make
make

echo "Finished building genomecomb extension"

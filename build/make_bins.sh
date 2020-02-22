#!/bin/bash

# This script builds the genomecomb binaries using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)
# -c|-clean|--clean: 1 to do "make clean" before building (default), use 0 to not run make clean
# -d|-debug|--debug: 1 to make binaries with debug info (-g) 0

# The Holy build box environment requires docker, make sure it is installed
# e.g. on ubuntu and derivatives
# sudo apt install docker.io
# Also make sure you have permission to use docker
# sudo usermod -a -G docker $USER

## If it does not find the required zstd libraries in builddir, it will make it

# stop on error
set -e

# Prepare and start docker with Holy Build box
# ============================================

script="$(readlink -f "$0")"
dir="$(dirname "$script")"
source "${dir}/start_hbb.sh"

# Parse arguments
# ===============

clean=1
strip=1
debug=0
while [[ "$#" -gt 0 ]]; do case $1 in
	-c|-clean|--clean) clean="$2"; shift;;
	-s|-strip|--strip) strip="$2"; shift;;
	-d|-debug|--debug) debug="$2"; shift;;
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# Tk does not compile with these settings (X)
# only use HBB for glibc compat, not static libs
source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# set up environment
# ------------------

# Deps
# ----
cd /build
if [ ! -f /build/zstd-1.4.4/lib/libzstd.a ] ; then
	source /hbb_shlib/activate
	CFLAGS="-g -O2 -fPIC -fvisibility=hidden -I/hbb_shlib/include"
	CXXFLAGS="-g -O2 -fPIC -fvisibility=hidden -I/hbb_shlib/include"
	SHLIB_CFLAGS="-g -O2 -fPIC -fvisibility=hidden -I/hbb_shlib/include"
	curl -O -L https://github.com/facebook/zstd/releases/download/v1.4.4/zstd-1.4.4.tar.gz
	tar xvzf zstd-1.4.4.tar.gz
	cd zstd-1.4.4
	make
	source /hbb_exe/activate
fi

# Build
# -----

# Compile
cd /io/src
if [ "$clean" = 1 ] ; then
	make clean
fi

if [ "$debug" = 1 ] ; then
	COPT="-g" CPATH="/build/zstd-1.4.4/lib:$CPATH" LIBRARY_PATH="/build/zstd-1.4.4/lib:$LIBRARY_PATH" make
else
	CPATH="/build/zstd-1.4.4/lib:$CPATH" LIBRARY_PATH="/build/zstd-1.4.4/lib:$LIBRARY_PATH" make
fi

if [ "$strip" = 1 ] && [ "$debug" != 1 ] ; then
	strip /io/bin$ARCH/*
fi

echo "Finished building genomecomb binaries"

#!/bin/bash

# This script builds the genomecomb extension using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build in (default ~/build/tcl$arch)
# -v|-tclversion|--tclversion: tcl version (default 8.5.19)
# -d|-debug|--debug: 1 to make binaries with debug info (-g) 0

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
debug=0
while [[ "$#" -gt 0 ]]; do case $1 in
	-v|-tclversion|--tclversion) tclversion="$2"; shift;;
	-c|-clean|--clean) clean="$2"; shift;;
	-d|-debug|--debug) debug="$2"; shift;;
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

yuminstall centos-release-scl
sudo yum upgrade -y
yuminstall devtoolset-9
# use source instead of scl enable so it can run in a script
source /opt/rh/devtoolset-9/enable

# locations
base=/build
tcldir=$base/tcl$tclversion
tkdir=$base/tk$tclversion
dirtcldir=$base/dirtcl$tclversion-$arch
destdir=$dirtcldir/exts

# put dirtcl tclsh in PATH
mkdir $base/bin || true
cd $base/bin
ln -sf $dirtcldir/tclsh8.5 .
ln -sf $dirtcldir/tclsh8.5 tclsh
PATH=$base/bin:$PATH

# Build
# -----

# Compile
mkdir -p /io/$arch
cd /io/$arch
../build/version.tcl
if [ "$clean" = 1 ] ; then
	make distclean || true
	if [ "$debug" = 1 ] ; then
		../configure --enable-symbols --disable-threads --prefix="$dirtcldir" --exec-prefix="$dirtcldir"
	else
		../configure --disable-threads --prefix="$dirtcldir"
	fi
fi

# make
make

echo "Finished building genomecomb extension"

#!/bin/bash

# This script builds portable pigz binaries using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)
# -a|-all|--all: if 1 (default) all binaries are (re)build, if 0, only the ones missing are build

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
source "${dir}/start_hbb3.sh"

# Parse arguments
# ===============

all=1
while [[ "$#" -gt 0 ]]; do case $1 in
	-a|-all|--all) all="$2"; shift;;
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

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall zlib
yuminstall zlib-devel
yuminstall gcc-c++
yuminstall centos-release-scl
yuminstall cmake3
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-11
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-11 bash
source /opt/rh/devtoolset-11/enable

for dir in lib include bin share ; do
	echo $dir
	mkdir /build/$dir || true
	sudo rmdir /usr/local/$dir || true
	sudo rm /usr/local/$dir || true
	sudo ln -s /build/$dir /usr/local/$dir
done

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

# zlib
# ---
if [ $all = 1 ] || [ ! -f /build/lib/libz.so.1.2.11 ] ; then
	zlibversion=1.3.1
	download https://zlib.net/zlib-$zlibversion.tar.gz
	cd /build/zlib-$zlibversion
	make distclean
	env CFLAGS="-O3 -fPIC" \
            ./configure --prefix=/build
	make
	make install
fi

# ultra
# -----
ultraversion=1.2.1
cd /build
wget https://github.com/TravisWheelerLab/ULTRA/archive/refs/tags/$ultraversion.tar.gz
tar xvzf $ultraversion.tar.gz
rm $ultraversion.tar.gz
cd /build/ULTRA-$ultraversion
cmake3 .
make

cp ultra /io/extern$ARCH/ultra-$ultraversion-linux-x86_64
strip /io/extern$ARCH/ultra-$ultraversion-linux-x86_64
cd /io/extern$ARCH
ln -sf ultra-$ultraversion-linux-x86_64 ultra-$ultraversion
ln -sf ultra-$ultraversion-linux-x86_64 ultra

echo "Finished installing extra external binary extern$ARCH/ultra"

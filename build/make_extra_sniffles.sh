#!/bin/bash

# This script builds portable sniffles binaries using the Holy build box environment
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
source "${dir}/start_hbb.sh"

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
yuminstall xz

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

# sniffles
# --------
if [ $all = 1 ] || [ ! -f /io/extra$ARCH/snifles ] ; then
	yuminstall cmake
	# yuminstall openmpi
	yuminstall gcc-c++
	# sniffles
	# --------
	snifflesversion=1.0.11
	download https://github.com/fritzsedlazeck/Sniffles/archive/$snifflesversion.tar.gz
	mv $snifflesversion.tar.gz sniffles-$snifflesversion.tar.gz
	cd /build/Sniffles-$snifflesversion
	rm -rf /build/Sniffles-$snifflesversion/build
	mkdir /build/Sniffles-$snifflesversion/build
	cd /build/Sniffles-$snifflesversion/build
	cmake ..
	make
	cp ../bin/sniffles-core-$snifflesversion/sniffles /io/extra$ARCH/sniffles-$snifflesversion
	cd /io/extra$ARCH
	ln -sf sniffles-$snifflesversion sniffles
fi

echo "Finished building extra external binary extra$ARCH/sniffles-$snifflesversion"

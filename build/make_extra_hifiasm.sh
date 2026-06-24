#!/bin/bash

# This script builds portable hifiasm binaries using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)
# -a|-all|--all: if 1 (default) all binaries are (re)build, if 0, only the ones missing in the extern dir are build

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

hifiasmversion=0.25.0

all=1
extra=1
while [[ "$#" -gt 0 ]]; do case $1 in
	-a|-all|--all) all="$2"; shift;;
	-e|-extra|--extra) extra="$2"; shift;;
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall gcc-c++
yuminstall centos-release-scl
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

cd /build

git clone https://github.com/chhylp123/hifiasm
mv hifiasm hifiasm-build-$hifiasmversion
cd hifiasm-build-$hifiasmversion
git checkout $hifiasmversion
make

cp -a hifiasm /build/hifiasm-$hifiasmversion-$arch
cd /build
ln -s hifiasm-$hifiasmversion-$arch hifiasm

rm hifiasm-$hifiasmversion-$arch.tar.gz || true
tar cvzf hifiasm-$hifiasmversion-$arch.tar.gz hifiasm-$hifiasmversion-$arch hifiasm
rm -rf /io/extra$ARCH/hifiasm-$hifiasmversion-$arch
cp -ra hifiasm-$hifiasmversion-$arch hifiasm hifiasm.py hifiasm3 gtf2db hifiasm_gtf2db hifiasm3_gtf2db /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building hifiasm-$hifiasmversion-$arch"

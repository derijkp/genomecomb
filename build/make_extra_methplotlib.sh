#!/bin/bash

# This script builds portable methplotlib binaries using the Holy build box environment
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
source "${dir}/start_hbb.sh"

# Parse arguments
# ===============

methplotlibversion=0.20.1

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
yuminstall devtoolset-9
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-9 bash
source /opt/rh/devtoolset-9/enable


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

# mamba
# -----
cd /build
export mambaversion=22.11.1-4
curl -L -O "https://github.com/conda-forge/miniforge/releases/download/$mambaversion/Mambaforge-$mambaversion-Linux-x86_64.sh"
unset PYTHONPATH
rm -rf /home/build/mambaforge
bash Mambaforge-$mambaversion-Linux-x86_64.sh -b

# bioconda
# --------

PATH=/home/build/mambaforge/bin:$PATH

mamba init bash
. ~/.bash_profile

# methplotlib
# --------
cd /build

mamba create -y -n methplotlib
mamba activate methplotlib
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y methplotlib=$methplotlibversion

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm methplotlib.tar.gz || true
conda pack -n methplotlib -o methplotlib.tar.gz
rm -rf methplotlib-$methplotlibversion-$arch.old || true
mv methplotlib-$methplotlibversion-$arch methplotlib-$methplotlibversion-$arch.old || true

mkdir /build/methplotlib-$methplotlibversion-$arch
cd /build/methplotlib-$methplotlibversion-$arch
tar xvzf ../methplotlib.tar.gz

cd /build/methplotlib-$methplotlibversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/methplotlib ${1+"$@"}
' > methplotlib
chmod ugo+x methplotlib

cd /build
ln -sf methplotlib-$methplotlibversion-$arch/methplotlib .
ln -sf methplotlib-$methplotlibversion-$arch/methplotlib methplotlib-$methplotlibversion
rm methplotlib-$methplotlibversion-$arch.tar.gz || true
tar cvzf methplotlib-$methplotlibversion-$arch.tar.gz methplotlib-$methplotlibversion-$arch methplotlib methplotlib-$methplotlibversion
rm -rf /io/extra$ARCH/methplotlib-$methplotlibversion-$arch
cp -ra methplotlib-$methplotlibversion-$arch methplotlib methplotlib-$methplotlibversion /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building methplotlib-$methplotlibversion-$arch"

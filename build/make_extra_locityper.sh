#!/bin/bash

# This script builds portable locityper binaries using the Holy build box environment
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

locityperversion=1.3.4
conda_locityperversion=1.3.4

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

# locityper
# --------
cd /build

mamba create -y -n locityper
mamba activate locityper
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y locityper=$conda_locityperversion python=3.12

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm locityper.tar.gz || true
conda pack -n locityper -o locityper.tar.gz
rm -rf locityper-$locityperversion-$arch.old || true
mv locityper-$locityperversion-$arch locityper-$locityperversion-$arch.old || true

mkdir /build/locityper-$locityperversion-$arch
cd /build/locityper-$locityperversion-$arch
tar xvzf ../locityper.tar.gz

cd /build/locityper-$locityperversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/locityper ${1+"$@"}
' > locityper
chmod ugo+x locityper

cd /build
ln -sf locityper-$locityperversion-$arch/locityper .
rm locityper-$locityperversion-$arch.tar.gz || true
tar cvzf locityper-$locityperversion-$arch.tar.gz locityper locityper.py locityper3 locityper-$locityperversion-$arch gtf2db locityper_gtf2db locityper3_gtf2db
rm -rf /io/extra$ARCH/locityper-$locityperversion-$arch
cp -ra locityper-$locityperversion-$arch locityper /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building locityper-$locityperversion-$arch"

# Locityper enables targeted genotyping of complex polymorphic genes
# https://pmc.ncbi.nlm.nih.gov/articles/PMC12597825/
# https://github.com/tprodanov/locityper

# Long-read reconstruction of many diverse haplotypes with devider
# https://pmc.ncbi.nlm.nih.gov/articles/PMC12642997/
# https://github.com/bluenote-1577/devider

https://github.com/bluenote-1577/floria

#!/bin/bash

# This script builds portable flye binaries using the Holy build box environment
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

flyeversion=2.9.2

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

# miniconda
# ---------
export minicondaversion=py38_4.9.2
cd /build
wget -c https://repo.anaconda.com/miniconda/Miniconda3-$minicondaversion-Linux-x86_64.sh
unset PYTHONPATH
rm -rf /build/miniconda-$minicondaversion
bash Miniconda3-$minicondaversion-Linux-x86_64.sh -b -p /build/miniconda-$minicondaversion

# bioconda
# --------

PATH=/build/miniconda-$minicondaversion/bin:$PATH

conda init bash
. ~/.bash_profile

# flye
# -----
cd /build

conda create -y -n flye
conda activate flye
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y flye=$flyeversion python=3.9

conda deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm flye.tar.gz || true
conda pack -n flye -o flye.tar.gz
rm -rf flye-$flyeversion-$arch.old || true
mv flye-$flyeversion-$arch flye-$flyeversion-$arch.old || true

mkdir /build/flye-$flyeversion-$arch
cd /build/flye-$flyeversion-$arch
tar xvzf ../flye.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/flye ${1+"$@"}
' > flye
chmod ugo+x flye

cd /build
ln -sf flye-$flyeversion-$arch/flye .
rm flye-$flyeversion-$arch.tar.gz || true
tar cvzf flye-$flyeversion-$arch.tar.gz flye-$flyeversion-$arch flye
rm -rf /io/extra$ARCH/flye-$flyeversion-$arch
cp -ra flye-$flyeversion-$arch flye /io/extra$ARCH

echo "Finished building flye-$flyeversion-$arch"

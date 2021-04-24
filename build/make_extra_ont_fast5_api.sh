#!/bin/bash

# This script builds portable ont_fast5 binaries using the Holy build box environment
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
yuminstall centos-release-scl
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-8
yuminstall rh-python36
# scl enable devtoolset-8 bash
scl enable devtoolset-8 rh-python36 bash


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

# miniconda
# ---------
cd /build
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
unset PYTHONPATH
bash Miniconda3-latest-Linux-x86_64.sh -b -p /build/miniconda

# bioconda
# --------

PATH=/build/miniconda/bin:$PATH

conda install -y -c conda-forge conda-pack

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda init bash
. ~/.bash_profile

conda create -y -n ont-fast5-api
conda activate ont-fast5-api
conda install -y ont-fast5-api

conda pack -n ont-fast5-api -o ont-fast5-api.tar.gz
mkdir ont-fast5-api
cd ont-fast5-api
tar xvzf ../ont-fast5-api.tar.gz
cd ..
rm ont-fast5-api.tar.gz
tar cvzf ont-fast5-api.tar.gz ont-fast5-api

conda deactivate

echo "Finished building ont-fast5-api"

fi # end of extra
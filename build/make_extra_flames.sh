#!/bin/bash

# This script builds portable flames binaries using the Holy build box environment
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

# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02525-6
# based on https://github.com/LuyiTian/FLAMES, but
# https://github.com/LuyiTian/FLAMES-1 forked from https://github.com/OliverVoogd/FLAMES

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall gcc-c++
yuminstall centos-release-scl
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-8
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-8 rh-python36 bash
source /opt/rh/devtoolset-8/enable


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
wget https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh
unset PYTHONPATH
rm -rf /home/build/miniconda || true
bash Miniconda2-py27_4.8.3-Linux-x86_64.sh -b -p /home/build/miniconda

# bioconda
# --------

PATH=/home/build/miniconda/bin:$PATH

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda init bash
. ~/.bash_profile

# flames
# -----
cd /build

flamesversion=0.1

conda create -y -n flames
conda activate flames

# install dependencies
conda install -y python=2.7
conda install -y samtools pysam minimap2 numpy editdistance

git clone https://github.com/LuyiTian/FLAMES.git
cd /build/FLAMES/src
# g++ -std=c++11 $STATICLIB_CFLAGS $CFLAGS $LDFLAGS -lz -O2 -o match_cell_barcode ssw/ssw_cpp.cpp ssw/ssw.c match_cell_barcode.cpp kseq.h edit_dist.cpp
g++ -std=c++11 -O2 -fvisibility=hidden -I/hbb_exe/include -L/hbb_exe/lib -static-libstdc++ -O2 -o match_cell_barcode ssw/ssw_cpp.cpp ssw/ssw.c match_cell_barcode.cpp kseq.h edit_dist.cpp -lz -lm -lc
mv match_cell_barcode /build/FLAMES

# make portable appdir
conda install -y conda-pack
cd /build
rm flames.tar.gz || true
conda pack -n flames -o flames.tar.gz
rm -rf flames-$flamesversion.old
mv flames-$flamesversion flames-$flamesversion.old || true
mkdir flames-$flamesversion
cd flames-$flamesversion
tar xvzf ../flames.tar.gz
# add FLAMES
cp -ra /build/FLAMES .

# make excutables in appdir root that will take use the appdir env
cd /build/flames-$flamesversion

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$dir/FLAMES/python:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/FLAMES/python/bulk_long_pipeline.py ${1+"$@"}
' > bulk_long_pipeline.py
chmod ugo+x bulk_long_pipeline.py
rm flames_bulk_long_pipeline.py
ln -s bulk_long_pipeline.py flames_bulk_long_pipeline.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$dir/FLAMES/python:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/FLAMES/python/bulk_long_pipeline.py ${1+"$@"}
' > sc_long_pipeline.py
chmod ugo+x sc_long_pipeline.py
rm flames_sc_long_pipeline.py
ln -s sc_long_pipeline.py flames_sc_long_pipeline.py

rm ../flames.tar.gz
cd /build
tar cvzf flames-$flamesversion.tar.gz flames-$flamesversion
cp -ra flames-$flamesversion /io/extra$ARCH
cd /io/extra$ARCH/

conda deactivate

echo "Finished building flames-$flamesversion"

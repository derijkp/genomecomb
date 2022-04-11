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

# print all executed commands to the terminal
set -x

# Activate Holy Build Box environment.
source /hbb_exe/activate

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall centos-release-scl
yuminstall hdf5
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-9
yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-9 rh-python36 bash
source /opt/rh/devtoolset-9/enable
source /opt/rh/rh-python36/enable


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

# conda info -s

# ont-fast5-api
# -------------
apiversion=4.0.2
conda create -y -n ont-fast5-api
conda activate ont-fast5-api
conda config --add channels bioconda
conda install -y ont-fast5-api=$apiversion python=3.8 -c conda-forge

# /build/miniconda-py38_4.9.2/envs/ont-fast5-api/bin/python3.9 pip install

conda deactivate

# create package
# --------------
cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm ont-fast5-api.tar.gz || true
conda pack -n ont-fast5-api -o ont-fast5-api.tar.gz

rm -rf ont-fast5-api || true
rm -rf ont-fast5-api-$apiversion || true
mkdir ont-fast5-api-$apiversion
cd ont-fast5-api-$apiversion
tar xvzf ../ont-fast5-api.tar.gz
cp -a /usr/bin/h5* bin
cp -a /usr/lib64/libhdf5*.so.* lib

cd /build/ont-fast5-api-$apiversion/bin

# fix demux
mv demux_fast5 demux_fast5.ori || true
echo '#!/usr/bin/env python3.8' > demux_fast5
tail -7 demux_fast5.ori >> demux_fast5

for binary in *fast5* h5ls h5dump h5stat h5diff h5format_convert ; do

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
' > ../$binary
echo "\$dir/bin/$binary \${1+\"\$@\"}" >> ../$binary
chmod ugo+x ../$binary

done

cd /build
rm ont-fast5-api-$apiversion.tar.gz || true
tar cvzf ont-fast5-api-$apiversion.tar.gz ont-fast5-api-$apiversion

echo "Finished building ont-fast5-api"

# end of extra
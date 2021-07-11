#!/bin/bash

# This script builds portable medaka binaries using the Holy build box environment
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
yuminstall gcc-c++
yuminstall centos-release-scl
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-8
yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-8 rh-python36 bash
source /opt/rh/devtoolset-8/enable
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

cd /build

# miniconda
# ---------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
unset PYTHONPATH
rm -rf /build/miniconda || true
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

# medaka
# -----
cd /build

medakaversion=1.4.3

conda create -y -n medaka
conda activate medaka
conda install -y medaka=$medakaversion

rm medaka.tar.gz || true
conda pack -n medaka -o medaka.tar.gz
rm -rf medaka-$medakaversion.old
mv medaka-$medakaversion medaka-$medakaversion.old || true

mkdir /build/medaka-$medakaversion
cd /build/medaka-$medakaversion
tar xvzf ../medaka.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka ${1+"$@"}
' > medaka
chmod ugo+x medaka

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_consensus ${1+"$@"}
' > medaka_consensus
chmod ugo+x medaka_consensus

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_counts ${1+"$@"}
' > medaka_counts
chmod ugo+x medaka_counts

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_data_path ${1+"$@"}
' > medaka_data_path
chmod ugo+x medaka_data_path

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_haploid_variant ${1+"$@"}
' > medaka_haploid_variant
chmod ugo+x medaka_haploid_variant

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_variant ${1+"$@"}
' > medaka_variant
chmod ugo+x medaka_variant

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
medaka_version_report ${1+"$@"}
' > medaka_version_report
chmod ugo+x medaka_version_report

# some fixes
cp /io/extern-src/medaka_whatshap /build/medaka-$medakaversion/bin/whatshap

rm ../medaka.tar.gz
cd /build
tar cvzf medaka-$medakaversion.tar.gz medaka-$medakaversion
cp -ra medaka-$medakaversion /io/extra$ARCH
cd /io/extra$ARCH/

conda deactivate

echo "Finished building medaka"

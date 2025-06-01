#!/bin/bash

# This script builds portable annotsv binaries using the Holy build box environment
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

annotsvversion=3.4.4

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
# only use HBB for glibc compat, not static libs
# source /hbb_exe/activate

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

yuminstall tcl
yuminstall bzip2-devel
yuminstall libstdc++-devel

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

# clair3
# -----
cd /build

mamba create -y -n annotsv
mamba activate annotsv
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y annotsv=$annotsvversion

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm annotsv.tar.gz || true
conda pack -n annotsv -o annotsv.tar.gz
rm -rf AnnotSV-$annotsvversion-$arch.old || true
mv AnnotSV-$annotsvversion-$arch annotsv-$annotsvversion-$arch.old || true
mkdir /build/AnnotSV-$annotsvversion-$arch
cd /build/AnnotSV-$annotsvversion-$arch
tar xvzf ../annotsv.tar.gz

# install from source because conda install misses some files (application.properties)
wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.4.6.tar.gz
tar xvzf v3.4.6.tar.gz
rm v3.4.6.tar.gz
cd AnnotSV-3.4.6
make PREFIX=/build/AnnotSV-$annotsvversion-$arch install
make PREFIX=/build/AnnotSV-$annotsvversion-$arch install-human-annotation
cd /build
rm -rf /build/AnnotSV-3.4.6

cd /build/AnnotSV-$annotsvversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
export ANNOTSV=$dir
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/AnnotSV ${1+"$@"}
' > AnnotSV
chmod ugo+x AnnotSV

# make distribution
cd /build
ln -sf AnnotSV-$annotsvversion-$arch/AnnotSV .
ln -sf AnnotSV-$annotsvversion-$arch/AnnotSV AnnotSV-$annotsvversion
rm AnnotSV-$annotsvversion-$arch.tar.gz || true
tar cvzf AnnotSV-$annotsvversion-$arch.tar.gz AnnotSV-$annotsvversion-$arch AnnotSV
rm -rf /io/extra$ARCH/AnnotSV-$annotsvversion-$arch
cp -ra AnnotSV-$annotsvversion-$arch AnnotSV-$annotsvversion AnnotSV /io/extra$ARCH

echo "Running test"
cd /tmp
cp /build/AnnotSV-$annotsvversion-linux-x86_64/share/doc/AnnotSV/Example/test.bed .
/build/AnnotSV-$annotsvversion-linux-x86_64/AnnotSV -SVinputFile test.bed -genomeBuild GRCh37 -outputDir /tmp -outputFile test.annotated.tsv -svtBEDcol 4

echo "Finished building AnnotSV-$annotsvversion-$arch"


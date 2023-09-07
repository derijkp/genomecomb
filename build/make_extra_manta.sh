#!/bin/bash

# This script builds portable manta binaries using the Holy build box environment
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

mantaversion=1.6.0

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

# manta
# -----
cd /build

mamba create -y -n manta
mamba activate manta
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y manta=$mantaversion

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm manta.tar.gz || true
conda pack -n manta -o manta.tar.gz
rm -rf manta-$mantaversion-$arch.old || true
mv manta-$mantaversion-$arch manta-$mantaversion-$arch.old || true

mkdir /build/manta-$mantaversion-$arch
cd /build/manta-$mantaversion-$arch
tar xvzf ../manta.tar.gz

cd /build/manta-$mantaversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/configManta.py ${1+"$@"}
' > configManta.py
chmod ugo+x configManta.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/configManta.py ${1+"$@"}
' > manta
chmod ugo+x manta

cd /build
ln -sf manta-$mantaversion-$arch/configManta.py .
ln -sf manta-$mantaversion-$arch/manta .
rm manta-$mantaversion-$arch.tar.gz || true
tar cvzf manta-$mantaversion-$arch.tar.gz manta-$mantaversion-$arch configManta.py manta
rm -rf /io/extra$ARCH/manta-$mantaversion-$arch
cp -ra manta-$mantaversion-$arch configManta.py /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building manta-$mantaversion-$arch"

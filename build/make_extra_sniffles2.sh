#!/bin/bash

# This script builds portable sniffles2 binaries using the Holy build box environment
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

snifflesversion=2.8.0

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

# sniffles
# -----
cd /build

mamba create -y -n sniffles
mamba activate sniffles
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y sniffles=$snifflesversion python=3.10

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm sniffles.tar.gz || true
conda pack -n sniffles -o sniffles.tar.gz
rm -rf sniffles-$snifflesversion-$arch.old || true
mv sniffles-$snifflesversion-$arch sniffles-$snifflesversion-$arch.old || true
mkdir /build/sniffles-$snifflesversion-$arch
cd /build/sniffles-$snifflesversion-$arch
tar xvzf ../sniffles.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/sniffles ${1+"$@"}
' > sniffles
chmod ugo+x sniffles

# wrap up
rm ../sniffles-$snifflesversion-$arch.tar.gz || true
cd /build
tar cvzf sniffles-$snifflesversion-$arch.tar.gz sniffles-$snifflesversion-$arch
cp -ra sniffles-$snifflesversion-$arch /io/extra$ARCH
cd /io/extra$ARCH/
rm sniffles || true
ln -s sniffles-$snifflesversion-$arch/sniffles .
ln -s sniffles-$snifflesversion-$arch/sniffles sniffles-$snifflesversion

echo "build in $builddir/sniffles-$snifflesversion-$arch"
echo "installed in $srcdir/sniffles-$snifflesversion-$arch"
echo "Finished building sniffles2"

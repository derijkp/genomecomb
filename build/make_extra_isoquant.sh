#!/bin/bash

# This script builds portable isoquant binaries using the Holy build box environment
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

isoquantversion=3.3.0
conda_isoquantversion=3.2.0

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

# isoquant
# --------
cd /build

mamba create -y -n isoquant
mamba activate isoquant
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y isoquant=$conda_isoquantversion python=3.9

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm isoquant.tar.gz || true
conda pack -n isoquant -o isoquant.tar.gz
rm -rf isoquant-$isoquantversion-$arch.old || true
mv isoquant-$isoquantversion-$arch isoquant-$isoquantversion-$arch.old || true

mkdir /build/isoquant-$isoquantversion-$arch
cd /build/isoquant-$isoquantversion-$arch
tar xvzf ../isoquant.tar.gz

# to 3.3.0
cd /build/isoquant-$isoquantversion-$arch/share
wget https://github.com/ablab/IsoQuant/releases/download/v$isoquantversion/IsoQuant-$isoquantversion.tar.gz
tar xvzf IsoQuant-$isoquantversion.tar.gz
rm -rf isoquant-$conda_isoquantversion*
rm IsoQuant-$isoquantversion.tar.gz
mv IsoQuant-$isoquantversion isoquant-$isoquantversion
cd isoquant-$isoquantversion
cd ../bin
rm isoquant.py ; ln -s ../share/isoquant-3.3.0/isoquant.py .

## patch
cd /build/isoquant-$isoquantversion-$arch/share/isoquant-$isoquantversion*/src
cp long_read_counter.py long_read_counter.py.ori || true
cp /io/extern-src/isoquant/long_read_counter.py long_read_counter.py

cd /build/isoquant-$isoquantversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/isoquant.py ${1+"$@"}
' > isoquant.py
chmod ugo+x isoquant.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/isoquant.py ${1+"$@"}
' > isoquant
chmod ugo+x isoquant

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/share/isoquant-*/src/gtf2db.py ${1+"$@"}
' > gtf2db
chmod ugo+x gtf2db

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/share/isoquant-*/src/gtf2db.py ${1+"$@"}
' > isoquant_gtf2db
chmod ugo+x isoquant_gtf2db

cd /build
ln -sf isoquant-$isoquantversion-$arch/isoquant .
ln -sf isoquant-$isoquantversion-$arch/isoquant isoquant3
ln -sf isoquant-$isoquantversion-$arch/isoquant isoquant-$isoquantversion
ln -sf isoquant-$isoquantversion-$arch/isoquant.py .
ln -sf isoquant-$isoquantversion-$arch/gtf2db .
ln -sf isoquant-$isoquantversion-$arch/isoquant_gtf2db .
ln -sf isoquant-$isoquantversion-$arch/isoquant_gtf2db isoquant3_gtf2db
rm isoquant-$isoquantversion-$arch.tar.gz || true
tar cvzf isoquant-$isoquantversion-$arch.tar.gz isoquant isoquant.py isoquant3 isoquant-$isoquantversion-$arch gtf2db isoquant_gtf2db isoquant3_gtf2db
rm -rf /io/extra$ARCH/isoquant-$isoquantversion-$arch
cp -ra isoquant-$isoquantversion-$arch isoquant isoquant.py isoquant3 isoquant-$isoquantversion-$arch gtf2db isoquant_gtf2db isoquant3_gtf2db /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building isoquant-$isoquantversion-$arch"

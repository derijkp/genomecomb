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

# isoquant
# --------
cd /build

isoquantversion=3.1.0

conda create -y -n isoquant
conda activate isoquant
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y isoquant=$isoquantversion python=3.9

conda deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm isoquant.tar.gz || true
conda pack -n isoquant -o isoquant.tar.gz
rm -rf isoquant-$isoquantversion-$arch.old || true
mv isoquant-$isoquantversion-$arch isoquant-$isoquantversion-$arch.old || true

mkdir /build/isoquant-$isoquantversion-$arch
cd /build/isoquant-$isoquantversion-$arch
tar xvzf ../isoquant.tar.gz

# patch
cp /io/extern-src/isoquant/long_read_counter.py share/isoquant-$isoquantversion-0/src/long_read_counter.py

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
$dir/share/isoquant-*-0/src/gtf2db.py ${1+"$@"}
' > gtf2db
chmod ugo+x gtf2db

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/share/isoquant-*-0/src/gtf2db.py ${1+"$@"}
' > isoquant_gtf2db
chmod ugo+x isoquant_gtf2db

cd /build
ln -s isoquant-$isoquantversion-$arch/isoquant .
ln -s isoquant-$isoquantversion-$arch/isoquant isoquant3
ln -s isoquant-$isoquantversion-$arch/isoquant isoquant-$isoquantversion
ln -s isoquant-$isoquantversion-$arch/isoquant.py .
ln -s isoquant-$isoquantversion-$arch/gtf2db .
ln -s isoquant-$isoquantversion-$arch/isoquant_gtf2db .
ln -s isoquant-$isoquantversion-$arch/isoquant_gtf2db isoquant3_gtf2db
rm isoquant-$isoquantversion-$arch.tar.gz || true
tar cvzf isoquant-$isoquantversion-$arch.tar.gz isoquant isoquant.py isoquant3 isoquant-$isoquantversion-$arch gtf2db isoquant_gtf2db isoquant3_gtf2db
cp -ra isoquant-$isoquantversion-$arch isoquant isoquant.py isoquant3 isoquant-$isoquantversion-$arch gtf2db isoquant_gtf2db isoquant3_gtf2db /io/extra$ARCH
cd /io/extra$ARCH/
rm isoquant || true

echo "Finished building isoquant"

#!/bin/bash

# This script builds portable lumpy binaries using the Holy build box environment
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

# print all executed commands to the terminal
set -x

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# Tk does not compile with these settings (X)
# only use HBB for glibc compat, not static libs
source /hbb_exe/activate

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
# yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# instead of: scl enable devtoolset-9 bash
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

PATH=/build/miniconda-$minicondaversion/bin:$PATH

conda init bash
. ~/.bash_profile

# lumpy
# -----
cd /build

lumpyversion=0.3.1

conda create -y -n lumpy
conda activate lumpy
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y python=2.7 lumpy-sv=$lumpyversion

conda deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm lumpy.tar.gz || true
conda pack -n lumpy -o lumpy.tar.gz
rm -rf lumpy-$lumpyversion-$arch.old || true
mv lumpy-$lumpyversion-$arch lumpy-$lumpyversion-$arch.old || true

#wget https://github.com/arq5x/lumpy-sv/archive/refs/tags/v$lumpyversion.tar.gz
#mv v$lumpyversion.tar.gz lumpy-v$lumpyversion.tar.gz

mkdir /build/lumpy-$lumpyversion-$arch
cd /build/lumpy-$lumpyversion-$arch
tar xvzf ../lumpy.tar.gz
#tar xvzf ../lumpy-v$lumpyversion.tar.gz
#mkdir scripts/bamkit

# add libraries
# somehow the binaries are linked to libcrypto.so.1.0.0, as we don't have this exact one in the hbb,
# copy the 1.0* version we do have, and create 1.0.0 as a softlink
cp -a -f /usr/lib64/libcrypto.so.1.0* /build/lumpy-$lumpyversion-$arch/lib
cd /build/lumpy-$lumpyversion-$arch/lib
ln -s libcrypto.so.1.0* libcrypto.so.1.0.0

# add executable wrappers
cd /build/lumpy-$lumpyversion-$arch

echo 'script="$(readlink -f "$0")"
dir="$(dirname "$script")"
dir="$(dirname "$dir")"
LUMPY_HOME=$dir/share/lumpy-sv-0.3.1-3
LUMPY=$dir/bin/lumpy
HEXDUMP=hexdump
SAMBLASTER=samblaster
SAMBAMBA=sambamba
SAMTOOLS=samtools
PYTHON=python
PAIREND_DISTRO=$dir/share/lumpy-sv-0.3.1-3/scripts/pairend_distro.py
BAMGROUPREADS=$dir/bin/bamgroupreads.py
BAMFILTERRG=$dir/bin/bamfilterrg.py
BAMLIBS=$dir/bin/bamlibs.py
' > bin/lumpyexpress.config

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/lumpyexpress ${1+"$@"}
' > lumpyexpress
chmod ugo+x lumpyexpress

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/lumpy ${1+"$@"}
' > lumpy
chmod ugo+x lumpy

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/lumpy_filter ${1+"$@"}
' > lumpy_filter
chmod ugo+x lumpy_filter

# wrap up
rm ../lumpy-$lumpyversion-$arch.tar.gz || true
cd /build
tar cvzf lumpy-$lumpyversion-$arch.tar.gz lumpy-$lumpyversion-$arch
cp -ra lumpy-$lumpyversion-$arch /io/extra$ARCH
cd /io/extra$ARCH/
rm lumpy || true
ln -s lumpy-$lumpyversion-$arch/lumpy .
ln -s lumpy-$lumpyversion-$arch/lumpyexpress .
ln -s lumpy-$lumpyversion-$arch/lumpy_filter .
ln -s lumpy-$lumpyversion-$arch/lumpy lumpy-$lumpyversion
ln -s lumpy-$lumpyversion-$arch/lumpyexpress lumpyexpress-$lumpyversion

echo "build in $builddir/lumpy-$lumpyversion-$arch"
echo "installed in $srcdir/lumpy-$lumpyversion-$arch"
echo "Finished building lumpy"

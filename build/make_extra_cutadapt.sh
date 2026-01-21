#!/bin/bash

# This script builds portable cutadapt binaries using the Holy build box environment
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

cutadaptversion=2.8

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
#source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
# yuminstall gcc-c++
yuminstall centos-release-scl
yuminstall zstd
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
export mambaversion=25.3.0-3
curl -L -O "https://github.com/conda-forge/miniforge/releases/download/$mambaversion/Miniforge3-$mambaversion-Linux-x86_64.sh"
unset PYTHONPATH
rm -rf /home/build/miniforge
bash Miniforge3-$mambaversion-Linux-x86_64.sh -b

# bioconda
# --------

PATH=/home/build/miniforge3/bin:$PATH

#mamba init bash
mamba shell init
. ~/.bash_profile


# cutadapt
# -----

cd /build
mamba create -y -n cutadapt -c bioconda -c conda-forge cutadapt=$cutadaptversion

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm cutadapt.tar.gz || true
conda pack -n cutadapt -o cutadapt.tar.gz
rm -rf cutadapt-$cutadaptversion-$arch.old || true
mv cutadapt-$cutadaptversion-$arch cutadapt-$cutadaptversion-$arch.old || true
mkdir /build/cutadapt-$cutadaptversion-$arch
cd /build/cutadapt-$cutadaptversion-$arch
tar xvzf ../cutadapt.tar.gz

cd /build/cutadapt-$cutadaptversion-$arch

mv bin/cutadapt bin/cutadapt.ori
cat << 'EOF' > bin/cutadapt
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys
from cutadapt.__main__ import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
chmod ugo+x bin/cutadapt

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
export PATH=$dir/bin:$PATH
export LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
export LANG=C
export LC_ALL=C
$dir/bin/cutadapt ${1+"$@"}
' > cutadapt
chmod ugo+x cutadapt

cd /build
ln -sf cutadapt-$cutadaptversion-$arch/cutadapt .
ln -sf cutadapt-$cutadaptversion-$arch/cutadapt cutadapt-$cutadaptversion
rm cutadapt-$cutadaptversion-$arch.tar.gz || true
tar cvzf cutadapt-$cutadaptversion-$arch.tar.gz cutadapt-$cutadaptversion-$arch cutadapt
rm -rf /io/extra$ARCH/cutadapt-$cutadaptversion-$arch
cp -ra cutadapt-$cutadaptversion-$arch cutadapt-$cutadaptversion cutadapt /io/extra$ARCH

echo "Finished building cutadapt-$cutadaptversion-$arch"

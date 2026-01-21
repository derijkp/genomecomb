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
source "${dir}/start_hbb3.sh"

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

# pod5
# ----
pod5version=0.3.1
mamba create -y -n pod5
mamba activate pod5
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y python=3.9 -c conda-forge
mamba install -y pip -c conda-forge

~/mambaforge/envs/pod5/bin/pip install pod5

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm pod5.tar.gz || true
conda pack -n pod5 -o pod5.tar.gz
rm -rf pod5-$pod5version-$arch.old || true
mv pod5-$pod5version-$arch pod5-$pod5version-$arch.old || true
mkdir /build/pod5-$pod5version-$arch
cd /build/pod5-$pod5version-$arch
tar xvzf ../pod5.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/pod5 ${1+"$@"}
' > pod5
chmod ugo+x pod5

# package
cd /build
rm pod5.tar.gz
ln -sf pod5-$pod5version-$arch/pod5 pod5
tar cvzf pod5-$pod5version-$arch.tar.gz pod5-$pod5version-$arch pod5
cp -ra pod5-$pod5version-$arch pod5 /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building pod5"

# end of extra
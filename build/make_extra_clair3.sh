#!/bin/bash

# This script builds portable clair3 binaries using the Holy build box environment
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

# clair3
# -----
cd /build

clair3version=0.1.12

conda create -y -n clair3
conda activate clair3
conda install -y clair3=$clair3version

# make package

rm clair3.tar.gz || true
conda pack -n clair3 -o clair3.tar.gz
rm -rf clair3-$clair3version.old
mv clair3-$clair3version clair3-$clair3version.old || true

mkdir /build/clair3-$clair3version
cd /build/clair3-$clair3version
tar xvzf ../clair3.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/run_clair3.sh ${1+"$@"}
' > run_clair3.sh
chmod ugo+x run_clair3.sh

mv bin/whatshap bin/whatshap.ori
cat << 'EOF' > bin/whatshap
#!/bin/sh
'''exec' python "$0" "$@"
' '''
# -*- coding: utf-8 -*-
import re
import sys
from whatshap.__main__ import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
chmod ugo+x bin/whatshap

cd /build/clair3-$clair3version
mkdir models
cd models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar xvzf clair3_models.tar.gz
rm clair3_models.tar.gz

rm ../clair3.tar.gz
cd /build
tar cvzf clair3-$clair3version.tar.gz clair3-$clair3version
cp -ra clair3-$clair3version /io/extra$ARCH
cd /io/extra$ARCH/
ln -s clair3-0.1-r5/run_clair3.sh .

conda deactivate

echo "Finished building clair3"

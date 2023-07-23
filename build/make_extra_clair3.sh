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

# clair3version=0.1.12
clair3version=1.0.3

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
yuminstall devtoolset-9
# yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-8 rh-python36 bash
source /opt/rh/devtoolset-9/enable
#source /opt/rh/rh-python36/enable


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

# clair3
# -----
cd /build

conda create -y -n clair3
conda activate clair3
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y clair3=$clair3version python=3.9

conda deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm clair3.tar.gz || true
conda pack -n clair3 -o clair3.tar.gz
rm -rf clair3-$clair3version-$arch.old || true
mv clair3-$clair3version-$arch clair3-$clair3version-$arch.old || true
mkdir /build/clair3-$clair3version-$arch
cd /build/clair3-$clair3version-$arch
tar xvzf ../clair3.tar.gz

# should be solved in v1.0.0
## patch (for error: `np.int` was a deprecated)
#catch {file copy bin/preprocess/CreateTensorPileupFromCffi.py bin/preprocess/CreateTensorPileupFromCffi.py.ori}
#set c [file_read bin/preprocess/CreateTensorPileupFromCffi.py.ori]
#regsub -all {np.int([^0-9])} $c {int\1} c2
#file_write bin/preprocess/CreateTensorPileupFromCffi.py $c2

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/run_clair3.sh ${1+"$@"}
' > run_clair3.sh
chmod ugo+x run_clair3.sh

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/run_clair3.sh ${1+"$@"}
' > clair3
chmod ugo+x clair3

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

cd /build/clair3-$clair3version-$arch
mkdir models
cd models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz
tar xvzf clair3_models.tar.gz
rm clair3_models.tar.gz

cd /build
ln -sf clair3-$clair3version-$arch/clair3 .
ln -sf clair3-$clair3version-$arch/run_clair3.sh .
rm clair3-$clair3version-$arch.tar.gz || true
tar cvzf clair3-$clair3version-$arch.tar.gz clair3-$clair3version-$arch clair3 run_clair3.sh
rm -rf /io/extra$ARCH/clair3-$clair3version-$arch
cp -ra clair3-$clair3version-$arch clair3-$clair3version-$arch clair3 run_clair3.sh /io/extra$ARCH

echo "Finished building clair3-$clair3version-$arch"

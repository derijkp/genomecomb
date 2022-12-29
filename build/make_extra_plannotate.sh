#!/bin/bash

# This script builds portable pLannotate binaries using the Holy build box environment
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
yuminstall devtoolset-9
yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-8 rh-python36 bash
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

# plannotate
# ----------
cd /build

plannotateversion=1.2.0

conda create -y -n plannotate python=3.9
conda activate plannotate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y plannotate=$plannotateversion
conda install -c conda-forge streamlit
conda update -y streamlit

conda deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack
rm plannotate.tar.gz || true
conda pack -n plannotate -o plannotate.tar.gz
rm -rf plannotate-$plannotateversion.old
mv plannotate-$plannotateversion plannotate-$plannotateversion.old || true

mkdir /build/plannotate-$plannotateversion
cd /build/plannotate-$plannotateversion
tar xvzf ../plannotate.tar.gz

# QnD hack to get base working (not html)
cp /io/extern-src/plannotate/pLannotate.py extra/plannotate-1.2.0/lib/python3.9/site-packages/plannotate/pLannotate.py

mv bin/plannotate bin/plannotate.ori
cat << 'EOF' > bin/plannotate
#!/bin/sh
'''exec' python "$0" "$@"
' '''
# -*- coding: utf-8 -*-
import re
import sys
from plannotate.pLannotate import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
chmod ugo+x bin/plannotate

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/plannotate ${1+"$@"}
' > plannotate
chmod ugo+x plannotate

rm ../plannotate.tar.gz
cd /build
tar cvzf plannotate-$plannotateversion.tar.gz plannotate-$plannotateversion
cp -ra plannotate-$plannotateversion /io/extra$ARCH
cd /io/extra$ARCH/
ln -s plannotate-$plannotateversion/plannotate .

echo "Finished building plannotate"

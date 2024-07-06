#!/bin/bash

# This script builds portable viennarna binaries using the Holy build box environment
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

viennarnaversion=2.6.4-0

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

# viennarna
# --------
cd /build

mamba create -y -n viennarna
mamba activate viennarna
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y viennarna=$viennarnaversion

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm viennarna.tar.gz || true
conda pack -n viennarna -o viennarna.tar.gz
rm -rf viennarna-$viennarnaversion-$arch.old || true
mv viennarna-$viennarnaversion-$arch viennarna-$viennarnaversion-$arch.old || true

mkdir /build/viennarna-$viennarnaversion-$arch
cd /build/viennarna-$viennarnaversion-$arch
tar xvzf ../viennarna.tar.gz

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAfold ${1+"$@"}
' > RNAfold
chmod ugo+x RNAfold

ln -s RNAfold rnafold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNA2Dfold ${1+"$@"}
' > RNA2Dfold
chmod ugo+x RNA2Dfold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAaliduplex ${1+"$@"}
' > RNAaliduplex
chmod ugo+x RNAaliduplex

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAalifold ${1+"$@"}
' > RNAalifold
chmod ugo+x RNAalifold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAcofold ${1+"$@"}
' > RNAcofold
chmod ugo+x RNAcofold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAdistance ${1+"$@"}
' > RNAdistance
chmod ugo+x RNAdistance

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAdistance ${1+"$@"}
' > RNAdistance
chmod ugo+x RNAdistance

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAdos ${1+"$@"}
' > RNAdos
chmod ugo+x RNAdos

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAduplex ${1+"$@"}
' > RNAduplex
chmod ugo+x RNAduplex

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAeval ${1+"$@"}
' > RNAeval
chmod ugo+x RNAeval

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAheat ${1+"$@"}
' > RNAheat
chmod ugo+x RNAheat

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAinverse ${1+"$@"}
' > RNAinverse
chmod ugo+x RNAinverse

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNALalifold ${1+"$@"}
' > RNALalifold
chmod ugo+x RNALalifold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNALfold ${1+"$@"}
' > RNALfold
chmod ugo+x RNALfold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAmultifold ${1+"$@"}
' > RNAmultifold
chmod ugo+x RNAmultifold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNApaln ${1+"$@"}
' > RNApaln
chmod ugo+x RNApaln

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNApdist ${1+"$@"}
' > RNApdist
chmod ugo+x RNApdist

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAparconv ${1+"$@"}
' > RNAparconv
chmod ugo+x RNAparconv

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAPKplex ${1+"$@"}
' > RNAPKplex
chmod ugo+x RNAPKplex

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAplex ${1+"$@"}
' > RNAplex
chmod ugo+x RNAplex

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAplfold ${1+"$@"}
' > RNAplfold
chmod ugo+x RNAplfold

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAplot ${1+"$@"}
' > RNAplot
chmod ugo+x RNAplot

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNApvmin ${1+"$@"}
' > RNApvmin
chmod ugo+x RNApvmin

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAsnoop ${1+"$@"}
' > RNAsnoop
chmod ugo+x RNAsnoop

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAsubopt ${1+"$@"}
' > RNAsubopt
chmod ugo+x RNAsubopt

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/RNAup ${1+"$@"}
' > RNAup
chmod ugo+x RNAup

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/AnalyseSeqs ${1+"$@"}
' > AnalyseSeqs
chmod ugo+x AnalyseSeqs

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/AnalyseDists ${1+"$@"}
' > AnalyseDists
chmod ugo+x AnalyseDists

cd /build
ln -sf viennarna-$viennarnaversion-$arch/RNAfold rnafold
ln -sf viennarna-$viennarnaversion-$arch/RNA* .
ln -sf viennarna-$viennarnaversion-$arch/AnalyseSeqs .
ln -sf viennarna-$viennarnaversion-$arch/AnalyseDists .
rm viennarna-$viennarnaversion-$arch.tar.gz || true
tar cvzf viennarna-$viennarnaversion-$arch.tar.gz rnafold RNA* AnalyseSeqs AnalyseDists
rm -rf /io/extra$ARCH/viennarna-$viennarnaversion-$arch
cp -ra viennarna-$viennarnaversion-$arch rnafold RNA* AnalyseSeqs AnalyseDists /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building viennarna-$viennarnaversion-$arch"


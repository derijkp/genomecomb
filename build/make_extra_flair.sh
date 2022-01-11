#!/bin/bash

# This script builds portable flair binaries using the Holy build box environment
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

# flair
# -----
cd /build

flairversion=1.5

conda create -y -n flair
conda activate flair
conda install -y flair=$flairversion

# make portable appdir
rm flair.tar.gz || true
conda pack -n flair -o flair.tar.gz
rm -rf flair-$flairversion.old
mv flair-$flairversion flair-$flairversion.old || true
mkdir flair-$flairversion
cd flair-$flairversion
tar xvzf ../flair.tar.gz
cp /io/build/flair_files/plot_isoform_usage.patched.py bin/bin

# make excutables in appdir root that will use the appdir env
cd /build/flair-$flairversion

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/flair.py ${1+"$@"}
' > flair.py
chmod ugo+x flair.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/bin/bam2Bed12.py ${1+"$@"}
' > bam2Bed12.py
chmod ugo+x bam2Bed12.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/bin/plot_isoform_usage.py ${1+"$@"}
' > plot_isoform_usage.py
chmod ugo+x plot_isoform_usage.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/bin/predictProductivity.py ${1+"$@"}
' > predictProductivity.py
chmod ugo+x predictProductivity.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/bin/mark_intron_retention.py ${1+"$@"}
' > mark_intron_retention.py
chmod ugo+x mark_intron_retention.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/bin/diff_iso_usage.py ${1+"$@"}
' > diff_iso_usage.py
chmod ugo+x diff_iso_usage.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/bin/diffsplice_fishers_exact.py ${1+"$@"}
' > diffsplice_fishers_exact.py
chmod ugo+x diffsplice_fishers_exact.py

# package
rm ../flair.tar.gz
cd /build
tar cvzf flair-$flairversion.tar.gz flair-$flairversion
cp -ra flair-$flairversion /io/extra$ARCH
cd /io/extra$ARCH/

conda deactivate

echo "Finished building flair"

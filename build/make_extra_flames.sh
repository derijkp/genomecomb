#!/bin/bash

# This script builds portable flames binaries using the Holy build box environment
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
yuminstall devtoolset-8
# yuminstall rh-python36
# use source instead of scl enable so it can run in a script
# scl enable devtoolset-8 rh-python36 bash
source /opt/rh/devtoolset-8/enable
# source /opt/rh/rh-python36/enable


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

# flames
# -----
cd /build

flamescommit=413e09c0e57b554662bfa341d00f191023ac1fa8
flamesversion=c1_413e09c

# conda install

conda create -y -n flames \
    python=2.7 samtools pysam minimap2 numpy editdistance \
    -c bioconda -c conda-forge

# make portable appdir

# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
conda install -y -c conda-forge conda-pack

cd /build
rm flames.tar.gz || true
conda pack -n flames -o flames.tar.gz
rm -rf flames-$flamesversion.old || true
mv flames-$flamesversion flames-$flamesversion.old || true
mkdir /build/flames-$flamesversion
cd /build/flames-$flamesversion
tar xvzf ../flames.tar.gz

# add flames code

wget -c https://github.com/LuyiTian/FLAMES/archive/$flamescommit.zip
unzip $flamescommit.zip
rm $flamescommit.zip
cd /build

# compile single cell code
cd /build/flames-$flamesversion/FLAMES-$flamescommit/src
g++ -std=c++11 -lz -O2 -o ../match_cell_barcode ssw/ssw_cpp.cpp ssw/ssw.c match_cell_barcode.cpp kseq.h edit_dist.cpp
cd /build/flames-$flamesversion/bin
ln -s ../FLAMES-$flamescommit/match_cell_barcode .
cd /build/flames-$flamesversion/FLAMES-$flamescommit/python
ln -s ../FLAMES-$flamescommit/match_cell_barcode .

# add wrapper scripts

cd /build/flames-$flamesversion
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred

chmod ugo+x gtfToGenePred
chmod ugo+x gff3ToGenePred
mv gtfToGenePred bin
mv gff3ToGenePred bin

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
pythondir=`ls -d $dir/FLAMES-*/python`
PATH=$pythondir:$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
bulk_long_pipeline.py ${1+"$@"}
' > bulk_long_pipeline.py
chmod ugo+x bulk_long_pipeline.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
pythondir=`ls -d $dir/FLAMES-*/python`
PATH=$pythondir:$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
sc_long_pipeline.py ${1+"$@"}
' > sc_long_pipeline.py
chmod ugo+x sc_long_pipeline.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
pythondir=`ls -d $dir/FLAMES-*/python`
PATH=$pythondir:$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
if [ $# -eq 0 ] ; then
	echo "no subcommand given to flames: should be one of: match_cell_barcode `ls $pythondir | grep \.py\$ | tr "\n" " "`"
	exit 1
fi
command=$1
shift 1
$command ${1+"$@"}
' > flames
chmod ugo+x flames

# copy final appdir to genomecomb extra dir

rm /build/flames.tar.gz
cd /build
# tar cvzf flames-$flamesversion.tar.gz flames-$flamesversion
cp -ra flames-$flamesversion /io/extra$ARCH
cd /io/extra$ARCH/
ln -s flames-$flamesversion/flames .

conda deactivate

echo "Finished building flames"

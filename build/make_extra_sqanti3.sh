#!/bin/bash

# This script builds portable sqanti3 binaries using the Holy build box environment
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

#conda install -y -c conda-forge conda-pack

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r
conda init bash
. ~/.bash_profile

# sqanti3
# -------
cd /build

sqanti3version=4.2

conda create -y -n sqanti3
conda activate sqanti3

wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v4.2.tar.gz
tar -xvf v4.2.tar.gz
cd SQANTI3-4.2

# install dependencies
conda install -y cdna_cupcake
conda install -y ucsc-gtftogenepred

# install dependencies in SQANTI3.conda_env.yml
conda install -y samtools scipy star slamem bcbiogff bedtools 
conda install -y python>=3.7.6 biopython bioconductor-noiseq bx-python desalt gffread gmap kallisto minimap2 numpy openssl pandoc perl psutil pybedtools pysam
conda install -y r>=3.6.0 r-dplyr r-ggplot2 r-ggplotify r-gridbase r-gridextra r-htmltools r-reshape r-scales r-rmarkdown r-stringi r-dt r-plotly r-plyr
conda install -y ultra-bioinformatics

# make portable appdir
conda install -y conda-pack
cd /build
rm sqanti3.tar.gz || true
conda pack -n sqanti3 -o sqanti3.tar.gz
rm -rf sqanti3-$sqanti3version.old
mv sqanti3-$sqanti3version sqanti3-$sqanti3version.old || true
mkdir sqanti3-$sqanti3version
cd /build/sqanti3-$sqanti3version
tar xvzf ../sqanti3.tar.gz
cp -ra /build/SQANTI3-4.2 .

# replace sqanti3 gtfToGenePred with hbb compiled one (that should work on more systems)
cd /build/sqanti3-$sqanti3version
mv ./SQANTI3-4.2/utilities/gtfToGenePred ./SQANTI3-4.2/utilities/gtfToGenePred.sqanti3
cp -al ./bin/gtfToGenePred ./SQANTI3-4.2/utilities/gtfToGenePred


# make import from cdna_cupcake work
cd /build/sqanti3-$sqanti3version
git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd /build/sqanti3-$sqanti3version/cDNA_Cupcake
#pip install Cython
#python setup.py build
# clean a bit
rm -rf /build/sqanti3-$sqanti3version/cDNA_Cupcake/.git

# make excutables in appdir root that will take use the appdir env
cd /build/sqanti3-$sqanti3version

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
PYTHONPATH=$dir/cDNA_Cupcake/sequence $dir/SQANTI3-4.2/sqanti3_qc.py ${1+"$@"}
' > sqanti3_qc.py
chmod ugo+x sqanti3_qc.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
PYTHONPATH=$dir/cDNA_Cupcake/sequence $dir/SQANTI3-4.2/sqanti3_qc.py ${1+"$@"}
' > sqanti3_qc
chmod ugo+x sqanti3_qc

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
PYTHONPATH=$dir/cDNA_Cupcake/sequence $dir/SQANTI3-4.2/sqanti3_RulesFilter.py ${1+"$@"}
' > sqanti3_RulesFilter.py
chmod ugo+x sqanti3_RulesFilter.py

# R and Rscript do not come out of conda in state that is working in appdir
mv /build/sqanti3-$sqanti3version/bin/R /build/sqanti3-$sqanti3version/bin/R.conda || true
mv /build/sqanti3-$sqanti3version/bin/Rscript /build/sqanti3-$sqanti3version/bin/Rscript.conda || true
cp /io/build/sqanti3_files/R* /build/sqanti3-$sqanti3version/bin

rm /build/sqanti3-$sqanti3version.tar.gz
cd /build
tar cvzf sqanti3-$sqanti3version.tar.gz sqanti3-$sqanti3version
cp -ra sqanti3-$sqanti3version /io/extra$ARCH
cd /io/extra$ARCH/

conda deactivate

echo "Finished building sqanti3"

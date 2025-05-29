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
source "${dir}/start_hbb3.sh"

# Parse arguments
# ===============

clair3version=1.1.1

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
yuminstall gcc-c++
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


# clair3
# -----

cd /build
mamba create -y -n clair3 -c bioconda -c conda-forge clair3=$clair3version

# version 1.1.0 did not want to install due to conda dependency problems, but this is solved in 1.1.1
# I am leaving the solution (I had a hard time finding) here (in comments) in case I need it again later
# problem: clair3 required samtools 1.15.1 and openssl 1.1.1w. Somehow a requirement for a more recent (than 1.15.1) htslib gets pulled in (don't know where) that requires openssl 3 leading to the conflict
# solution: install deps via conda, adding samtools from genomecomb, and install clair3 from source

#cd /build

# mamba deactivate
# mamba env remove -n clair3
# mamba search clair3
# mamba repoquery depends -c bioconda -c conda-forge clair3=1.1.1 --recursive --tree

#mamba create -y -n clair3 -c bioconda -c conda-forge \
#    openssl=1.1.1w \
#    xz=5.2.10 \
#    bzip2=1.0.8 \
#    cffi=1.14.4 \
#    libcurl=7.87.0 \
#    libgcc=13 \
#    libstdcxx=13 \
#    libzlib=1.2.13 \
#    numpy \
#    parallel \
#    pigz \
#    pypy3.6 \
#    pytables \
#    python=3.9 \
#    python_abi=3.9 \
#    tensorflow=2.8.0 \
#    whatshap=1.7 \
#    zstd
#mamba activate clair3
#cp /io/extern/samtools $CONDA_PREFIX/bin
#
#cd /build
#rm -rf Clair3-$clair3version || true
#wget https://github.com/HKU-BAL/Clair3/archive/refs/tags/v$clair3version.tar.gz
#tar xvzf v1.1.1.tar.gz
#rm v$clair3version.tar.gz
#cd /build/Clair3-$clair3version
#
#make PREFIX=${CONDA_PREFIX}
#cp -rv clair3 preprocess postprocess scripts shared $CONDA_PREFIX/bin
#cp clair3.py $CONDA_PREFIX/bin/
#cp run_clair3.sh $CONDA_PREFIX/bin/
#cp longphase libclair3* $CONDA_PREFIX/bin
#
#cd ${CONDA_PREFIX}/bin
#wget https://downloads.python.org/pypy/pypy3.10-v7.3.19-linux64.tar.bz2 && tar -jxvf ${CONDA_PREFIX}/bin/pypy3.10-v7.3.19-linux64.tar.bz2
#ln -sf pypy3.10-v7.3.19-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy3
#
#mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm clair3.tar.gz || true
conda pack -n clair3 -o clair3.tar.gz
rm -rf clair3-$clair3version-$arch.old || true
mv clair3-$clair3version-$arch clair3-$clair3version-$arch.old || true
mkdir /build/clair3-$clair3version-$arch
cd /build/clair3-$clair3version-$arch
tar xvzf ../clair3.tar.gz

# add models
mkdir /build/clair3-$clair3version-$arch/models || true
cd /build/clair3-$clair3version-$arch/models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz
tar xvzf clair3_models.tar.gz
rm clair3_models.tar.gz

## (extra) ONT models
mkdir /build/clair3-$clair3version-$arch/models || true

# git clone https://github.com/nanoporetech/rerio
cd /build

rm -rf rerio || true
git clone https://github.com/nanoporetech/rerio
cp -ral rerio/clair3_models/* /build/clair3-$clair3version-$arch/models

cd /build/clair3-$clair3version-$arch/models
for model in \
	r1041_e82_400bps_sup_v500.tar.gz \
	r1041_e82_400bps_hac_v500.tar.gz \
	r1041_e82_400bps_sup_v410.tar.gz \
	r1041_e82_400bps_hac_v410.tar.gz \
	r1041_e82_400bps_sup_v430.tar.gz \
	r1041_e82_400bps_hac_v430.tar.gz \
	r1041_e82_400bps_sup_v420.tar.gz \
	r1041_e82_400bps_hac_v420.tar.gz \
	r1041_e82_260bps_sup_v400.tar.gz \
	r1041_e82_260bps_hac_v400.tar.gz \
	r1041_e82_260bps_fast_g632.tar.gz \
	r1041_e82_260bps_sup_g632.tar.gz \
	r1041_e82_400bps_hac_g632.tar.gz \
	r1041_e82_400bps_fast_g615.tar.gz \
	r1041_e82_400bps_sup_g615.tar.gz \
	r1041_e82_400bps_hac_g615.tar.gz \
	r1041_e82_260bps_hac_g632.tar.gz \
	r1041_e82_400bps_fast_g632.tar.gz \
	r104_e81_sup_g5015.tar.gz \
	r104_e81_hac_g5015.tar.gz
do
	wget https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/$model
	tar xvzf $model
	rm $model
done

cd /build/clair3-$clair3version-$arch

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

cd /build
ln -sf clair3-$clair3version-$arch/clair3 .
ln -sf clair3-$clair3version-$arch/clair3 clair3-$clair3version
ln -sf clair3-$clair3version-$arch/run_clair3.sh .
rm clair3-$clair3version-$arch.tar.gz || true
tar cvzf clair3-$clair3version-$arch.tar.gz clair3-$clair3version-$arch clair3 run_clair3.sh
rm -rf /io/extra$ARCH/clair3-$clair3version-$arch
cp -ra clair3-$clair3version-$arch clair3-$clair3version clair3 run_clair3.sh /io/extra$ARCH

echo "Finished building clair3-$clair3version-$arch"

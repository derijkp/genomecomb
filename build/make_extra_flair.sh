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

flairversion=2.0
# build has patches, check if they still work/apply before just changing versions

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

# flair
# -----
cd /build

mamba create -y -n flair
mamba activate flair
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y flair=$flairversion python=3.9

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm flair.tar.gz || true
conda pack -n flair -o flair.tar.gz
rm -rf flair-$flairversion-$arch.old || true
mv flair-$flairversion-$arch flair-$flairversion-$arch.old || true
mkdir /build/flair-$flairversion-$arch
cd /build/flair-$flairversion-$arch
tar xvzf ../flair.tar.gz

# patches
# avoid collapse error on missing _ in names due to chrom not in junc_to_tn
mv lib/python3.9/site-packages/flair/identify_gene_isoform.py lib/python3.9/site-packages/flair/identify_gene_isoform.py.ori
cp /io/build/flair_files/identify_gene_isoform.patched.py lib/python3.9/site-packages/flair/identify_gene_isoform.py
# different color handling, not updated to new flair version
cp /io/build/flair_files/plot_isoform_usage.patched.py bin/plot_isoform_usage.patched

# make excutables in appdir root that will use the appdir env
cd /build/flair-$flairversion-$arch

mv bin/flair bin/flair.ori
cat << 'EOF' > bin/flair
#!/bin/sh
'''exec' python "$0" "$@"
' '''
# -*- coding: utf-8 -*-
import re
import sys
from flair.flair import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
chmod ugo+x bin/flair

mv bin/f2py bin/f2py.ori
cat << 'EOF' > bin/f2py
#!/bin/sh
'''exec' python "$0" "$@"
' '''
# -*- coding: utf-8 -*-
import re
import sys
from numpy.f2py.f2py2e import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
chmod ugo+x bin/f2py

mv bin/f2py3 bin/f2py3.ori
cp -al bin/f2py bin/f2py3

mv bin/f2py3.9 bin/f2py3.9.ori
cp -al bin/f2py bin/f2py3.9

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/flair ${1+"$@"}
' > flair.py
chmod ugo+x flair.py

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/flair ${1+"$@"}
' > flair
chmod ugo+x flair

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/bam2Bed12 ${1+"$@"}
' > bam2Bed12
chmod ugo+x bam2Bed12

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/plot_isoform_usage ${1+"$@"}
' > plot_isoform_usage
chmod ugo+x plot_isoform_usage

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/predictProductivity ${1+"$@"}
' > predictProductivity
chmod ugo+x predictProductivity

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/mark_intron_retention ${1+"$@"}
' > mark_intron_retention
chmod ugo+x mark_intron_retention

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/diff_iso_usage ${1+"$@"}
' > diff_iso_usage
chmod ugo+x diff_iso_usage

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/bin/diffsplice_fishers_exact ${1+"$@"}
' > diffsplice_fishers_exact
chmod ugo+x diffsplice_fishers_exact

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
python $dir/lib/python3.9/site-packages/flair/identify_gene_isoform.py ${1+"$@"}
' > identify_gene_isoform.py
chmod ugo+x identify_gene_isoform.py

# package
cd /build
rm flair.tar.gz
ln -sf flair-$flairversion-$arch/flair.py flair.py
ln -sf flair-$flairversion-$arch/flair flair
ln -sf flair-$flairversion-$arch/flair flair-$flairversion
tar cvzf flair-$flairversion-$arch.tar.gz flair-$flairversion-$arch flair.py flair
cp -ra flair-$flairversion-$arch flair.py flair /io/extra$ARCH
cd /io/extra$ARCH/

conda deactivate

echo "Finished building flair"

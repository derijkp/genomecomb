#!/bin/bash

# This script builds portable crispat binaries using the Holy build box environment
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

crispatversion=0.9.8

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
# source /hbb_exe/activate

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

# crispat
# -------
cd /build

mamba create -y -n crispat
mamba activate crispat
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba install -y python=3.10
mamba install -y openBLAS
mamba install -y cmake xsimd pkgconfig
pip install "crispat==$crispatversion"

mamba deactivate

# make package
# ------------

cd /build
# installing conda-pack in the beginning causes further commands to fail (network/ssl), so we do it here at the end
mamba install -y -c conda-forge conda-pack

rm crispat.tar.gz || true
conda pack -n crispat -o crispat.tar.gz
rm -rf crispat-$crispatversion-$arch.old || true
mv crispat-$crispatversion-$arch crispat-$crispatversion-$arch.old || true
mkdir /build/crispat-$crispatversion-$arch
cd /build/crispat-$crispatversion-$arch
tar xvzf ../crispat.tar.gz

cd /build/crispat-$crispatversion-$arch

cat << 'EOF' > bin/crispat
#!/bin/sh
'''exec' python "$0" "$@"
' '''
import argparse
import crispat
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns 
def run_crispat(input_csv, output_dir, method):
    """Run CRISPRat analysis on the given CSV file using the specified method."""
    # Check if the method is available
    if not hasattr(crispat, method):
        print(f"Error: Method '{method}' not found in CRISPRat.")
        return
    # Check the input CSV
    # data = pd.read_csv(input_csv)
    # data.iloc[0:5, 0:5]
    # Create anndata object
    output_dir = output_dir +'/'
    crispat.create_anndata_from_csv(input_csv, output_dir)
    ad.read_h5ad(output_dir + '/gRNA_counts.h5ad')
    # Get the method and apply it
    print('Running crispat using method ',method)
    if method == "cellranger" or method == "gauss":
        crispat.ga_gauss(output_dir + 'gRNA_counts.h5ad', output_dir + 'guide_assignments/gauss/', inference = 'em')
    elif method == "cellrangernonzero" or method == "gaussnonzero":
        crispat.ga_gauss(output_dir + 'gRNA_counts.h5ad', output_dir + 'guide_assignments/gaussnonzero/', nonzero = True, inference = 'em')
    elif method == "poisson_gauss":
        crispat.ga_poisson_gauss(output_dir + 'gRNA_counts.h5ad', output_dir + 'guide_assignments/poisson_gauss/')
    elif method == "ratios":
        print('Writing', output_dir + 'guide_assignments/ratios/')
        crispat.ga_ratio(output_dir + 'gRNA_counts.h5ad', [0.3, 0.5, 0.7], output_dir + 'guide_assignments/ratios/')
    else:
        print("Unknown crispat method: ", method)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CRISPRat analysis from command line.")
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory.")
    parser.add_argument("-m", "--method", required=True, help="CRISPRat analysis method to use.")
    args = parser.parse_args()
    run_crispat(args.input, args.output, args.method)
EOF
chmod ugo+x bin/crispat

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/crispat ${1+"$@"}
' > crispat
chmod ugo+x crispat

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/bin/python ${1+"$@"}
' > python
chmod ugo+x python

cd /build
ln -sf crispat-$crispatversion-$arch/crispat .
ln -sf crispat-$crispatversion-$arch/crispat crispat-$crispatversion
rm crispat-$crispatversion-$arch.tar.gz || true
tar cvzf crispat-$crispatversion-$arch.tar.gz crispat-$crispatversion-$arch crispat
rm -rf /io/extra$ARCH/crispat-$crispatversion-$arch
cp -ra crispat-$crispatversion-$arch crispat-$crispatversion crispat /io/extra$ARCH

echo "Finished building crispat-$crispatversion-$arch"

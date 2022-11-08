#!/bin/bash

# This script creates a portable nextflow appdir
# While it does not compile nextflow, Holy build box is still used to have same environment as other tools
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

all=1
extra=1
while [[ "$#" -gt 0 ]]; do case $1 in
	-a|-all|--all) all="$2"; shift;;
	-e|-extra|--extra) extra="$2"; shift;;
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

if [ $arch != "linux-x86_64" ]; then
	echo "Only linux-x86_64 supported"
	exit 1
fi

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
#yuminstall git
yuminstall wget
#yuminstall gcc-c++
#yuminstall centos-release-scl
sudo yum upgrade -y
yuminstall devtoolset-11
## use source instead of scl enable so it can run in a script
## scl enable devtoolset-11 bash
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

# nextflow
# --------
cd /build

nextflowversion=22.10.0

mkdir /build/nextflow-$nextflowversion-$arch
cd /build/nextflow-$nextflowversion-$arch


# java 18
curl -O https://download.java.net/java/GA/jdk18/43f95e8614114aeaa8e8a5fcf20a682d/36/GPL/openjdk-18_linux-x64_bin.tar.gz
tar xvf openjdk-18_linux-x64_bin.tar.gz
rm openjdk-18_linux-x64_bin.tar.gz

wget https://github.com/nextflow-io/nextflow/releases/download/v$nextflowversion/nextflow-$nextflowversion-all
chmod ugo+x nextflow-$nextflowversion-all

wget https://github.com/nextflow-io/nextflow/archive/refs/tags/v$nextflowversion.tar.gz
tar xvzf v$nextflowversion.tar.gz
cp ./nextflow-22.10.0/README.md .

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
export JAVA_HOME=$dir/jdk-18/
export PATH=$dir/bin:$dir/jdk-18/bin:$PATH
export NXF_JAVA_HOME=$dir/jdk-18/
$dir/nextflow-22.10.0-all ${1+"$@"}
' > nextflow
chmod ugo+x nextflow

# make package
# ------------

cd /build
tar cvzf nextflow-$nextflowversion-$arch.tar.gz nextflow-$nextflowversion-$arch
cp -ra nextflow-$nextflowversion-$arch /io/extra$ARCH
cd /io/extra$ARCH/
rm nextflow || true
ln -s nextflow-$nextflowversion-$arch/nextflow .

echo "Finished building nextflow"

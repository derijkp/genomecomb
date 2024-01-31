#!/bin/bash

# This script builds portable gatk3 binaries using the Holy build box environment
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

gatk3version=3.8-1-0

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
yuminstall wget
yuminstall java

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

# gatk3
# -----

# make package
# ------------

cd /build
mkdir /build/gatk3-$gatk3version-$arch
cd /build/gatk3-$gatk3version-$arch

cp -raL /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.275.b01-0.el6_10.x86_64 /build/gatk3-$gatk3version-$arch
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
tar xvjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
rm GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
$dir/java-1.8.0-openjdk-1.8.0.275.b01-0.el6_10.x86_64/jre/bin/java -XX:ParallelGCThreads=1 -Xms512m -Xmx8g -jar $dir/GenomeAnalysisTK-*/GenomeAnalysisTK.jar ${1+"$@"}
' > gatk3
chmod ugo+x gatk3

cd /build
ln -sf gatk3-$gatk3version-$arch/gatk3 .
ln -sf gatk3-$gatk3version-$arch/gatk3 gatk3-$gatk3version
rm gatk3-$gatk3version-$arch.tar.gz || true
tar cvzf gatk3-$gatk3version-$arch.tar.gz gatk3-$gatk3version-$arch gatk3 gatk3-$gatk3version
rm -rf /io/extra$ARCH/gatk3-$gatk3version-$arch
cp -ra gatk3-$gatk3version-$arch gatk3 gatk3-$gatk3version /io/extra$ARCH
cd /io/extra$ARCH/

echo "Finished building gatk3-$gatk3version-$arch"

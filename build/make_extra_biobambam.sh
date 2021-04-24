#!/bin/bash

# This script builds portable biobambam binaries using the Holy build box environment
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)
# -a|-all|--all: if 1 (default) all binaries are (re)build, if 0, only the ones missing are build

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
while [[ "$#" -gt 0 ]]; do case $1 in
	-a|-all|--all) all="$2"; shift;;
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
#yuminstall xz

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

# libmaus (used by biobambam)
# ---------------------------
if [ $all = 1 ] || [ ! -f /build/lib/libmaus2.a ] ; then
	# also a library, needs -fPIC, so compile as lib
	source /hbb_shlib/activate
	# download https://github.com/gt1/libmaus2/archive/2.0.349-release-20170704123913.tar.gz libmaus2-2.0.349-release-20170704123913.tar.gz
	download https://gitlab.com/german.tischler/libmaus2/-/archive/2.0.611-release-20190408112550/libmaus2-2.0.611-release-20190408112550.tar.gz libmaus2-2.0.611-release-20190408112550.tar.gz
	cd /build/libmaus2-2.0.611-release-20190408112550
	make distclean || true
	env CFLAGS="$STATICLIB_CFLAGS -I/build/include" LDFLAGS="$LDFLAGS -L/build/lib" \
		./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	source /hbb_exe/activate
fi

# biobambam
# ---------
if [ $all = 1 ] || [ ! -f /io/extra$ARCH/bamsort ] ; then
	# download https://github.com/gt1/biobambam2/archive/2.0.73-release-20170620145717.tar.gz biobambam2-2.0.73-release-20170620145717.tar.gz
	download https://gitlab.com/german.tischler/biobambam2/-/archive/2.0.95-release-20190320141403/biobambam2-2.0.95-release-20190320141403.tar.gz biobambam2-2.0.95-release-20190320141403.tar.gz
	cd /build/biobambam2-2.0.95-release-20190320141403
	make distclean || true
	./configure --prefix=/build --with-libmaus2=/build
	make
	dest=/io/extra$ARCH/biobambam2
	mkdir $dest
	cp src/bammarkduplicates2 src/bamsort src/bamtofastq $dest
	strip $dest/bammarkduplicates2 $dest/bamsort $dest/bamtofastq
	cp README.md COPYING $dest
	cd /io/extra
	ln -sf biobambam2/bam* .
	
fi

echo "Finished building extra external binaries in extra$ARCH/biobambam2"

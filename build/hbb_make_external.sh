#!/bin/bash

# This script builds the genomecomb binaries using the Holy build box environment
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
    else
        tar xvzf $filename
    fi
}

# zstd (for lib)
# --------------
cd /build
if [ ! -f "/build/zstd-1.3.8/lib/libzstd.a" ] ; then
	source /hbb_shlib/activate
	wget -c https://github.com/facebook/zstd/releases/download/v1.3.8/zstd-1.3.8.tar.gz
	tar xvzf zstd-1.3.8.tar.gz
	cd /build/zstd-1.3.8
	make
	source /hbb_exe/activate
fi

# zstdmt
# ------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/zstd-mt ] ; then
	download https://github.com/mcmilk/zstdmt/archive/0.6.tar.gz zstdmt-0.6.tar.gz
	cd /build/zstdmt-0.6/programs
	make clean
	make loadsource
	make zstd-mt
	rm -f /io/extern$ARCH/zstd-mt
	cp zstd-mt /io/extern$ARCH
	strip /io/extern$ARCH/zstd-mt
fi

# lz4
# ---
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/lz4 ] ; then
	download https://github.com/lz4/lz4/archive/v1.8.3.tar.gz lz4-1.8.3.tar.gz
	cd /build/lz4-1.8.3
	make clean
	make
	cp lz4 /io/extern$ARCH
	strip /io/extern$ARCH/lz4
	mkdir /io/extern$ARCH/lz4_docs
	cp LICENSE /io/extern$ARCH/lz4_docs/lz4.LICENSE
fi

# razip
# -----
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/razip ] ; then
	cd /build
	tar xvzf /io/extern-src/razip.tar.gz
	cd /build/razip
	make clean
	make
	cp razip /io/extern$ARCH
	strip /io/extern$ARCH/razip
fi

# gnusort8
# --------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/gnusort8 ] ; then
	download https://ftp.gnu.org/gnu/coreutils/coreutils-8.9.tar.gz
	cd /build/coreutils-8.9
	if [ ! -f src/sort.c.orig ] ; then
		cp src/sort.c src/sort.c.orig
	fi
	cp /io/extern-src/gnusort8_sort.c src/sort.c
	make distclean
	./configure
	make
	cp src/sort /io/extern$ARCH/gnusort8
	strip /io/extern$ARCH/gnusort8
fi

# pigz
# ----
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/pigz ] ; then
	download https://github.com/madler/pigz/archive/v2.4.tar.gz pigz-2.4.tar.gz
	cd /build/pigz-2.4
	if [ ! -f Makefile.ori ] ; then
		cp Makefile Makefile.ori
	fi
	grep -v 'LDFLAGS=' Makefile.ori > Makefile
	make clean
	make
	cp pigz /io/extern$ARCH
	strip /io/extern$ARCH/pigz
fi

# bzlib, used by htslib
# ---------------------
if [ ! -f "/build/lib/libbz2.a" ] ; then
	source /hbb_shlib/activate
	download https://sourceforge.net/projects/bzip2/files/latest/bzip2-1.0.6.tar.gz
	cd /build/bzip2-1.0.6
	if [ ! -f Makefile.ori ] ; then
		cp Makefile Makefile.ori
	fi
	grep -v 'LDFLAGS=' Makefile.ori | grep -v 'CFLAGS=' > Makefile
	make clean
	env CFLAGS="$STATICLIB_CFLAGS -D_FILE_OFFSET_BITS=64" CXXFLAGS="$STATICLIB_CXXFLAGS -D_FILE_OFFSET_BITS=64" \
		make
	cp /build/bzip2-1.0.6/libbz2.a /build/lib
	cp /build/bzip2-1.0.6/*.h /build/include
	source /hbb_exe/activate
fi

# xz/lzma, used by htslib
# -----------------------
if [ ! -f "/build/lib/liblzma.a" ] ; then
	source /hbb_shlib/activate
	download https://tukaani.org/xz/xz-5.2.3.tar.gz
	cd /build/xz-5.2.3
	make distclean
	env CFLAGS="$STATICLIB_CFLAGS" CXXFLAGS="$STATICLIB_CXXFLAGS" \
  		./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	source /hbb_exe/activate
fi

# ncurses, used by samtools
# -------------------------
if [ ! -f "/build/lib/liblzma.a" ] ; then
	source /hbb_shlib/activate
	download http://ftp.gnu.org/gnu/ncurses/ncurses-6.1.tar.gz
	cd /build/ncurses-6.1
	make distclean
	env CFLAGS="$STATICLIB_CFLAGS" CXXFLAGS="$STATICLIB_CXXFLAGS" \
  		./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	source /hbb_exe/activate
fi

# htslib (tabix, bgzip, lib for samtools)
# ---------------------------------------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/tabix ] || [ ! -f /io/extern$ARCH/bgzip ] || [ ! -f /build/lib/libhts.a ] ; then
	htsversion=1.9
	# also a library, needs -fPIC, so compile as lib
	source /hbb_shlib/activate
	download https://github.com/samtools/htslib/releases/download/$htsversion/htslib-$htsversion.tar.bz2
	cd /build/htslib-$htsversion
	make distclean
	env CFLAGS="$STATICLIB_CFLAGS -I/build/include" LDFLAGS="$LDFLAGS -L/build/lib" \
		./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	rm /build/lib/libhts.so*
	cp /build/bin/tabix /io/extern$ARCH/tabix
	cp /build/bin/bgzip /io/extern$ARCH/bgzip
	strip /io/extern$ARCH/tabix
	strip /io/extern$ARCH/bgzip
	source /hbb_exe/activate
fi

# samtools
# --------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/samtools ] ; then
	samversion=1.9
	download https://github.com/samtools/samtools/releases/download/$samversion/samtools-$samversion.tar.bz2
	cd /build/samtools-$samversion
	make distclean
	./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	cp /build/bin/samtools /io/extern$ARCH
	strip /io/extern$ARCH/samtools
fi

# bcftools
# --------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/bcftools ] ; then
	bcfversion=1.9
	download https://github.com/samtools/bcftools/releases/download/$bcfversion/bcftools-$bcfversion.tar.bz2
	cd /build/bcftools-$bcfversion
	make distclean
	./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	cp /build/bin/bcftools /io/extern$ARCH
	strip /io/extern$ARCH/bcftools
fi

# bwa
# ---
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/bwa ] ; then
	bwaversion=0.7.17
	download https://sourceforge.net/projects/bio-bwa/files/bwa-$bwaversion.tar.bz2
	cd /build/bwa-$bwaversion
	make
	cp /build/bin/bwa /io/extern$ARCH
	strip /io/extern$ARCH/bwa
fi

# libmaus (used by biobambam)
# ---------------------------
if [ $all = 1 ] || [ ! -f /build/lib/libmaus2.a ] ; then
	# also a library, needs -fPIC, so compile as lib
	source /hbb_shlib/activate
	download https://github.com/gt1/libmaus2/archive/2.0.349-release-20170704123913.tar.gz libmaus2-2.0.349-release-20170704123913.tar.gz
	cd /build/libmaus2-2.0.349-release-20170704123913
	make distclean
	env CFLAGS="$STATICLIB_CFLAGS -I/build/include" LDFLAGS="$LDFLAGS -L/build/lib" \
		./configure --prefix=/build --disable-shared --enable-static
	make
	make install
	source /hbb_exe/activate
fi

# biobambam
# ---------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH/bamsort ] ; then
	download https://github.com/gt1/biobambam2/archive/2.0.73-release-20170620145717.tar.gz biobambam2-2.0.73-release-20170620145717.tar.gz
	cd /build/biobambam2-2.0.73-release-20170620145717
	make distclean
	./configure --prefix=/build --with-libmaus2=/build
	make
	dest=/io/extern$ARCH/biobambam2
	mkdir $dest
	cp src/bammarkduplicates2 src/bamsort src/bamtofastq $dest
	strip $dest/bammarkduplicates2 $dest/bamsort $dest/bamtofastq
	cp README.md COPYING $dest
	cd /io/extern
	ln -sf biobambam2/bam* .
	
fi

# ea-utils
# --------
if [ $all = 1 ] || [ ! -f /io/extern$ARCH ] ; then
	cd /build
	wget -c https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/ea-utils/ea-utils.1.1.2-537.tar.gz
	tar xvzf ea-utils.1.1.2-537.tar.gz
	cd /build/ea-utils.1.1.2-537
	if [ ! -f Makefile.ori ] ; then 
		cp Makefile Makefile.ori
	fi
	if [ ! -f fastq-lib.cpp.ori ] ; then 
		cp fastq-lib.cpp fastq-lib.cpp.ori
	fi
	if [ ! -f fastq-mcf.c.ori ] ; then 
		cp fastq-mcf.c fastq-mcf.c.ori
	fi
	cp /io/extern-src/ea-utils-changes/*.c .
	cp /io/extern-src/ea-utils-changes/Makefile .
	make clean
	make fastq-mcf
	make fastq-stats
	dest=/io/extern$ARCH/ea-utils
	mkdir $dest || true
	cp fastq-mcf fastq-stats $dest
	strip $dest/fastq-mcf $dest/fastq-stats
	cp /io/extern-src/ea-utils-changes/README $dest
	cd /io/extern
	ln -sf ea-utils/fastq-* .
fi

echo "Finished building external binaries"

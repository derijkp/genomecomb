#!/bin/bash

# This script builds a distributable directory contained version of R using the Holy build box environment
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
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# X software does not compile with these settings
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
yuminstall libX11-devel
yuminstall libXt-devel
yuminstall tcl
#yuminstall compat-gcc-34-g77
yuminstall gcc-gfortran
yuminstall gcc-c++
yuminstall readline-devel
yuminstall xz
#yuminstall java
#yuminstall texlive
yuminstall texlive-latex

# force static nuilds by removing so for libpng, libtiff
# sudo rm -f /usr/lib64/libpng*

## need more recent curl
#sudo echo '[CityFan]
#name=City Fan Repo
#baseurl=http://www.city-fan.org/ftp/contrib/yum-repo/rhel$releasever/$basearch/
#enabled=1
#gpgcheck=0' > /tmp/temp
#sudo mv -f /tmp/temp /etc/yum.repos.d/city-fan.repo
#sudo yum clean all
# sudo yum install -y libcurl
# sudo yum install -y binutils

#for dir in lib include bin share ; do
#	echo $dir
#	mkdir /build/$dir || true
#	sudo rmdir /usr/local/$dir || true
#	sudo rm /usr/local/$dir || true
#	sudo ln -s /build/$dir /usr/local/$dir
#done

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

# zlib
# ----
cd /build
download https://zlib.net/zlib-1.2.11.tar.gz
cd /build/zlib-1.2.11
make distclean || true
./configure -static
make CFLAGS="-O3 -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -fPIC"
sudo make install
sudo rm /usr/local/lib/libz.so*

# bzip2
# -----
cd /build
download https://sourceforge.net/projects/bzip2/files/latest/bzip2-1.0.6.tar.gz
cd /build/bzip2-1.0.6
make clean || true
make CFLAGS="-fPIC"
sudo make install

# lzma
# ----
cd /build
download https://tukaani.org/xz/xz-5.2.3.tar.gz
cd /build/xz-5.2.3
make distclean || true
./configure -enable-static -disable-shared
make CFLAGS="-fPIC"
sudo make install

# pcre
# ----
cd /build
download ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.43.tar.gz
cd /build/pcre-8.43
make distclean || true
./configure --enable-utf -enable-static -disable-shared
make CFLAGS="-fPIC"
sudo make install
#sudo rm /usr/local/lib/libpcre.so*

# openssl
# -------
cd /build
download https://www.openssl.org/source/openssl-1.1.1b.tar.gz
cd /build/openssl-1.1.1b
make clean || true
make distclean || true
./config -enable-static
make CFLAGS="-fPIC"
sudo make install
#sudo rm /usr/local/lib64/libssl.so* /usr/local/lib64/libcrypto.so*

# libcurl
# -------
cd /build
download https://curl.haxx.se/download/curl-7.22.0.tar.gz
cd /build/curl-7.22.0
make distclean || true

LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    CPPFLAGS="-I/usr/local/include -I/usr/local/include/openssl" \
    LIBS="-ldl" \
    ./configure -with-ssl -enable-static -disable-shared

PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure -with-ssl=/usr/local -enable-static -disable-shared

PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    ./configure -with-ssl=/usr/local -enable-static -disable-shared

env PKG_CONFIG_PATH=/usr/local/lib64/pkgconfig ./configure --with-ssl -enable-static -disable-shared

PATH=/usr/local/bin $PATH PKG_CONFIG_PATH=/usr/local/lib64/pkgconfig LIBS="-ldl" ./configure -with-ssl -enable-static -disable-shared

LIBS="-ldl" ./configure --with-ssl=/usr/local/ssl -enable-static -disable-shared

make CFLAGS="-fPIC"
sudo make install
#sudo rm /usr/local/lib/libcurl.so*

# libpng
# ------
cd /build
download https://sourceforge.net/projects/libpng/files/libpng12/1.2.59/libpng-1.2.59.tar.gz/download
cd /build/libpng-1.2.59
make distclean || true
LDFLAGS="-L/usr/local/lib" CPPFLAGS="-I/usr/local/include" ./configure -enable-static -disable-shared
make CFLAGS="-fPIC"
sudo make install

# libtiff
# -------
cd /build
download https://download.osgeo.org/libtiff/tiff-4.0.10.tar.gz
cd /build/tiff-4.0.10
make distclean || true
./configure -enable-static -disable-shared
make CFLAGS="-fPIC"
sudo make install

# pixman
# -------
cd /build
download https://www.cairographics.org/releases/pixman-0.36.0.tar.gz
cd /build/pixman-0.36.0
make distclean || true
LDFLAGS="-L/usr/local/lib" CPPFLAGS="-I/usr/local/include" ./configure --enable-static --disable-shared
make CFLAGS="-fPIC"
sudo make install

# cairo
# -------
cd /build
download https://cairographics.org/snapshots/cairo-1.9.14.tar.gz
cd /build/cairo-1.9.14
make distclean || true
PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure --enable-static --disable-shared
make CFLAGS="-fPIC"
sudo make install

# R
# -
version=3.5.3
download https://cran.r-project.org/src/base/R-3/R-$version.tar.gz
cd /build/R-$version
make distclean || true
# make sure the right curl lib is used
PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    ./configure --disable-java --enable-R-shlib --prefix=/build/diR-$version-$arch
PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    make
make install

echo '
puts [list set argv $argv]
foreach {file basedir locate} $argv break
if {![file exists $file.ori]} {file copy $file $file.ori}
set f [open $file.ori]
set c [read $f]
close $f
regsub {# Shell wrapper for R executable.} $c {# Shell wrapper for R executable.
script="$(readlink -f "$0")"
R_BASEDIR="$(dirname "$(dirname "$script")")"} c
regsub -all $basedir $c {${R_BASEDIR}} c
if {$locate ne ""} {
	set pattern {$(dirname "$(dirname "$script")")}
	set pos [string first $pattern $c]
	set c [string range $c 0 [expr {$pos-1}]]$locate[string range $c [expr {$pos+[string length $pattern]}] end]
}
set o [open $file w]
puts $o $c
close $f
exec chmod ugo+x $file
' > /tmp/convert
tclsh /tmp/convert /build/diR-$version-$arch/bin/R /build/diR-$version-$arch
tclsh /tmp/convert /build/diR-$version-$arch/lib64/R/bin/R /build/diR-$version-$arch '$(dirname "$(dirname "$(dirname "$(dirname "$script")")")")'

echo "Finished building diR"

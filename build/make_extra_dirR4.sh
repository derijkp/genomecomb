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
source "${dir}/start_hbb3.sh"

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
yuminstall libXft-devel
yuminstall pango-devel
yuminstall libjpeg-devel
yuminstall tcl
yuminstall tcl-devel
yuminstall tk-devel
#yuminstall compat-gcc-34-g77
yuminstall gcc-gfortran
yuminstall gcc-c++
yuminstall readline-devel
yuminstall xz
yuminstall libtool
yuminstall pcre2-devel
#yuminstall java
#yuminstall texlive
yuminstall texlive-latex
yuminstall devtoolset-11
## use source instead of scl enable so it can run in a script
## scl enable devtoolset-11 bash
source /opt/rh/devtoolset-11/enable
yuminstall xz

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
        wget --no-check-certificate -c -O $filename $url
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

export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
export PATH=/usr/local/bin:$PATH

## zlib
## ----
#cd /build
#download https://zlib.net/zlib-1.2.11.tar.gz
#cd /build/zlib-1.2.11
#make distclean || true
#./configure -static
#make CFLAGS="-O3 -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -fPIC"
#sudo make install

# bzip2
# -----
cd /build
download https://sourceforge.net/projects/bzip2/files/latest/bzip2-1.0.6.tar.gz
cd /build/bzip2-1.0.6
make clean || true
make CFLAGS="-fPIC"
sudo make install

## lzma
## ----
#cd /build
#download https://tukaani.org/xz/xz-5.2.3.tar.gz
#cd /build/xz-5.2.3
#make distclean || true
#./configure -enable-static -disable-shared
#make CFLAGS="-fPIC"
#sudo make install
#
## pcre
## ----
#cd /build
#download ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.43.tar.gz
#cd /build/pcre-8.43
#make distclean || true
#./configure --enable-utf -enable-static -disable-shared
#make CFLAGS="-fPIC"
#sudo make install
##sudo rm /usr/local/lib/libpcre.so*
#
## pcre2
## ----
#cd /build
#download https://github.com/PhilipHazel/pcre2/releases/download/pcre2-10.39/pcre2-10.39.tar.gz
#cd /build/pcre2-10.39
#make distclean || true
#./configure --enable-utf -enable-static -disable-shared
#make CFLAGS="-fPIC"
#sudo make install

# openssl
# -------
cd /build
download https://www.openssl.org/source/openssl-1.1.1m.tar.gz
cd /build/openssl-1.1.1m
make clean || true
make distclean || true
./config -enable-static
make CFLAGS="-fPIC"
sudo make install
#sudo rm /usr/local/lib64/libssl.so* /usr/local/lib64/libcrypto.so*

# libcurl
# -------
cd /build
# download https://curl.haxx.se/download/curl-7.22.0.tar.gz
download https://curl.se/download/curl-7.80.0.tar.gz
cd /build/curl-7.80.0
make distclean || true

# PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure -with-ssl=/usr/local -enable-static -disable-shared
PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure -with-ssl=/usr/local -enable-static
make CFLAGS="-fPIC"
sudo make install
#sudo rm /usr/local/lib/libcurl.so*

## libpng
## ------
#cd /build
#download https://sourceforge.net/projects/libpng/files/libpng16/1.6.37/libpng-1.6.37.tar.gz
#cd /build/libpng-1.6.37
#make distclean || true
#LDFLAGS="-L/usr/local/lib" CPPFLAGS="-I/usr/local/include" ./configure -enable-static
#make CFLAGS="-fPIC"
#sudo make install
#
## libtiff
## -------
#cd /build
#download https://download.osgeo.org/libtiff/tiff-4.0.10.tar.gz
#cd /build/tiff-4.0.10
#make distclean || true
#./configure -enable-static
#make CFLAGS="-fPIC"
#sudo make install
#
## pixman
## -------
#cd /build
#download https://www.cairographics.org/releases/pixman-0.36.0.tar.gz
#cd /build/pixman-0.36.0
#make distclean || true
#LDFLAGS="-L/usr/local/lib" CPPFLAGS="-I/usr/local/include" ./configure --enable-static
#make CFLAGS="-fPIC"
#sudo make install
#
## cairo
## -----
#cd /build
#download https://cairographics.org/snapshots/cairo-1.9.14.tar.gz
#cd /build/cairo-1.9.14
#make distclean || true
#PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure --enable-static --without-python
#make CFLAGS="-fPIC"
#sudo make install

# libxml2
# -------
cd /build
download ftp://xmlsoft.org/libxml2/libxml2-2.9.9.tar.gz
cd /build/libxml2-2.9.9

make distclean || true
PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure --enable-static --without-python
make CFLAGS="-fPIC"
sudo make install

## libjpeg
## -------
#cd /build
#download https://sourceforge.net/projects/libjpeg/files/libjpeg/6b/jpegsrc.v6b.tar.gz
#cd /build/jpeg-6b
#
#make distclean || true
#PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure
#sed -i s#./libtool#libtool#g Makefile
#make CFLAGS="-fPIC"
#sudo make install

# R
# -
# useful hints in getting it compiled found at https://tdhock.github.io/blog/2017/compiling-R/
#version=3.5.3
#majorversion=3
version=4.1.2
majorversion=4

cd /build
rm -rf /build/R-$version
download https://cran.r-project.org/src/base/R-$majorversion/R-$version.tar.gz
cd /build/R-$version
make distclean || true

# make sure the manually installed libs (in /usr/local) are used
cd /build/R-$version
PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    CPPFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    ./configure  --enable-static --disable-java --enable-R-shlib --prefix=/build/dirR-$version-$arch

PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    CPPFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    make

PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    CPPFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    make install

# copy dl libraries to R local lib dir
cp -a -f /usr/lib64/libgfortran*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libgomp*.so.1* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/*pango* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /lib64/libreadline*.so.6* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libquadmath.so.0* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /build/zlib-1.2.11/libz.so.1* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/local/lib64/libssl*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/local/lib64/libcrypto*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/local/lib/libcurl*.so* /build/dirR-$version-$arch/lib64/R/lib
# cp -a -f /usr/local/lib/libpng*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libpng*.so* /build/dirR-$version-$arch/lib64/R/lib
# cp -a -f /usr/local/lib/libtiff*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libtiff*.so* /build/dirR-$version-$arch/lib64/R/lib
#cp -a -f /usr/local/lib/libpixman*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libpixman*.so* /build/dirR-$version-$arch/lib64/R/lib
#cp -a -f /usr/local/lib/libcairo*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libcairo*.so* /build/dirR-$version-$arch/lib64/R/lib
#cp -a -f /usr/local/lib/libxml2*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libxml2*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libjpeg*.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tcl* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/dirR-$version-$arch/lib64/R/lib
cp -a -fr /usr/share/tcl8.5 /build/dirR-$version-$arch/lib64/R/share
cp -a -fr /usr/share/tk8.5 /build/dirR-$version-$arch/lib64/R/share
#cp -a -f /usr/lib64/libpng12.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/liblzma* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicuuc.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicui18n.so* /build/dirR-$version-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicudata.so* /build/dirR-$version-$arch/lib64/R/lib

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

tclsh /tmp/convert /build/dirR-$version-$arch/bin/R /build/dirR-$version-$arch
cp /build/dirR-$version-$arch/bin/R.ori /build/dirR-$version-$arch/R.ori
tclsh /tmp/convert /build/dirR-$version-$arch/R /build/dirR-$version-$arch '$(dirname "$script")'
tclsh /tmp/convert /build/dirR-$version-$arch/lib64/R/bin/R /build/dirR-$version-$arch '$(dirname "$(dirname "$(dirname "$(dirname "$script")")")")'

# packages
# --------

#if [ ! -f /build/dirR-$version-$arch/lib64/R/etc/Makeconf.ori ] ; then
#	cp /build/dirR-$version-$arch/lib64/R/etc/Makeconf /build/dirR-$version-$arch/lib64/R/etc/Makeconf.ori
#fi
#sed 's/CXX11 =/CXX11 = g++ -std=c++11/g' /build/dirR-$version-$arch/lib64/R/etc/Makeconf.ori \
#    | sed 's#CXX11FLAGS =  $(LTO)#CXX11FLAGS = -I/usr/local/include -fPIC#g' \
#    > /build/dirR-$version-$arch/lib64/R/etc/Makeconf

/build/dirR-$version-$arch/R --vanilla -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("rmarkdown", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("shiny", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("googleVis", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("DT", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("pheatmap", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("curl", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("plotly", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("jpeg", repos="http://cran.us.r-project.org")'

/build/dirR-$version-$arch/R --vanilla -e 'install.packages("xml2", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("XML", repos="http://cran.us.r-project.org")'

/build/dirR-$version-$arch/R --vanilla -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("DESeq2")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("PCAtools")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("pcaExplorer")'

/build/dirR-$version-$arch/R --vanilla -e $'BiocManager::install("Rhtslib")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("Rsamtools")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("ggbio")'

/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("bambu")'

# leafcutter
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("rstan", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("rstantools", repos="http://cran.us.r-project.org")'
yuminstall gsl-devel
cp -a /usr/lib64/libgsl*.so* /build/dirR-$version-$arch/lib64/R/lib
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("DirichletMultinomial")'
/build/dirR-$version-$arch/R --vanilla -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

/build/dirR-$version-$arch/R --vanilla -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")'

# other
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("MOFA2")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("mixOmics")'
# Rhtslib does not compile for some reason, so not for now
# /build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("chromVAR")'
/build/dirR-$version-$arch/R --vanilla -e 'BiocManager::install("OmicCircos")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("RCircos", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("ggdendro", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("GGally", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("ComplexUpset", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("corrplot", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("eulerr", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("vioplot", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("gplots", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("gtools", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("ggrepel", repos="http://cran.us.r-project.org")'
# /build/dirR-$version-$arch/R --vanilla -e 'install.packages("SRAdb", repos="http://cran.us.r-project.org")'
# /build/dirR-$version-$arch/R --vanilla -e 'install.packages("factoextra", repos="http://cran.us.r-project.org")'
/build/dirR-$version-$arch/R --vanilla -e 'install.packages("MatrixEQTL", repos="http://cran.us.r-project.org")'


echo "Finished building dirR4"

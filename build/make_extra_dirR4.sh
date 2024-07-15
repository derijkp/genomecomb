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

#dirRversion=3.5.3
#majorversion=3
dirRversion=4.2.1
majorversion=4

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
yuminstall atlas
#yuminstall java
#yuminstall texlive
yuminstall texlive-latex
yuminstall fontconfig
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

cd /build
rm -rf /build/R-$dirRversion
download https://cran.r-project.org/src/base/R-$majorversion/R-$dirRversion.tar.gz
cd /build/R-$dirRversion
make distclean || true

# make sure the manually installed libs (in /usr/local) are used
cd /build/R-$dirRversion
PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    CPPFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    ./configure  --enable-static --disable-java --enable-R-shlib --prefix=/build/dirR-$dirRversion-$arch

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
cp -a -f /usr/lib64/libgfortran*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libgomp*.so.1* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/*pango* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /lib64/libreadline*.so.6* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libquadmath.so.0* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /build/zlib-1.2.11/libz.so.1* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/local/lib64/libssl*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/local/lib64/libcrypto*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/local/lib/libcurl*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
# cp -a -f /usr/local/lib/libpng*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libpng*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
# cp -a -f /usr/local/lib/libtiff*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libtiff*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
#cp -a -f /usr/local/lib/libpixman*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libpixman*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libcairo*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libfribidi*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
#cp -a -f /usr/local/lib/libxml2*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libxml2*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libjpeg*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tcl* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -fr /usr/share/tcl8.5 /build/dirR-$dirRversion-$arch/lib64/R/share
cp -a -fr /usr/share/tk8.5 /build/dirR-$dirRversion-$arch/lib64/R/share
#cp -a -f /usr/lib64/libpng12.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/liblzma* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicuuc.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicui18n.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libicudata.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libpcre2-8.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f /usr/lib64/libtinfo.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
# undefined symbol: g_log_structured_standard
cp -a -f /usr/lib64/libglib* /build/dirR-$dirRversion-$arch/lib64/R/lib
# undefined symbol: FcWeightFromOpenTypeDouble
cp -a -f /usr/lib64/libfontconfig* /build/dirR-$dirRversion-$arch/lib64/R/lib
# undefined symbol: hb_variation_from_string
cp -a -f /usr/lib64/libharfbuzz* /build/dirR-$dirRversion-$arch/lib64/R/lib
# undefined symbol: FT_Get_Var_Design_Coordinates
cp -a -f /usr/lib64/libfreetype* /build/dirR-$dirRversion-$arch/lib64/R/lib
# other deps
cd /usr/lib64
cp -a -f libgobject* libthai* libpcre.so libpcre.so.1 libpcre.so.1.2.0 \
    libffi* libgraphite2* libexpat* libuuid* libbz2* \
    libpixman* libEGL* libxcb-shm.* libxcb.so* libxcb-render.* \
    libXrender.* libX11* libXext* libGL* libXau* \
	/build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f libstdc++.so* libgcc_s* \
	/build/dirR-$dirRversion-$arch/lib64/R/lib
# docker test ubuntu:20.04
yuminstall libSM
cp -a -f libSM.* libgobject* \
	libXt.so.* libXft.so.* atlas/libsatlas.so.* libjbig.so.* libICE.so.* \
	/build/dirR-$dirRversion-$arch/lib64/R/lib

cd /usr/share
cp -a -f fonts fontconfig \
	/build/dirR-$dirRversion-$arch/lib64/R/share
mkdir -f /build/dirR-$dirRversion-$arch/lib64/R/conf || true
cp -a -f /etc/fonts/fonts.conf \
	/build/dirR-$dirRversion-$arch/lib64/R/conf

cd /build/R-$dirRversion

echo '
puts [list set argv $argv]
foreach {file basedir locate} $argv break
if {![file exists $file.ori]} {file copy $file $file.ori}
set f [open $file.ori]
set c [read $f]
close $f
regsub {# Shell wrapper for R executable.} $c {# Shell wrapper for R executable.
script="$(readlink -f "$0")"
R_BASEDIR="$(dirname "$(dirname "$script")")"
export R_LIBS="${R_BASEDIR}/lib64/R/library"
export R_LIBS_USER="${R_BASEDIR}/lib64/R/library"
export R_LIBS_SITE="${R_BASEDIR}/lib64/R/library"
} c
regsub -all $basedir $c {${R_BASEDIR}} c
regsub {export R_DOC_DIR
} $c {export R_DOC_DIR
export TCL_LIBRARY="${R_BASEDIR}/lib64/R/share/tcl8.5"
export TK_LIBRARY="${R_BASEDIR}/lib64/R/share/tk8.5"
} c
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

tclsh /tmp/convert /build/dirR-$dirRversion-$arch/bin/R /build/dirR-$dirRversion-$arch
cp /build/dirR-$dirRversion-$arch/bin/R.ori /build/dirR-$dirRversion-$arch/R.ori
tclsh /tmp/convert /build/dirR-$dirRversion-$arch/R /build/dirR-$dirRversion-$arch '$(dirname "$script")'
tclsh /tmp/convert /build/dirR-$dirRversion-$arch/lib64/R/bin/R /build/dirR-$dirRversion-$arch '$(dirname "$(dirname "$(dirname "$(dirname "$script")")")")'

# make an Rscript analog using only R and bash
cat << 'EOF' > Rscript
#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
# Initialize variables
R_OPTS=()
R_SCRIPT=""
SCRIPT_ARGS=()
# Parse the arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|--version|--vanilla|--no-save|--no-restore|-save|--no-environ|--no-site-file|--no-init-file|--restore)
            R_OPTS+=("$1")
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            if [ -z "$R_SCRIPT" ]; then
                R_SCRIPT="$1"
                shift
            else
                SCRIPT_ARGS+=("$1")
                shift
            fi
            ;;
    esac
done
# Check if the R script file is specified
if [ -z "$R_SCRIPT" ]; then
    echo "Usage: $0 [options] <R_script> [args...]"
    exit 1
fi
# Check if the R script file exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "R script '$R_SCRIPT' not found!"
    exit 1
fi
# run the script with arguments
"${dir}/R" "${R_OPTS[@]}" --slave -e "args <- commandArgs(trailingOnly = TRUE); source('$R_SCRIPT')" --args "${SCRIPT_ARGS[@]}"
EOF
chmod ugo+x Rscript

# packages
# --------

#if [ ! -f /build/dirR-$dirRversion-$arch/lib64/R/etc/Makeconf.ori ] ; then
#	cp /build/dirR-$dirRversion-$arch/lib64/R/etc/Makeconf /build/dirR-$dirRversion-$arch/lib64/R/etc/Makeconf.ori
#fi
#sed 's/CXX11 =/CXX11 = g++ -std=c++11/g' /build/dirR-$dirRversion-$arch/lib64/R/etc/Makeconf.ori \
#    | sed 's#CXX11FLAGS =  $(LTO)#CXX11FLAGS = -I/usr/local/include -fPIC#g' \
#    > /build/dirR-$dirRversion-$arch/lib64/R/etc/Makeconf

yuminstall libtiff-devel
cp -a -f /usr/lib64/libtiff.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
/build/dirR-$dirRversion-$arch/R --vanilla -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'

/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("askpass", version = "1.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("bit", version = "4.0.4", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("blob", version = "1.2.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("brio", version = "1.1.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("cowplot", version = "1.1.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("cpp11", version = "0.4.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("DBI", version = "1.1.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("dbplyr", version = "2.3.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("dtplyr", version = "1.2.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("dqrng", version = "0.3.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("bslib", version = "0.5.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("callr", version = "3.7.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("checkmate", version = "2.2.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("cli", version = "3.6.1", repos = "http://cran.us.r-project.org")'

/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("tidyverse", version = "2.0.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("rmarkdown", version = "2.23", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("shiny", version = "1.7.4", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("googleVis", version = "0.7.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("DT", version = "0.28", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("RColorBrewer", version = "1.1-3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("pheatmap", version = "1.0.12", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("curl", version = "5.0.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("plotly", version = "4.10.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("jpeg", version = "0.1-10", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("xml2", version = "1.3.4", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("XML", version = "3.99-0.14", repos = "http://cran.us.r-project.org")'
yuminstall cmake3
yuminstall NLopt NLopt-devel
cp -a -f /usr/lib64/libnlopt*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib


/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("nloptr", version = "2.0.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("lme4", version = "1.1-33", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("pbkrtest", version = "0.5.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("car", version = "3.1-2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("rstatix", version = "0.7.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("ggpubr", version = "0.6.0", repos = "http://cran.us.r-project.org")'

/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("BiocManager", version = "1.30.21", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("DESeq2")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("edgeR")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'install.packages("https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz",repos=NULL,type="source");'
# /build/dirR-$dirRversion-$arch/R --vanilla -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",repos=NULL,type="source");'
# get java part with: wget http://hartleys.github.io/QoRTs/QoRTs.jar

/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("PCAtools")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("pcaExplorer")'
/build/dirR-$dirRversion-$arch/R --vanilla -e $'BiocManager::install("Rhtslib")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("Rsamtools")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("ggbio")'

/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("bambu")'

# leafcutter
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("rstan", version = "2.21.8", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("rstantools", version = "2.3.1", repos = "http://cran.us.r-project.org")'
yuminstall gsl-devel
cp -a /usr/lib64/libgsl*.so* /build/dirR-$dirRversion-$arch/lib64/R/lib
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("DirichletMultinomial")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("DropletUtils")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter", version = "0.2.9")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("Seurat", version = "4.3.0.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("sleepwalk", version = "0.3.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("SCINA", version = "1.2.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("HGNChelper", version = "0.8.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("openxlsx", version = "4.2.5.2", repos = "http://cran.us.r-project.org")'

# other
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("MOFA2")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("mixOmics")'
# Rhtslib does not compile for some reason, so not for now
# /build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("chromVAR")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("OmicCircos")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("RCircos", version = "1.2.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("ggdendro", version = "0.1.23", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("GGally", version = "2.1.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("ComplexUpset", version = "1.3.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("corrplot", version = "0.92", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("eulerr", version = "7.0.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("vioplot", version = "0.4.0", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("gplots", version = "3.1.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("gtools", version = "3.9.4", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("ggrepel", version = "0.9.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("MatrixEQTL", version = "2.3", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'install.packages("https://github.com/dzhang32/ggtranscript/archive/refs/tags/v0.99.3.tar.gz", repos=NULL, type="source")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("survminer", version = "0.4.9", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("scDblFinder")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", version = "2.0.3")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("argparser", version = "0.7.1", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("scSorter", version = "0.0.2", repos = "http://cran.us.r-project.org")'
/build/dirR-$dirRversion-$arch/R --vanilla -e 'BiocManager::install("cleaver")'

cd /build/R-$dirRversion
yuminstall hdf5-devel
wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_14_1/src/hdf5-1.14.1-2.tar.gz
tar xvzf hdf5-1.14.1-2.tar.gz
cd /build/R-$dirRversion/hdf5-1.14.1-2
./configure
make
make install
sudo cp -ra hdf5/bin/* /usr/bin
sudo cp -ra hdf5/lib/* /usr/lib
sudo cp -ra hdf5/share/* /usr/share/
sudo cp -ra hdf5/include/* /usr/include/
cp -ra hdf5/lib/* /build/dirR-$dirRversion-$arch/lib64/R/lib
/build/dirR-$dirRversion-$arch/R --vanilla -e 'require(devtools); install_version("hdf5r", version = "1.3.8", repos = "http://cran.us.r-project.org")'
cp -a -f hdf5/lib/* /build/dirR-$dirRversion-$arch/lib64/R/lib
cp -a -f hdf5/bin/* /build/dirR-$dirRversion-$arch/lib64/R/bin

# "install" sc-type
cd /build/dirR-$dirRversion-$arch/lib64/R/library
wget https://github.com/IanevskiAleksandr/sc-type/archive/refs/tags/v1.0.tar.gz
tar xvzf v1.0.tar.gz
rm v1.0.tar.gz

cd /build
ln -sf dirR-$dirRversion-$arch/R dirR-4.2.1
ln -sf dirR-$dirRversion-$arch/R dirR
ln -sf dirR-$dirRversion-$arch/R R
rm dirR-$dirRversion-$arch.tar.gz || true
tar cvzf dirR-$dirRversion-$arch.tar.gz dirR-$dirRversion-$arch dirR-4.2.1 dirR R
cp -ra dirR-$dirRversion-$arch dirR-4.2.1 dirR R /io/extra$ARCH

echo "Finished building dirR4 in $builddir/dirR-$dirRversion-$arch"
echo "Copied to $srcdir/extra$ARCH"

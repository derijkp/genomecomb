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
# source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# Build
# =====

# set up environment
# ------------------
yuminstall git
yuminstall wget
yuminstall centos-release-scl
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-8
yuminstall rh-python36
# scl enable devtoolset-8 bash
scl enable devtoolset-8 rh-python36 bash


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

# vbz_compression
# ===============

yuminstall gcc-c++
#yuminstall zstd
#yuminstall libzstd-devel
yuminstall libffi
yuminstall libffi-devel
#yuminstall hdf5
#yuminstall hdf5-devel
sudo python3 -m pip install wheel
PATH=/hbb/bin:$PATH
# PATH=/build/cmake-3.16.4-Linux-x86_64/bin:$PATH

# zstd (for lib)
# --------------
cd /build
if [ ! -f "/build/zstd-1.4.4/lib/libzstd.a" ] ; then
	source /hbb_shlib/activate
	wget -c https://github.com/facebook/zstd/releases/download/v1.4.4/zstd-1.4.4.tar.gz
	tar xvzf zstd-1.4.4.tar.gz
	cd /build/zstd-1.4.4
	CFLAGS="-O3 -fPIC" make
	cd /build/zstd-1.4.4/contrib/seekable_format/examples
	CFLAGS="-O3 -fPIC" make
	source /hbb_exe/activate
fi

## cmake
## ----
#wget https://github.com/Kitware/CMake/releases/download/v3.16.4/cmake-3.16.4-Linux-x86_64.sh
#bash cmake-3.16.4-Linux-x86_64.sh --skip-license --include-subdir

# hdf5
# ----
# hack to get recent enough hdf5 version available to cmake (cmake ignores HDF5_ROOT)
cd /build
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.16/bin/linux-centos6-x86_64-gcc447/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static.tar.gz
tar xvzf hdf5-1.8.16-linux-centos6-x86_64-gcc447-static.tar.gz
# sudo cp -ra /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/* /usr

##cd /build
### This wget does not actually work, need to get a real url
##wget https://www.hdfgroup.org/package/hdf5-1-10-6-tar-gz
##ln -s hdf5-1-10-6-tar-gz hdf5-1-10-6.tar.gz
##tar xvzf hdf5-1.10.6.tar.gz
##cd /build/hdf5-1.10.6
### ./configure --enable-static --disable-shared --disable-hd --prefix=$HOME
##./configure --enable-static --disable-shared
##make
##make install

# google/benchmark
# ----------------
cd /build
wget https://github.com/google/benchmark/archive/v1.5.0.tar.gz
tar xvzf v1.5.0.tar.gz
rm v1.5.0.tar.gz
cd /build/benchmark-1.5.0
mkdir build
cd /build/benchmark-1.5.0/build
cmake -D BENCHMARK_ENABLE_TESTING=OFF ../
make
sudo make install

# vbz_compression
# ---------------

cd /build
wget https://github.com/nanoporetech/vbz_compression/archive/v1.0.0.zip
mv v1.0.0.zip vbz_compression-1.0.0.zip
unzip vbz_compression-1.0.0.zip
cd /build/vbz_compression-1.0.0/third_party/streamvbyte
wget https://github.com/lemire/streamvbyte/archive/v0.4.1.tar.gz
tar xvzf v0.4.1.tar.gz
mv streamvbyte-0.4.1/* .

# I don't know how to get cmake to add needed libraries, so do it this
# terribly hacky way; but since we're in a container it should work
mv /build/zstd-1.4.4/lib/libzstd.so /build/zstd-1.4.4/lib/libzstd.so.save || true

mkdir /build/vbz_compression-1.0.0/build
cd /build/vbz_compression-1.0.0/build
rm -rf /build/vbz_compression-1.0.0/build/*

hdf5dir=/build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static
cmake \
    -D CMAKE_PREFIX_PATH='/build/zstd-1.4.4;/build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static' \
    -D ZSTD_INCLUDE_DIR=/build/zstd-1.4.4/lib \
    -D ZSTD_LIBRARY=/build/zstd-1.4.4/lib/libzstd.a \
    -D HDF5_C_LIBRARY=$hdf5dir \
    -D HDF5_C_LIBRARY=$hdf5dir/lib \
    -D HDF5_C_INCLUDE_DIR=$hdf5dir/include \
    -D ENABLE_CONAN=OFF \
    ..

make VERBOSE=1
cd /build/vbz_compression-1.0.0/build/vbz/perf
/opt/rh/devtoolset-8/root/usr/bin/c++  -g   CMakeFiles/vbz_perf_test.dir/vbz_perf.cpp.o  -o ../../bin/vbz_perf_test ../../lib/libvbz.a /usr/local/lib64/libbenchmark.a ../../streamvbyte/src/streamvbyte-build/./libstreamvbyte_static.a /build/zstd-1.4.4/lib/libzstd.a -lpthread -lrt
cd /build/vbz_compression-1.0.0/build
make VERBOSE=1
cd /build/vbz_compression-1.0.0/build/vbz_plugin/test
/opt/rh/devtoolset-8/root/usr/bin/c++  -g   CMakeFiles/vbz_hdf_plugin_test.dir/vbz_hdf_plugin_test.cpp.o CMakeFiles/vbz_hdf_plugin_test.dir/main.cpp.o  -o ../../bin/vbz_hdf_plugin_test -Wl,-rpath,/build/vbz_compression-1.0.0/build/bin ../../bin/libvbz_hdf_plugin.so /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libhdf5.a ../../lib/libhdf_test_utils.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libhdf5.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libz.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libsz.a -ldl
cd /build/vbz_compression-1.0.0/build
make VERBOSE=1
cd /build/vbz_compression-1.0.0/build/vbz_plugin/perf
/opt/rh/devtoolset-8/root/usr/bin/c++  -g   CMakeFiles/vbz_hdf_perf_test.dir/vbz_hdf_perf.cpp.o  -o ../../bin/vbz_hdf_perf_test -Wl,-rpath,/build/vbz_compression-1.0.0/build/bin /usr/local/lib64/libbenchmark.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libhdf5.a ../../lib/libhdf_test_utils.a ../../bin/libvbz_hdf_plugin.so -lpthread /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libhdf5.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libz.a /build/hdf5-1.8.16-linux-centos6-x86_64-gcc447-static/lib/libsz.a -ldl -lrt
cd /build/vbz_compression-1.0.0/build
make VERBOSE=1

echo "Finished building vbz_compression"

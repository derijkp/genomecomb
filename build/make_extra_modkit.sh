#!/bin/bash

# This script downloads a build of modkit, and puts it in a portable directory
# 32 bit option is ignored
# options:
# -b|-bits|--bits: 32 for 32 bits build (default 64)
# -d|-builddir|--builddir: top directory to build external software in (default ~/build/bin-$arch)

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

# settings
# ========

modkitversion=0.4.2

# Script run within Holy Build box
# ================================

echo "Entering Holy Build Box environment"

# Activate Holy Build Box environment.
# does not compile with these settings (X)
# only use HBB for glibc compat, not static libs
# source /hbb_exe/activate

# print all executed commands to the terminal
set -x

# Build
# =====

# set up environment
# ------------------
yuminstall wget
# yuminstall openssl-devel
yuminstall perl-devel
yuminstall perl-IPC-Cmd

# yuminstall xz
sudo yum install curl gcc -y
yuminstall centos-release-scl
sudo yum upgrade -y
# sudo yum list all | grep devtoolset
yuminstall devtoolset-11
# use source instead of scl enable so it can run in a script
source /opt/rh/devtoolset-11/enable


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

if [ ! "$arch" = "linux-x86_64"] ; then
	echo "Only linux-x86_64 supported currently"
	exit 1
fi


## openssl
## -------
#cd /build
#download https://www.openssl.org/source/openssl-3.4.0.tar.gz
#cd /build/openssl-3.4.0
#make clean || true
#make distclean || true
#./config -enable-static
##./config
#make CFLAGS="-fPIC"
#sudo make install
##sudo rm /usr/local/lib64/libssl.so* /usr/local/lib64/libcrypto.so*

# rust/cargo
# ----------
cd /build
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

export PATH="$HOME/.cargo/bin:$PATH"
source $HOME/.cargo/env
rustup install 1.83.0

# modkit
# ------
cd /build

wget https://github.com/nanoporetech/modkit/archive/refs/tags/v${modkitversion}.tar.gz
tar xvzf v${modkitversion}.tar.gz
rm v${modkitversion}.tar.gz

rm -rf /build/modkit-${modkitversion}-build || true
mv /build/modkit-${modkitversion} /build/modkit-${modkitversion}-build
cd /build/modkit-${modkitversion}-build
rm -rf /build/modkit-${modkitversion}-$arch || true
mkdir /build/modkit-${modkitversion}-$arch
cargo install --root /build/modkit-${modkitversion}-$arch --path .
cp -ra README.md LICENCE.txt docs /build/modkit-${modkitversion}-$arch
mv /build/modkit-${modkitversion}-$arch/bin/modkit /build/modkit-${modkitversion}-$arch
rmdir /build/modkit-${modkitversion}-$arch/bin
cd /build
ln -sf modkit-${modkitversion}-$arch/modkit modkit
ln -sf modkit-${modkitversion}-$arch/modkit modkit-$modkitversion
tar cvzf modkit-${modkitversion}-$arch.tar.gz modkit-${modkitversion}-$arch modkit modkit-${modkitversion}

echo "Finished building modkit-${modkitversion}-$arch in $builddir"


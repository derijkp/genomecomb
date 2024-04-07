#!/bin/bash

# This script makes a portable java
# using the Holy build box environment to create the portable binaries.
# The build directory contains a (patched) tcc that will find out from which dir it is run, and use the included 
# musl library for compiling portable, statically linked executables. 
# The tcc-system will use system libraries (in /usr/lib/$arch-linux-gnu) instead
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
source "${dir}/start_hbb.sh"

# settings
# ========

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
yuminstall java

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

# java 1.8
cd /build
cp -raL /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.275.b01-0.el6_10.x86_64 /build/java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64
cp -ra /usr/lib64/libfreetype.so* /build/java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64/jre/lib/amd64
ln -sf java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64/jre/bin/java java1.8
tar cvzf java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64.tar.gz java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64 java1.8

wget https://github.com/ojdkbuild/ojdkbuild/releases/download/1.8.0.121-1/java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64.zip
unzip java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64.zip
rm java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64.zip
cp -ra /usr/lib64/libfreetype.so* /build/java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64/jre/lib/amd64
cp -ra /usr/lib64/libfontconfig.so* /build/java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64/jre/lib/amd64
echo 'version=1
filename.Arial=fonts/Arial.ttf
filename.Arial_Bold=/fonts/ArialBold.ttf
filename.Arial_Italic=/fonts/ArialItalic.ttf
filename.Arial_Bold_Italic=/fonts/ArialBoldItalic.ttf
' > /build/java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64/jre/lib/fontconfig.properties
ln -sf java-1.8.0-openjdk-1.8.0.121-0.b13.el6_8.x86_64/jre/bin/java java1.8

# java 22
wget https://download.java.net/java/GA/jdk22/830ec9fcccef480bb3e73fb7ecafe059/36/GPL/openjdk-22_linux-x64_bin.tar.gz
tar xvzf openjdk-22_linux-x64_bin.tar.gz
rm openjdk-22_linux-x64_bin.tar.gz
mv jdk-22 openjdk-22-linux-x86_64
cp -ra /usr/lib64/libfreetype.so* /build/openjdk-22-linux-x86_64/lib
echo 'version=1
filename.Arial=/home/frankie/fonts/Arial.ttf
filename.Arial_Bold=/home/frankie/fonts/ArialBold.ttf
filename.Arial_Italic=/home/frankie/fonts/ArialItalic.ttf
filename.Arial_Bold_Italic=/home/frankie/fonts/ArialBoldItalic.ttf
' > /build/openjdk-22-linux-x86_64/lib/fontconfig.properties
ln -sf openjdk-22-linux-x86_64/bin/java .
ln -sf openjdk-22-linux-x86_64/bin/java java22
tar cvzf openjdk-22-linux-x86_64.tar.gz openjdk-22-linux-x86_64 java java22

echo "made $builddir/java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64.tar.gz"
echo "made $builddir/openjdk-22-linux-x86_64.tar.gz"
echo "Finished building java"


#!/bin/bash

# This script builds a portable tiny (6.4MiB on my system) build environment based on tcc and musl,
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
yuminstall git
yuminstall wget

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


# tcc
# ---
cd /build

# release version does not work with musl, used master version (at given patch) instead (next)
tccpatch=cc40305a121246116d5e08037e94942a9c0ec1fb
tccversion=master2020-12-25

cd /build
rm -rf tinycc.old || true
mv tinycc tinycc.old || true
git clone https://github.com/TinyCC/tinycc.git
cd /build/tinycc
git reset --hard $tccpatch

# patch libtcc.c to use env variables TCC_CRT_PATH and TCC_LIB_PATH
if [ ! -f "libtcc.c.ori" ] ; then
    echo "copy libtcc.c to libtcc.c.ori"
    cp libtcc.c libtcc.c.ori
fi
echo "patching libtcc.c"
cp libtcc.c.ori libtcc.c
pattern='tcc_split_path(s, &s->crt_paths, &s->nb_crt_paths, CONFIG_TCC_CRTPREFIX);'
replace=`echo '{
    char *path;
    path = getenv("TCC_CRT_PATH");
    if(path != NULL) {
        tcc_split_path(s, \\&s->crt_paths, \\&s->nb_crt_paths, path);
    } else {
        tcc_split_path(s, \\&s->crt_paths, \\&s->nb_crt_paths, CONFIG_TCC_CRTPREFIX);
    }
    }' | awk 1 ORS='\\\\n'`
sed -i "s/$pattern/$replace/" libtcc.c

pattern='tcc_set_lib_path(s, CONFIG_TCCDIR);'
replace=`echo '{
    char *path;
    path = getenv("TCC_LIB_PATH");
    if(path != NULL) {
        tcc_set_lib_path(s, path);
    } else {
        tcc_set_lib_path(s, CONFIG_TCCDIR);
    }
    }' | awk 1 ORS='\\\\n'`
sed -i "s/$pattern/$replace/" libtcc.c

./configure --enable-cross --strip-binaries --prefix=/build/tcc-$tccversion-$arch --disable-shared --enable-static
make
make test
make install

# musl
# ----

download http://musl.libc.org/releases/musl-1.2.1.tar.gz
cd /build/musl-1.2.1
./configure --prefix=/build/musl-1.2.1-$arch --disable-shared --enable-static
make
make install

# make final portable tcc toolchain
# ---------------------------------
rm -rf /io/extra$ARCH/tcc-$tccversion-$arch
mkdir -p /io/extra$ARCH/tcc-$tccversion-$arch/bin
cp -ra /build/tcc-$tccversion-$arch/bin/tcc /io/extra$ARCH/tcc-$tccversion-$arch/bin
cp -ra /build/tcc-$tccversion-$arch/include /io/extra$ARCH/tcc-$tccversion-$arch
cp -ra /build/tcc-$tccversion-$arch/lib /io/extra$ARCH/tcc-$tccversion-$arch
cp -ra /build/tcc-$tccversion-$arch/share /io/extra$ARCH/tcc-$tccversion-$arch

cp -ra /build/musl-1.2.1-$arch /io/extra$ARCH/tcc-$tccversion-$arch/lib

# make "executable" shell script that sets the env variables
echo '#!/bin/sh
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
TCC_LIB_PATH=$dir/lib/tcc \
TCC_CRT_PATH=$dir/lib/musl-1.2.1-$arch/lib \
C_INCLUDE_PATH=$dir/lib/musl-1.2.1-$arch/include:$C_INCLUDE_PATH \
LIBRARY_PATH=$dir/lib/tcc:$dir/lib/musl-1.2.1-$arch/lib:$LIBRARY_PATH \
    $dir/bin/tcc ${1+"$@"}
' > /io/extra$ARCH/tcc-$tccversion-$arch/tcc
chmod ugo+x /io/extra$ARCH/tcc-$tccversion-$arch/tcc

# make alternative "executable" shell script that uses the system libraries iso musl
echo '#!/bin/sh
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
SYSTEM_LIB=/usr/lib/$arch-linux-gnu
TCC_LIB_PATH=$dir/lib/tcc \
TCC_CRT_PATH=$SYSTEM_LIB \
C_INCLUDE_PATH=$dir/lib/tcc/include:$C_INCLUDE_PATH \
LIBRARY_PATH=$dir/lib/tcc:$SYSTEM_LIB:$LIBRARY_PATH \
    $dir/bin/tcc ${1+"$@"}
' > /io/extra$ARCH/tcc-$tccversion-$arch/tcc-system
chmod ugo+x /io/extra$ARCH/tcc-$tccversion-$arch/tcc

# link to tcc
cd /io/extra$ARCH/
ln -sf tcc-$tccversion-$arch/tcc tcc
ln -sf tcc-$tccversion-$arch/tcc-system tcc-system

echo "Finished building tcc-$tccversion-$arch"


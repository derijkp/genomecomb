#!/bin/bash

# This script downloads portable julia binaries and adds commands
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

# Parse arguments
# ===============

juliaversion=1.12.6
arch=linux-x86_64

cd $dir/../extra
wget https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-$juliaversion-$arch.tar.gz
tar xvzf julia-$juliaversion-$arch.tar.gz
rm julia-$juliaversion-$arch.tar.gz
mv julia-$juliaversion julia-$juliaversion-$arch
cd julia-$juliaversion-$arch

echo '#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
export LANG=C
export LC_ALL=C
$dir/bin/julia $dir/share/julia/juliac/juliac.jl ${1+"$@"}
' > juliac
chmod ugo+x juliac

cd ..
ln -sf julia-$juliaversion-$arch/bin/julia .
ln -sf julia-$juliaversion-$arch/bin/julia julia-$juliaversion
ln -sf julia-$juliaversion-$arch/juliac .
ln -sf julia-$juliaversion-$arch/juliac juliac-$juliaversion

echo "Finished building julia-$juliaversion-$arch"

exit 0

# how to use/test juliac
echo '@main function main(args)
    println("Hello, World!")
    return 0
end
' > helloworld.jl

./juliac --output-exe hello helloworld.jl 

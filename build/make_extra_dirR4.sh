#!/bin/bash

# This script builds a distributable directory contained version of R using the Holy build box environment
# options:
# -a|-arch|--arch: 64, x86_64 or linux-x86_64 for 64 bit Linux build (default); ix86, 32 or linux-ix86 for 32 bits Linux build; win, windows-x86_64 or mingw-w64 for Windows 64 bit build
# -b|-bits|--bits: select 32 or 64 bits Linux build (default 64 = same as -arch x86_64)
# -d|-builddir|--builddir: top directory to build in (default ~/build/bin-$arch)

# The Holy build box environment requires docker, make sure it is installed
# e.g. on ubuntu and derivatives
# sudo apt install docker.io
# Also make sure you have permission to use docker
# sudo usermod -a -G docker $USER

# stop on error
# best not to set this when debugging a build script by running interactively (copy/paste parts)
set -e

# Prepare and start docker with Holy Build box
# ============================================

script="$(readlink -f "$0")"
dir="$(dirname "$script")"
source "${dir}/start_hbb3.sh"

# Parse arguments
# ===============

Rversion=4.2.1
dirRversion=4.2.1-3
majorversion=4
dest_appdir=dirR-$dirRversion-$arch

all=1
while [[ "$#" -gt 0 ]]; do case $1 in
	*) echo "Unknown parameter: $1"; exit 1;;
esac; shift; done

if [ "$clear" = "1" ] ; then
    rm -rf /build/$dest_appdir || true
fi

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

# bzip2
# -----
cd /build
download https://sourceforge.net/projects/bzip2/files/latest/bzip2-1.0.6.tar.gz
cd /build/bzip2-1.0.6
make clean || true
make CFLAGS="-fPIC"
sudo make install

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

# libxml2
# -------
cd /build
download ftp://xmlsoft.org/libxml2/libxml2-2.9.9.tar.gz
cd /build/libxml2-2.9.9

make distclean || true
PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH ./configure --enable-static --without-python
make CFLAGS="-fPIC"
sudo make install

# build R
# -------
# useful hints in getting it compiled found at https://tdhock.github.io/blog/2017/compiling-R/

cd /build
rm -rf /build/R-$Rversion
download https://cran.r-project.org/src/base/R-$majorversion/R-$Rversion.tar.gz
cd /build/R-$Rversion
make distclean || true

# make sure the manually installed libs (in /usr/local) are used
cd /build/R-$Rversion
PATH=/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
    CFLAGS="-I/usr/local/include -fPIC" \
    CPPFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-L/usr/local/lib -L/usr/local/lib64" \
    ./configure  --enable-static --disable-java --enable-R-shlib --prefix=/build/$dest_appdir

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
cp -a -f /usr/lib64/libgfortran*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libgomp*.so.1* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/*pango* /build/$dest_appdir/lib64/R/lib
cp -a -f /lib64/libreadline*.so.6* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libquadmath.so.0* /build/$dest_appdir/lib64/R/lib
cp -a -f /build/zlib-1.2.11/libz.so.1* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/local/lib64/libssl*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/local/lib64/libcrypto*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/local/lib/libcurl*.so* /build/$dest_appdir/lib64/R/lib
# cp -a -f /usr/local/lib/libpng*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libpng*.so* /build/$dest_appdir/lib64/R/lib
# cp -a -f /usr/local/lib/libtiff*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libtiff*.so* /build/$dest_appdir/lib64/R/lib
#cp -a -f /usr/local/lib/libpixman*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libpixman*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libcairo*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libfribidi*.so* /build/$dest_appdir/lib64/R/lib
#cp -a -f /usr/local/lib/libxml2*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libxml2*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libjpeg*.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/*tcl* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/*tk* /build/$dest_appdir/lib64/R/lib
cp -a -fr /usr/share/tcl8.5 /build/$dest_appdir/lib64/R/share
cp -a -fr /usr/share/tk8.5 /build/$dest_appdir/lib64/R/share
#cp -a -f /usr/lib64/libpng12.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/liblzma* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libicuuc.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libicui18n.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libicudata.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libpcre2-8.so* /build/$dest_appdir/lib64/R/lib
cp -a -f /usr/lib64/libtinfo.so* /build/$dest_appdir/lib64/R/lib
# undefined symbol: g_log_structured_standard
cp -a -f /usr/lib64/libglib* /build/$dest_appdir/lib64/R/lib
# undefined symbol: FcWeightFromOpenTypeDouble
cp -a -f /usr/lib64/libfontconfig* /build/$dest_appdir/lib64/R/lib
# undefined symbol: hb_variation_from_string
cp -a -f /usr/lib64/libharfbuzz* /build/$dest_appdir/lib64/R/lib
# undefined symbol: FT_Get_Var_Design_Coordinates
cp -a -f /usr/lib64/libfreetype* /build/$dest_appdir/lib64/R/lib
# other deps
cd /usr/lib64
cp -a -f libgobject* libthai* libpcre.so libpcre.so.1 libpcre.so.1.2.0 \
    libffi* libgraphite2* libexpat* libuuid* libbz2* \
    libpixman* libEGL* libxcb-shm.* libxcb.so* libxcb-render.* \
    libXrender.* libX11* libXext* libGL* libXau* \
	/build/$dest_appdir/lib64/R/lib
cp -a -f libstdc++.so* libgcc_s* \
	/build/$dest_appdir/lib64/R/lib
# docker test ubuntu:20.04
yuminstall libSM
cp -a -f libSM.* libgobject* \
	libXt.so.* libXft.so.* atlas/libsatlas.so.* libjbig.so.* libICE.so.* \
	/build/$dest_appdir/lib64/R/lib

cd /usr/share
cp -a -f fonts fontconfig \
	/build/$dest_appdir/lib64/R/share
mkdir -f /build/$dest_appdir/lib64/R/conf || true
cp -a -f /etc/fonts/fonts.conf \
	/build/$dest_appdir/lib64/R/conf

cd /build/R-$Rversion

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

tclsh /tmp/convert /build/$dest_appdir/bin/R /build/$dest_appdir
cp /build/$dest_appdir/bin/R.ori /build/$dest_appdir/R.ori
tclsh /tmp/convert /build/$dest_appdir/R /build/$dest_appdir '$(dirname "$script")'
tclsh /tmp/convert /build/$dest_appdir/lib64/R/bin/R /build/$dest_appdir '$(dirname "$(dirname "$(dirname "$(dirname "$script")")")")'

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

# List all installed packages (and their versions) in order of dependency
# -----------------------------------------------------------------------
# The following script helps in determining the versions and order in which to install packages
# based on installed packages of an existing R dir
# it will output the commands needed to install all packages at the specific versions
# base packages and packages not installed through the default R package install (e.g. Bioconductor packages) are not included

# make installed_packages.R
cat << 'EOF' > installed_packages.R
# Topologically sort packages by their dependencies
topological_sort <- function(pkgs, deps) {
  # Create a set of sorted packages
  sorted <- character(0)
  # Recursive function to process packages
  visit <- function(pkg) {
    if (!pkg %in% sorted) {
      # Process dependencies first
      for (dep in deps[[pkg]]) {
        if (dep %in% names(deps)) visit(dep)
      }
      # Add package to sorted list
      sorted <<- c(sorted, pkg)
    }
  }
  # Process all packages
  for (pkg in pkgs) {
    visit(pkg)
  }
  
  return(sorted)
}
# list packages
list.installed <- function() {
    options(repos = c(CRAN = "http://cloud.r-project.org"))
    library(tools)
    # Get all installed packages
    installed <- installed.packages()
    # Filter out base packages
    non_base <- installed[is.na(installed[, "Priority"]) | installed[, "Priority"] == "recommended",]
    recommended <- rownames(installed[!is.na(installed[, "Priority"]) & installed[, "Priority"] == "recommended",])
    # Get the package names and versions
    package_names <- rownames(non_base)
    package_versions <- non_base[, "Version"]
    # Get dependency information
    dependencies <- package_dependencies(package_names, recursive = TRUE)
    # Perform the topological sort
    sorted_packages <- topological_sort(package_names, dependencies)
    # Create a data frame of sorted packages and their versions
    sorted_installed <- data.frame(
      Package = sorted_packages,
      Version = package_versions[match(sorted_packages, package_names)],
      stringsAsFactors = FALSE
    )
    # filter out bioconductor packages
    available <- rownames(available.packages(filter=NULL,CRAN = "http://cloud.r-project.org"))
    # recommended packages can fall out (Matrix), add them
    available = c(available, recommended)
    sorted_installed = sorted_installed[sorted_installed[,1] %in% available,]
    # View the result
    sorted_installed = unname(sorted_installed)
    # output
    cat('export R_REMOTES_UPGRADE="never"',"\n")
    # print(sorted_installed)
    for (i in 1:nrow(sorted_installed)) {
        cat(paste("/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version(\"",sorted_installed[i,1],"\", version = \"",sorted_installed[i,2],"\", repos=\"http://cran.us.r-project.org\")'\n",sep=""))
    }
}
# run function
list.installed()
EOF
# run with
# ./Rscript installed_packages.R
# or
# ./R --vanilla < installed_packages.R
# where R or Rscript is the version of R you want to check installed packages

# packages
# --------

# install some libraries needed
yuminstall libtiff-devel
cp -a -f /usr/lib64/libtiff.so* /build/$dest_appdir/lib64/R/lib
# for nloptr
yuminstall cmake3
yuminstall NLopt NLopt-devel
cp -a -f /usr/lib64/libnlopt*.so* /build/$dest_appdir/lib64/R/lib

# install remotes to install specific versions
# needed for remotes::install_version
/build/$dest_appdir/R --vanilla -e 'install.packages("remotes", repos="http://cran.us.r-project.org")'

export R_REMOTES_UPGRADE="never" 
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("abind", version = "1.4-5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("argparser", version = "0.7.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sys", version = "3.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("askpass", version = "1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("assertthat", version = "0.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("backports", version = "1.4.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("base64enc", version = "0.1-3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("beeswarm", version = "0.4.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("BH", version = "1.81.0-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("BiocManager", version = "1.30.21", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("bit", version = "4.0.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("bit64", version = "4.0.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("bitops", version = "1.0-7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rlang", version = "1.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cli", version = "3.6.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("glue", version = "1.6.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lifecycle", version = "1.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("vctrs", version = "0.6.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("blob", version = "1.2.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("boot", version = "1.3-28", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("brew", version = "1.0-7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("brio", version = "1.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("generics", version = "0.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("magrittr", version = "2.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fansi", version = "1.0.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("utf8", version = "1.2.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pillar", version = "1.9.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("R6", version = "2.5.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pkgconfig", version = "2.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tibble", version = "3.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("withr", version = "2.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tidyselect", version = "1.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dplyr", version = "1.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("purrr", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("stringi", version = "1.7.12", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("stringr", version = "1.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cpp11", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tidyr", version = "1.3.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("broom", version = "1.0.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fastmap", version = "1.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cachem", version = "1.0.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("digest", version = "0.6.37", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("htmltools", version = "0.5.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("jquerylib", version = "0.1.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("jsonlite", version = "1.8.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("memoise", version = "2.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("mime", version = "0.12", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fs", version = "1.6.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rappdirs", version = "0.3.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sass", version = "0.4.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("bslib", version = "0.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ca", version = "0.71.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Cairo", version = "1.6-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ps", version = "1.7.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("processx", version = "3.8.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("callr", version = "3.7.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("carData", version = "3.0-5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Formula", version = "1.2-5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("MASS", version = "7.3-57", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lattice", version = "0.20-45", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("nlme", version = "3.1-157", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Matrix", version = "1.5-4.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("mgcv", version = "1.8-40", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("nnet", version = "7.3-17", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Rcpp", version = "1.0.11", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("minqa", version = "1.2.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("nloptr", version = "2.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppEigen", version = "0.3.3.9.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lme4", version = "1.1-33", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("numDeriv", version = "2016.8-1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gtable", version = "0.3.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("isoband", version = "0.2.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("farver", version = "2.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("labeling", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("colorspace", version = "2.1-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("munsell", version = "0.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RColorBrewer", version = "1.1-3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("viridisLite", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("scales", version = "1.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggplot2", version = "3.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cowplot", version = "1.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Deriv", version = "4.1.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("modelr", version = "0.1.11", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("microbenchmark", version = "1.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("doBy", version = "4.6.24", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pbkrtest", version = "0.5.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("SparseM", version = "1.81", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("MatrixModels", version = "0.5-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("survival", version = "3.3-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("quantreg", version = "5.94", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("car", version = "3.1-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("caTools", version = "1.18.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rematch", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cellranger", version = "1.1.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("checkmate", version = "2.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("class", version = "7.3-20", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("clipr", version = "0.8.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("cluster", version = "2.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("codetools", version = "0.2-18", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("commonmark", version = "1.9.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("patchwork", version = "1.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ComplexUpset", version = "1.3.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("conflicted", version = "1.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("corpcor", version = "1.6.10", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("corrplot", version = "0.92", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("crayon", version = "1.5.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("openssl", version = "2.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("curl", version = "6.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("credentials", version = "1.3.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lazyeval", version = "0.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("crosstalk", version = "1.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("data.table", version = "1.14.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("DBI", version = "1.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dbplyr", version = "2.3.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("deldir", version = "1.0-6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gridExtra", version = "2.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("viridis", version = "0.6.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dendextend", version = "1.16.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("desc", version = "1.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rstudioapi", version = "0.14", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("zip", version = "2.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gert", version = "1.9.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gitcreds", version = "0.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("httr2", version = "1.0.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ini", version = "0.3.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gh", version = "1.3.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rprojroot", version = "2.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("whisker", version = "0.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("yaml", version = "2.3.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("usethis", version = "2.1.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ellipsis", version = "0.3.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("later", version = "1.3.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("promises", version = "1.2.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("httpuv", version = "1.6.11", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("xtable", version = "1.8-4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fontawesome", version = "0.5.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sourcetools", version = "0.1.7-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shiny", version = "1.7.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("miniUI", version = "0.1.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pkgbuild", version = "1.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("evaluate", version = "0.21", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("downlit", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("systemfonts", version = "1.0.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("textshaping", version = "0.3.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ragg", version = "1.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("xfun", version = "0.39", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("highr", version = "0.10", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("knitr", version = "1.43", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tinytex", version = "0.45", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rmarkdown", version = "2.23", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("xml2", version = "1.3.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pkgdown", version = "2.0.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pkgload", version = "1.3.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("htmlwidgets", version = "1.6.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("profvis", version = "0.3.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("prettyunits", version = "1.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sessioninfo", version = "1.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("xopen", version = "1.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rcmdcheck", version = "1.4.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("remotes", version = "2.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("roxygen2", version = "7.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rversions", version = "2.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("praise", version = "1.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("diffobj", version = "0.3.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("waldo", version = "0.5.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("testthat", version = "3.1.9", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("urlchecker", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("devtools", version = "2.4.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dichromat", version = "2.0-0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("distributional", version = "0.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("iterators", version = "1.0.14", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("foreach", version = "1.5.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("doMC", version = "1.3.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("doParallel", version = "1.0.17", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dotCall64", version = "1.1-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sitmo", version = "2.0.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dqrng", version = "0.3.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("DT", version = "0.28", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("dtplyr", version = "1.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("egg", version = "0.4.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ellipse", version = "0.4.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("GenSA", version = "1.1.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("polyclip", version = "1.10-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("polylabelr", version = "0.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppArmadillo", version = "0.11.2.4.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("eulerr", version = "7.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("exactRankTests", version = "0.8-35", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spam", version = "2.10-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("maps", version = "3.4.1.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fields", version = "15.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("filelock", version = "1.0.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("fitdistrplus", version = "1.1-8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("FNN", version = "1.1.3.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("forcats", version = "1.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("foreign", version = "0.8-82", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("formatR", version = "1.12", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lambda.r", version = "1.2.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("futile.options", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("futile.logger", version = "1.4.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("globals", version = "0.16.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("listenv", version = "0.8.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("parallelly", version = "1.32.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("future", version = "1.28.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("future.apply", version = "1.9.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("httr", version = "1.4.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gargle", version = "1.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gclus", version = "1.3.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("getopt", version = "1.20.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("plyr", version = "1.8.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("hms", version = "1.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("progress", version = "1.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("GGally", version = "2.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("vipor", version = "0.4.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggbeeswarm", version = "0.7.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggdendro", version = "0.1.23", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggrepel", version = "0.9.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggsci", version = "2.9", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggsignif", version = "0.6.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("polynom", version = "1.4-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rstatix", version = "0.7.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggpubr", version = "0.6.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("png", version = "0.1-7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggrastr", version = "1.0.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggridges", version = "0.5.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("markdown", version = "1.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("jpeg", version = "0.1-10", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gridtext", version = "0.1.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ggtext", version = "0.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("goftest", version = "1.2-3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("uuid", version = "1.1-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("googledrive", version = "2.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ids", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rematch2", version = "2.1.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("googlesheets4", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("googleVis", version = "0.7.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gtools", version = "3.9.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("KernSmooth", version = "2.23-20", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gplots", version = "3.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("gridBase", version = "0.4-7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tzdb", version = "0.3.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("vroom", version = "1.6.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("readr", version = "2.1.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("haven", version = "2.5.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("hdf5r", version = "1.3.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("plotly", version = "4.10.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("reshape2", version = "1.4.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("qap", version = "0.1-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("registry", version = "0.5-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("TSP", version = "1.2-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("permute", version = "0.9-7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("vegan", version = "2.6-8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("seriation", version = "1.3.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("webshot", version = "0.5.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("heatmaply", version = "1.3.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("here", version = "1.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("HGNChelper", version = "0.8.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rpart", version = "4.1.16", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("htmlTable", version = "2.4.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Hmisc", version = "5.1-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ica", version = "1.0-3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("igraph", version = "1.3.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("inline", version = "0.3.19", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("interp", version = "1.1-3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("intervals", version = "0.15.4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("irlba", version = "2.3.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("R.methodsS3", version = "1.8.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("R.oo", version = "1.25.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("R.utils", version = "2.12.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("jrc", version = "0.5.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("km.ci", version = "0.5-6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("KMsurv", version = "0.1-5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("latticeExtra", version = "0.6-30", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppTOML", version = "0.1.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("reticulate", version = "1.26", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("leiden", version = "0.4.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("zoo", version = "1.8-11", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lmtest", version = "0.9-40", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("locfit", version = "1.5-9.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("matrixStats", version = "1.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tensorA", version = "0.36.2.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("posterior", version = "1.6.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("loo", version = "2.6.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("timechange", version = "0.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("lubridate", version = "1.9.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("MatrixEQTL", version = "2.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("mvtnorm", version = "1.2-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("maxstat", version = "0.7-25", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rngtools", version = "1.5.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("NMF", version = "0.28", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("oompaBase", version = "3.2.9", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("oompaData", version = "3.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("openxlsx", version = "4.2.5.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("optparse", version = "1.7.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pbapply", version = "1.5-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pheatmap", version = "1.0.12", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("pkgmaker", version = "0.32.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("plogr", version = "0.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("progressr", version = "0.11.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("QuickJSR", version = "1.4.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RANN", version = "2.6.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RSpectra", version = "0.16-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rARPACK", version = "0.11-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RCircos", version = "1.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppAnnoy", version = "0.0.19", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppHNSW", version = "0.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppML", version = "0.3.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppParallel", version = "5.1.7", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RcppProgress", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RCurl", version = "1.98-1.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("readxl", version = "1.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("reprex", version = "2.0.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("reshape", version = "0.8.9", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("XML", version = "3.99-0.14", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rjson", version = "0.2.21", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("restfulr", version = "0.0.15", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("ROCR", version = "1.0-11", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("RSQLite", version = "2.2.17", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("StanHeaders", version = "2.26.27", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rstan", version = "2.21.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rstantools", version = "2.3.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rsvd", version = "1.0.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Rtsne", version = "0.16", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("selectr", version = "0.4-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("rvest", version = "1.0.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("scattermore", version = "0.8", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("SCINA", version = "1.2.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("scSorter", version = "0.0.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sctransform", version = "0.3.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sp", version = "1.5-0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("SeuratObject", version = "4.1.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.utils", version = "3.1-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.data", version = "3.0-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.univar", version = "3.1-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.geom", version = "3.2-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.random", version = "3.1-5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tensor", version = "1.5", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.sparse", version = "3.0-2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatstat.explore", version = "3.2-1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("uwot", version = "0.1.14", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("Seurat", version = "4.3.0.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shinyAce", version = "0.4.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shinyBS", version = "0.61.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shinycssloaders", version = "1.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shinydashboard", version = "0.7.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("shinyjs", version = "2.1.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sleepwalk", version = "0.3.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("sm", version = "2.2-5.7.1", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("snow", version = "0.4-4", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("spatial", version = "7.3-15", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("statmod", version = "1.5.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("survMisc", version = "0.5.6", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("survminer", version = "0.4.9", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("TailRank", version = "3.2.2", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("threejs", version = "0.3.3", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("tidyverse", version = "2.0.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("vioplot", version = "0.4.0", repos="http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'library(remotes) ; remotes::install_version("xgboost", version = "1.6.0.1", repos="http://cran.us.r-project.org")'


# bioconductor packages
/build/$dest_appdir/R --vanilla -e 'require(devtools); install_version("BiocManager", version = "1.30.21", dependencies=FALSE, repos = "http://cran.us.r-project.org")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("edgeR")'
/build/$dest_appdir/R --vanilla -e 'install.packages("https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz",repos=NULL,type="source");'
# /build/$dest_appdir/R --vanilla -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",repos=NULL,type="source");'
# get java part with: wget http://hartleys.github.io/QoRTs/QoRTs.jar

/build/$dest_appdir/R --vanilla -e 'BiocManager::install("PCAtools")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("pcaExplorer")'
/build/$dest_appdir/R --vanilla -e $'BiocManager::install("Rhtslib")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("Rsamtools")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("ggbio")'

/build/$dest_appdir/R --vanilla -e 'BiocManager::install("bambu")'

# leafcutter (needs rstan and rstantools)
yuminstall gsl-devel
cp -a /usr/lib64/libgsl*.so* /build/$dest_appdir/lib64/R/lib
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("DirichletMultinomial")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("DropletUtils")'
/build/$dest_appdir/R --vanilla -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter", version = "0.2.9")'

# other
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("MOFA2")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("mixOmics")'
# Rhtslib does not compile for some reason, so not for now
# /build/$dest_appdir/R --vanilla -e 'BiocManager::install("chromVAR")'
#/build/$dest_appdir/R --vanilla -e 'BiocManager::install("OmicCircos")'
/build/$dest_appdir/R --vanilla -e 'install.packages("https://github.com/dzhang32/ggtranscript/archive/refs/tags/v0.99.3.tar.gz", repos=NULL, type="source")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("scDblFinder")'
/build/$dest_appdir/R --vanilla -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", version = "2.0.3")'
/build/$dest_appdir/R --vanilla -e 'BiocManager::install("cleaver")'

# hdf5r
cd /build/R-$Rversion
yuminstall hdf5-devel
wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_14_1/src/hdf5-1.14.1-2.tar.gz
tar xvzf hdf5-1.14.1-2.tar.gz
cd /build/R-$Rversion/hdf5-1.14.1-2
./configure
make
make install
sudo cp -ra hdf5/bin/* /usr/bin
sudo cp -ra hdf5/lib/* /usr/lib
sudo cp -ra hdf5/share/* /usr/share/
sudo cp -ra hdf5/include/* /usr/include/
cp -ra hdf5/lib/* /build/$dest_appdir/lib64/R/lib
/build/$dest_appdir/R --vanilla -e 'require(devtools); install_version("hdf5r", version = "1.3.8", dependencies=FALSE, repos = "http://cran.us.r-project.org")'
cp -a -f hdf5/lib/* /build/$dest_appdir/lib64/R/lib
cp -a -f hdf5/bin/* /build/$dest_appdir/lib64/R/bin

# "install" sc-type
cd /build/$dest_appdir/lib64/R/library
wget https://github.com/IanevskiAleksandr/sc-type/archive/refs/tags/v1.0.tar.gz
tar xvzf v1.0.tar.gz
rm v1.0.tar.gz

# make links copy endproduct
cd /build
ln -sf $dest_appdir/R dirR-$dirRversion
ln -sf $dest_appdir/Rscript dirRscript-$dirRversion
ln -sf $dest_appdir/R dirR
ln -sf $dest_appdir/Rscript dirRscript
ln -sf $dest_appdir/R R
ln -sf $dest_appdir/Rscript Rscript
rm $dest_appdir.tar.gz || true
tar cvzf $dest_appdir.tar.gz $dest_appdir dirR-$dirRversion dirR R dirRscript-$dirRversion dirRscript Rscript
rm -rf /io/extra$ARCH/$dest_appdir || true
mkdir /io/extra$ARCH || true
cp -raf $dest_appdir dirR-$dirRversion dirR R /io/extra$ARCH

echo "Finished building dirR4 in $builddir/$dest_appdir"
echo "Copied to $srcdir/extra$ARCH"

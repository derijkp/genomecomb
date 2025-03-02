# This script is meant to be sourced in a Holy build box build script to start the HBB docker
# it interprets the options:
# -a|-arch|--arch: 64, x86_64 or linux-x86_64 for 64 bit Linux build (default); ix86, 32 or linux-ix86 for 32 bits Linux build; win, windows-x86_64 or mingw-w64 for Windows 64 bit build
# -b|-bits|--bits: select 32 or 64 bits Linux build (default 64 = same as -arch x86_64)
# -d|-builddir|--builddir: top directory to build in (default ~/build/bin-$arch)
# and presents all options (excluding -bbuilddir) to the program in the docker

# you can use the following to use in the build script
## Prepare and start docker with Holy Build box
## ============================================
#script="$(readlink -f "$0")"
#dir="$(dirname "$script")"
#source "${dir}/hbb_start.sh"

# The Holy build box environment requires docker, make sure it is installed
# e.g. on ubuntu and derivatives
# sudo apt install docker.io
# Also make sure you have permission to use docker
# sudo usermod -a -G docker $USER

# Prepare and start docker with Holy Build box
# ============================================

# check if we are already in holy build box
if [ ! -f /hbb_exe/activate ]; then
	# find directory of script
	script="$(readlink -f "$0")"
	dir="$(dirname "$script")"
	file="$(basename "$script")"
	if [ $(basename "$dir") = "build" ] ; then
		srcdir="$(dirname "$dir")"
		file="build/$file"
	else
		srcdir="$dir"
	fi

	echo "Running script $dir/$file"
	builddir=""
	arch=linux-x86_64
	argumentspos=1; 
	while [[ "$#" -gt 0 ]]; do case $1 in
		-b|-bits|--bits) arch="$2"; shift;;
		-a|-arch|--arch) arch="$2"; shift;;
		-d|-builddir|--builddir) builddir="$(readlink -f "$2")" ; shift;;
		*) arguments[$argumentspos]="$1"; argumentspos+=1 ; arguments[$argumentspos]="$2"; argumentspos+=1 ; shift;;
	esac; shift; done
	if [ "$bits" = "32" ] ; then arch=linux-ix86 ; fi
	if [ "$bits" = "64" ] ; then arch=linux-x86_64 ; fi
	if [ "$arch" = "32" ] || [ "$arch" = "ix86" ] ; then 
		arch=linux-ix86
	elif [ "$arch" = "64" ] || [ "$arch" = "x86_64" ]; then
		arch=linux-x86_64
	elif [ "$arch" = "win-ix86" ] || [ "$arch" = "mingw-w32" ]; then
		arch=windows-ix86
	elif [ "$arch" = "win" ] || [ "$arch" = "windows" ] || [ "$arch" = "win-x86_64" ] || [ "$arch" = "mingw-w64" ]; then
		arch=windows-x86_64
	fi
	if [ "$builddir" = "" ] ; then
		if [ "$arch" = "linux-x86_64" ] ; 	then
			builddir="$HOME/build/bin-x86_64"
		elif [ "$arch" = "linux-ix86" ] ; 	then
			builddir="$HOME/build/bin-ix86"
		elif [ "$arch" = "windows-x86_64" ] ; 	then
			builddir="$HOME/build/bin-windows-x86_64"
		elif [ "$arch" = "windows-ix86" ] ; 	then
			builddir="$HOME/build/bin-windows-ix86"
		else
			echo "unknown arch $arch"
			exit 1
		fi
	fi
	mkdir -p "$builddir"
	echo "Build $arch version"
	echo "builddir=$builddir"
	echo "srcdir=$srcdir"
	# run the script in holy build box
	uid=$(id -u)
	gid=$(id -g $uid)
	
	if [ "$arch" = "linux-ix86" ] ; 	then
		if docker image list | grep --quiet 'hbb32.*2.2.0'; then
			buildbox=hbb32:2.2.0
		else
			buildbox=phusion/holy-build-box-32:2.2.0
		fi
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$builddir:/build" "$buildbox" linux32 bash "/io/$file" "stage2" "$file" "$arch" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	else
		if docker image list | grep --quiet 'hbb64.*2.2.0'; then
			buildbox=hbb64:2.2.0
		else
			buildbox=phusion/holy-build-box-64:2.2.0
		fi
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$builddir:/build" "$buildbox" bash "/io/$file" "stage2" "$file" "$arch" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	fi
	exit
fi

if [ "$1" = "stage2" ] ; then
	# in stage 2 we will create the user build with sudo access
	# then restart the script (skipping to stage 3)
	file=$2
	arch=$3
	uid=$4
	gid=$5
	srcdir=$6
	builddir=$7
	# prepare the user build with sudo rights
	echo "installing sudo ($arch)"
	# to stop "checksum is invalid" errors when using yum in 32 bit docker
	if [ "$arch" = linux-ix86 ] ; then
		rm -f /etc/yum.repos.d/phusion_centos-6-scl-i386.repo
		echo "change repos to vault"
		curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo
		curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo
		if ! rpm --quiet --query yum-plugin-ovl; then
			yum install -q -y yum-plugin-ovl
		fi
	fi
	if ! rpm --quiet --query sudo; then
		yum install -q -y sudo
	fi
	echo "preparing user build with uid=$uid and gid=$gid"
	groupadd build --gid $gid
	useradd build --uid $uid --gid $gid -s /usr/bin/bash
	# usermod -a -G wheel build
	echo "build ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-build
	# default nr of processes (for user build) is sometimes not enough
	sudo sed -i 's/1024/10240/' /etc/security/limits.d/90-nproc.conf
	# (re)start script for stage 3: running the actual code
	sudo -u build bash /io/$file "stage3" ${@:2}
	exit
fi

# stage 3: run the actual script (first do some settings)

function yuminstall {
	echo "yuminstall $1"
	if ! rpm --quiet --query "$1"; then
		sudo yum install -y "$1"
	fi
}

# centos 6 is EOL, moved to vault: adapt the repos
if [ "$arch" = "linux-ix86" ] ; 	then

rm -f /etc/yum.repos.d/phusion_centos-6-scl-i386.repo
if ! cat /etc/yum.repos.d/CentOS-Base.repo | grep --quiet vault; then
	echo "change repos to vault"
	curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo
	curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo
fi


else

if ! cat /etc/yum.repos.d/CentOS-Base.repo | grep --quiet vault; then
	echo "change repos to vault"
	sudo curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo
	sudo curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo
	sudo curl https://www.getpagespeed.com/files/centos6-scl-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl.repo
	sudo curl https://www.getpagespeed.com/files/centos6-scl-rh-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl-rh.repo
fi

fi
# endof: centos 6 is EOL, moved to vault: adapt the repos

file=$2
# install yuminstall and env vars in .bashrc so it will be available if the new shell is started
mkdir -p /home/build

echo "
file=$2
arch=$3
echo arch_sh=$arch
uid=$4
gid=$5
srcdir=$6
builddir=$7
" > /home/build/.bashrc

echo 'if [ "$arch" = 'linux-ix86' ] ; then
	ARCH='-linux-ix86'
	arch=linux-ix86
	bits=32
elif [ "$arch" = 'windows-ix86' ] ; then
	ARCH='-windows-ix86'
	arch=windows-ix86
	bits=64
	sudo yum install -y mingw64-gcc mingw64-zlib mingw64-zlib-static
	sudo yum install -y mingw32-gcc mingw32-zlib mingw32-zlib-static
	# sudo yum install -y mingw64-g++ mingw64-libgnurx-static mingw64-boost mingw64-boost-static
	export CROSSTARGET="w64-mingw32"
	export TARGETARCHITECTURE="i686"
	#
	export TARGET="$TARGETARCHITECTURE-$CROSSTARGET"
	export HOST="$TARGET"
	export CROSSBASE="/usr/lib/gcc/$TARGET/4.9.2/"
	# export BASE=/build/deps-$TARGET
	export LDFLAGS="-L$CROSSBASE"
	export CPPFLAGS="-I$CROSSBASE/include"
	export CFLAGS="-I$CROSSBASE/include"
	export CROSS_COMPILE="$TARGET-"
	#export CC=${CROSS_COMPILE}gcc
	#export CXX=${CROSS_COMPILE}g++
	#export CPP=${CROSS_COMPILE}cpp
	#export AR=${CROSS_COMPILE}ar
	#export LD=${CROSS_COMPILE}ld
	#export RANLIB=${CROSS_COMPILE}ranlib
	export DIRTCL=/build/dirtcl-$TARGET
elif [ "$arch" = 'windows-x86_64' ] ; then
	ARCH='-windows-x86_64'
	arch=windows-x86_64
	bits=64
	sudo yum install -y mingw64-gcc mingw64-zlib mingw64-zlib-static
	# sudo yum install -y mingw64-g++ mingw64-libgnurx-static mingw64-boost mingw64-boost-static
	export CROSSTARGET="w64-mingw32"
	export TARGETARCHITECTURE="x86_64"
	#
	export TARGET="$TARGETARCHITECTURE-$CROSSTARGET"
	export HOST="$TARGET"
	export CROSSBASE="/usr/lib/gcc/$TARGET/4.9.2/"
	# export BASE=/build/deps-$TARGET
	export LDFLAGS="-L$CROSSBASE"
	export CPPFLAGS="-I$CROSSBASE/include"
	export CFLAGS="-I$CROSSBASE/include"
	export CROSS_COMPILE="$TARGET-"
	#export CC=${CROSS_COMPILE}gcc
	#export CXX=${CROSS_COMPILE}g++
	#export CPP=${CROSS_COMPILE}cpp
	#export AR=${CROSS_COMPILE}ar
	#export LD=${CROSS_COMPILE}ld
	#export RANLIB=${CROSS_COMPILE}ranlib
	export DIRTCL=/build/dirtcl-$TARGET
else
	ARCH=''
	arch=linux-x86_64
	bits=64
fi
function yuminstall {
	echo "yuminstall $1"
	if ! rpm --quiet --query "$1"; then
		sudo yum install -y "$1"
	fi
}
' >> /home/build/.bashrc
shift 7;

if [ $(basename "$file") = "start_hbb.sh" ] ; then
	# if run as start_hbb.sh directly, show a shell
	echo "shell started by start_hbb.sh"
	bash
	exit
fi

# if sourced in another script, continue executing this other script
# after sourcing the settings put into the bashrc
source /home/build/.bashrc

echo "Entering Holy Build Box environment; building using arch $arch, $bits bits, uid=$uid, gid=$gid"

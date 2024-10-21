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
	elif [ "$arch" = "win" ] || [ "$arch" = "mingw-w64" ]; then
		arch=windows-x86_64
	fi
	if [ "$builddir" = "" ] ; then
		if [ "$arch" = "linux-ix86" ] ; 	then
			builddir="$HOME/build/bin-ix86"
		else
			builddir="$HOME/build/bin-x86_64"
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
		if docker image list | grep --quiet hbb3-32; then
			buildbox=hbb3-32
		else
			buildbox=phusion/holy-build-box-32:3.0.2
		fi
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$builddir:/build" "$buildbox" linux32 bash "/io/$file" "stage2" "$file" "$bits" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	else
		if docker image list | grep --quiet hbb3-64; then
			buildbox=hbb3-64
		else
			buildbox=phusion/holy-build-box-64:3.0.2
		fi
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$builddir:/build" "$buildbox" bash "/io/$file" "stage2" "$file" "$bits" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	fi
	exit
fi

# centos 7 is EOL, moved to vault: adapt the repos

if ! cat /etc/yum.repos.d/CentOS-Base.repo | grep --quiet vault; then
	echo "change repos to vault"

echo '# /etc/yum.repos.d/CentOS-Base.repo
[base]
name=CentOS-$releasever - Base
baseurl=http://vault.centos.org/7.9.2009/os/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
[updates]
name=CentOS-$releasever - Updates
baseurl=http://vault.centos.org/7.9.2009/updates/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
[extras]
name=CentOS-$releasever - Extras
baseurl=http://vault.centos.org/7.9.2009/extras/$basearch/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
[centosplus]
name=CentOS-$releasever - Plus
baseurl=http://vault.centos.org/7.9.2009/centosplus/$basearch/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
[contrib]
name=CentOS-$releasever - Contrib
baseurl=http://vault.centos.org/7.9.2009/contrib/$basearch/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
' > /etc/yum.repos.d/CentOS-Base.repo

echo '[epel]
name=Extra Packages for Enterprise Linux 7 - $basearch
#baseurl=http://download.fedoraproject.org/pub/epel/7/$basearch
baseurl=https://archives.fedoraproject.org/pub/archive/epel/7/x86_64
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-7&arch=$basearch
failovermethod=priority
enabled=1
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7
[epel-debuginfo]
name=Extra Packages for Enterprise Linux 7 - $basearch - Debug
#baseurl=http://download.fedoraproject.org/pub/epel/7/$basearch/debug
mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-debug-7&arch=$basearch
failovermethod=priority
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7
gpgcheck=1
[epel-source]
name=Extra Packages for Enterprise Linux 7 - $basearch - Source
#baseurl=http://download.fedoraproject.org/pub/epel/7/SRPMS
baseurl=https://archives.fedoraproject.org/pub/archive/epel/7/SRPMS
#mirrorlist=https://mirrors.fedoraproject.org/metalink?repo=epel-source-7&arch=$basearch
failovermethod=priority
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-EPEL-7
gpgcheck=1
' > /etc/yum.repos.d/epel.repo

echo '# CentOS-SCLo-sclo.repo
#
# Please see http://wiki.centos.org/SpecialInterestGroup/SCLo for more
# information
[centos-sclo-sclo]
name=CentOS-7 - SCLo sclo
baseurl=http://vault.epel.cloud/centos/7/sclo/$basearch/sclo/
gpgcheck=1
enabled=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
[centos-sclo-sclo-testing]
name=CentOS-7 - SCLo sclo Testing
baseurl=http://buildlogs.centos.org/centos/7/sclo/$basearch/sclo/
gpgcheck=0
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
[centos-sclo-sclo-source]
name=CentOS-7 - SCLo sclo Sources
baseurl=http://vault.epel.cloud/centos/7/sclo/Source/sclo/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
[centos-sclo-sclo-debuginfo]
name=CentOS-7 - SCLo sclo Debuginfo
baseurl=http://debuginfo.centos.org/centos/7/sclo/$basearch/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
' > /etc/yum.repos.d/CentOS-SCLo-scl.repo

echo '# CentOS-SCLo-rh.repo
#
# Please see http://wiki.centos.org/SpecialInterestGroup/SCLo for more
# information

[centos-sclo-rh]
name=CentOS-7 - SCLo rh
baseurl=http://vault.epel.cloud/centos/7/sclo/$basearch/rh/
gpgcheck=1
enabled=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo

[centos-sclo-rh-testing]
name=CentOS-7 - SCLo rh Testing
baseurl=http://buildlogs.centos.org/centos/7/sclo/$basearch/rh/
gpgcheck=0
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo

[centos-sclo-rh-source]
name=CentOS-7 - SCLo rh Sources
baseurl=http://vault.epel.cloud/centos/7/sclo/Source/rh/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo

[centos-sclo-rh-debuginfo]
name=CentOS-7 - SCLo rh Debuginfo
baseurl=http://debuginfo.centos.org/centos/7/sclo/$basearch/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-SIG-SCLo
' > /etc/yum.repos.d/CentOS-SCLo-scl-rh.repo

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
	if ! rpm --quiet --query sudo; then
		yum install -q -y sudo
	fi
	echo "preparing user build with uid=$uid and gid=$gid"
	groupadd build --gid $gid
	useradd build --uid $uid --gid $gid
	# usermod -a -G wheel build
	echo "build ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-build
	# default nr of processes (for user build) is sometimes not enough
	sudo sed -i 's/4096/10240/' /etc/security/limits.d/20-nproc.conf
	# (re)start script for stage 3: running the actual code
	sudo -u build bash /io/$file "stage3" "$file" "$bits" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	exit
fi

# stage 3: run the actual script (first do some settings)

function yuminstall {
	echo "yuminstall $1"
	if ! rpm --quiet --query "$1"; then
		sudo yum install -y "$1"
	fi
}

file=$2

if [ $(basename "$file") = "start_hbb3.sh" ] ; then
	# install yuminstall in .bashrc so it will be available in the new shell started here
	mkdir -p /home/build
	echo "
	file=$2
	arch=$3
	uid=$4
	gid=$5
	srcdir=$6
	builddir=$7
	" >> /home/build/.bashrc
	echo 'if [ "$arch" = 'linux-ix86' ] ; then
		ARCH='-linux-ix86'
		arch=linux-ix86
		bits=32
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
	# if run as start_hbb.sh directly, show a shell
	echo "shell started by start_hbb3.sh"
	bash
	exit
fi

# if sourced in another script, continue executing this other script
if [ "$3" = 'linux-ix86' ] ; then
	ARCH='-linux-ix86'
	arch=linux-ix86
	bits=32
else
	ARCH=''
	arch=linux-x86_64
	bits=64
fi
uid=$4;
gid=$5;
srcdir=$6;
builddir=$7;
shift 7;

echo "Entering Holy Build Box environment; building using arch $arch, uid=$uid, gid=$gid"

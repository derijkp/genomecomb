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
		if docker image list | grep --quiet hbb32; then
			buildbox=hbb32
		else
			buildbox=phusion/holy-build-box-32:2.2.0
		fi
		docker run --net=host -t -i --rm -v "$srcdir:/io" -v "$builddir:/build" "$buildbox" linux32 bash "/io/$file" "stage2" "$file" "$arch" "$uid" "$gid" "$srcdir" "$builddir" ${arguments[*]}
	else
		if docker image list | grep --quiet hbb64; then
			buildbox=hbb64
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
		if ! rpm --quiet --query yum-plugin-ovl; then
			yum install -q -y yum-plugin-ovl
		fi
	fi
	if ! rpm --quiet --query sudo; then
		yum install -q -y sudo
	fi
	echo "preparing user build with uid=$uid and gid=$gid"
	groupadd build --gid $gid
	useradd build --uid $uid --gid $gid
	# usermod -a -G wheel build
	echo "build ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-build
	# default nr of processes (for user build) is sometimes not enough
	sudo sed -i 's/1024/10240/' /etc/security/limits.d/90-nproc.conf
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

# centos 6 is EOL, moved to vault: adapt the repos
if [ "$arch" = "linux-ix86" ] ; 	then

cd
echo '
[base]
name=CentOS-$releasever - Base
#mirrorlist=http://mirrorlist.centos.org/?release=$releasever&arch=i386&repo=os&infra=$infra
baseurl=http://vault.centos.org/centos/$releasever/os/i386/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-6

#released updates
[updates]
name=CentOS-$releasever - Updates
#mirrorlist=http://mirrorlist.centos.org/?release=$releasever&arch=i386&repo=updates&infra=$infra
baseurl=http://vault.centos.org/centos/$releasever/updates/i386/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-6

#additional packages that may be useful
[extras]
name=CentOS-$releasever - Extras
#mirrorlist=http://mirrorlist.centos.org/?release=$releasever&arch=i386&repo=extras&infra=$infra
baseurl=http://vault.centos.org/centos/$releasever/extras/i386/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-6

#additional packages that extend functionality of existing packages
[centosplus]
name=CentOS-$releasever - Plus
#mirrorlist=http://mirrorlist.centos.org/?release=$releasever&arch=i386&repo=centosplus&infra=$infra
baseurl=http://vault.centos.org/centos/$releasever/centosplus/i386/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-6

#contrib - packages by Centos Users
[contrib]
name=CentOS-$releasever - Contrib
#mirrorlist=http://mirrorlist.centos.org/?release=$releasever&arch=i386&repo=contrib&infra=$infra
baseurl=http://vault.centos.org/centos/$releasever/contrib/i386/
gpgcheck=1
enabled=0
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-6
' > temp
sudo mv -f temp /etc/yum.repos.d/CentOS-Base.repo

else

if ! cat /etc/yum.repos.d/CentOS-Base.repo | grep --quiet vault; then
	echo "change repos to vault"
	sudo curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo
	sudo curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo
	sudo curl https://www.getpagespeed.com/files/centos6-scl-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl.repo
	sudo curl https://www.getpagespeed.com/files/centos6-scl-rh-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl-rh.repo
fi

fi

file=$2

if [ $(basename "$file") = "start_hbb.sh" ] ; then
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
	echo "shell started by start_hbb.sh"
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

echo "Entering Holy Build Box environment; building using arch $arch, $bits bits, uid=$uid, gid=$gid"

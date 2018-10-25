#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Installs the tools used in CI/CD (tested for Debian 9.3, CentOS 7, Amazon Linux 2, Amazon Linux AMI).
##
EDGE_CURRENT_DIR=$(pwd)
EDGE_TMP_DIR=$(mktemp -d)
EDGE_N_BUILD_PROC=$(cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}')

EDGE_DIST=$(cat /etc/*-release | grep PRETTY_NAME)

cd ${EDGE_TMP_DIR}

###############
# Basic Tools #
###############
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  # recent build tools
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y gcc
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y gfortran
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y libiomp5 libiomp-dev
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y software-properties-common
  # other
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y wget
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y unzip
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y m4
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y dh-autoreconf
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y make
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y cmake
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y scons
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y git
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y libxml2-utils
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y python-pip python3-pip
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y openmpi-bin libopenmpi-dev
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y cppcheck
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y gmsh
elif [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  # recent build tools
  sudo yum install -y -q -e 0 centos-release-scl
  sudo yum install -y -q -e 0 devtoolset-7
  source /opt/rh/devtoolset-7/enable
  echo "source /opt/rh/devtoolset-7/enable" | sudo tee --append /etc/bashrc
  # other
  sudo yum install -y -q -e 0 wget
  sudo yum install -y -q -e 0 unzip
  sudo yum install -y -q -e 0 m4
  sudo yum install -y -q -e 0 dh-autoreconf
  sudo yum install -y -q -e 0 make
  sudo yum install -y -q -e 0 cmake
  sudo yum install -y -q -e 0 scons
  sudo yum install -y -q -e 0 git
  sudo yum install -y -q -e 0 libxml2-python.x86_64
  sudo yum install -y -q -e 0 python python34 python-devel python34-devel python-setuptools python34-setuptools
  sudo yum install -y -q -e 0 openmpi openmpi-devel
  echo "module load mpi" | sudo tee --append /etc/bashrc
  sudo yum install -y -q -e 0 cppcheck
  # TODO: no gmsh RPM available, move to custom install
elif [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]]
then
  sudo yum groupinstall -y -q -e 0 "Development Tools"
  sudo yum install -y -q -e 0 cmake
  sudo yum install -y -q -e 0 python python-pip python3 python3-pip
  sudo pip install -q --egg scons
  sudo yum install -y -q -e 0 openmpi openmpi-devel
  echo "export PATH=/usr/lib64/openmpi/bin/:$PATH" | sudo tee --append /etc/bashrc
  # TODO: no cppcheck RPM available
  # TODO: no gmsh RPM available
elif [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  # install custom gcc, as the OS doesnt ship it with OpenMP support
  svn co svn://gcc.gnu.org/svn/gcc/tags/gcc_8_2_0_release/ gcc > /dev/null
  cd gcc
  ./contrib/download_prerequisites
  ./configure --enable-languages=c,c++,fortran --disable-multilib > /dev/null
  make -j ${EDGE_N_BUILD_PROC} > /dev/null
  sudo make install > /dev/null
  cd ..
  rm -rf gcc

  # link gcc-8
  sudo unlink /usr/bin/gcc
  sudo unlink /usr/bin/g++
  sudo unlink /usr/bin/gfortran
  sudo ln -s /usr/local/bin/gcc /usr/bin
  sudo ln -s /usr/local/bin/g++ /usr/bin
  sudo ln -s /usr/local/bin/gfortran /usr/bin

  sudo yum install -y -q -e 0 python python34 python-devel python34-devel python-setuptools python34-setuptools
  sudo ln -s /usr/local/bin/easy_install* /bin
  sudo yum install -y -q -e 0 scons
  echo "module load mpi" | sudo tee --append /etc/bashrc
  # TODO: no cppcheck RPM available
  # TODO: no gmsh RPM available
fi

#########
# Clang #
#########
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  echo "deb http://apt.llvm.org/stretch/ llvm-toolchain-stretch-7 main"     | sudo tee /etc/apt/sources.list.d/clang.list
  echo "deb-src http://apt.llvm.org/stretch/ llvm-toolchain-stretch-7 main" | sudo tee --append /etc/apt/sources.list.d/clang.list
  sudo wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key|sudo apt-key add -
  sudo apt-get update
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y clang-7 lldb-7 lld-7
  sudo ln -s /usr/bin/clang-7   /usr/bin/clang
  sudo ln -s /usr/bin/clang++-7 /usr/bin/clang++
elif [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo yum install -y -q -e 0 llvm-toolset-7-clang
  sudo yum install -y -q -e 0 clang
  source /opt/rh/llvm-toolset-7/enable
  echo "source /opt/rh/llvm-toolset-7/enable" | sudo tee --append /etc/bashrc
elif [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]]
then
  echo "" > /dev/null # TODO: fix Clang support
elif [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  sudo yum install -y -q clang clang-develop
fi

############
# Valgrind #
############
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 libc6-dbg
  wget http://mirrors.kernel.org/sourceware/valgrind/valgrind-3.13.0.tar.bz2
  tar -xjf valgrind-3.13.0.tar.bz2
  cd valgrind-3.13.0
  ./configure
  make -j ${EDGE_N_BUILD_PROC}
  sudo make install > /dev/null
  cd ..
elif [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  echo "" > /dev/null # valgrind is already part of the CentOS tools
elif [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]]
then
  sudo yum install -y -q -e 0 valgrind
fi

###########
# Git LFS #
###########
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
  sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y git-lfs
elif [[ ${EDGE_DIST} == *"CentOS"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.rpm.sh | sudo bash
  sudo yum install -y -q -e 0 git-lfs
fi

##################
# Python modules #
##################
if [[ ${EDGE_DIST} == *"CentOS"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  sudo easy_install-3.4 pip
  sudo easy_install pip
  sudo ln -s /usr/local/bin/pip* /bin
  sudo pip install --upgrade pip setuptools
  sudo pip3 install --upgrade pip setuptools
fi
sudo pip -q install xmltodict matplotlib netCDF4 h5py
sudo pip3 -q install xmltodict matplotlib h5py

###########
# Vagrant #
###########
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  wget https://releases.hashicorp.com/vagrant/2.1.2/vagrant_2.1.2_x86_64.deb -O vagrant.deb
  sudo dpkg -i vagrant.deb
elif [[ ${EDGE_DIST} == *"CentOS"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  wget https://releases.hashicorp.com/vagrant/2.2.0/vagrant_2.2.0_x86_64.rpm -O vagrant.rpm
  sudo yum install -y -q -e 0 vagrant.rpm
fi

############
# Clean up #
############
if [[ ${EDGE_DIST} == *"Debian"* ]]
then
  sudo apt-get clean
elif [[ ${EDGE_DIST} == *"CentOS"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux 2"* ]] || [[ ${EDGE_DIST} == *"Amazon Linux AMI"* ]]
then
  sudo yum clean all
  sudo rm -rf /var/cache/yum
fi

sudo rm -rf ${EDGE_TMP_DIR}
cd ${EDGE_CURRENT_DIR}
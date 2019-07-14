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
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  # recent build tools
  sudo yum install -y -q -e 0 centos-release-scl
  sudo yum install -y -q -e 0 devtoolset-7
  sudo yum install -y -q -e 0 devtoolset-7-libasan-devel devtoolset-7-libubsan-devel
  source /opt/rh/devtoolset-7/enable
  echo "source /opt/rh/devtoolset-7/enable > /dev/null" | sudo tee --append /etc/bashrc
  # other
  sudo yum install -y -q -e 0 wget
  sudo yum install -y -q -e 0 unzip
  sudo yum install -y -q -e 0 m4
  sudo yum install -y -q -e 0 dh-autoreconf
  sudo yum install -y -q -e 0 make
  sudo yum install -y -q -e 0 cmake3
  sudo ln -s /usr/bin/cmake3 /usr/bin/cmake
  sudo yum install -y -q -e 0 git
  sudo yum install -y -q -e 0 libxml2-python.x86_64
  sudo yum install -y -q -e 0 python python34 python-devel python34-devel python-setuptools python34-setuptools python-pip python34-pip
  sudo yum install -y -q -e 0 cppcheck
  sudo yum install -y -q -e 0 glibc-static
  sudo yum install -y -q -e 0 libhugetlbfs libhugetlbfs-devel libhugetlbfs-utils
  sudo yum install -y -q -e 0 irqbalance
  sudo yum install -y -q -e 0 mesa-libGLU libXcursor libXft libXinerama
  sudo yum install -y -q -e 0 atop
  sudo yum install -y -q -e 0 GMT
  sudo yum install -y -q -e 0 ghostscript
  # UCVM and EDGEv dependencies
  sudo yum install -y -q -e 0 proj-static proj-devel
  # EDGEcut dependencies
  sudo yum install -y -q -e 0 gmp-devel mpfr-devel boost-devel
  sudo pip install meshio > /dev/null
  # GoCD dependencies
  sudo yum install -y -q -e 0 java-latest-openjdk
fi

########
# GMSH #
########
wget http://gmsh.info/bin/Linux/gmsh-4.1.5-Linux64.tgz -O gmsh.tgz
sudo tar -xf gmsh.tgz -C /usr --strip-components=1

#########
# Clang #
#########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo yum install -y -q -e 0 llvm-toolset-7-clang llvm-toolset-7-libomp llvm-toolset-7-libomp-devel
  sudo yum install -y -q -e 0 clang
  source /opt/rh/llvm-toolset-7/enable
  echo "source /opt/rh/llvm-toolset-7/enable > /dev/null" | sudo tee --append /etc/bashrc
fi

###########
# OpenMPI #
###########
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.0.tar.bz2 -O openmpi.tar.bz2
mkdir openmpi
tar -xjf openmpi.tar.bz2 -C openmpi --strip-components=1
cd openmpi
./configure --enable-static=yes --enable-mpi1-compatibility > /dev/null
make -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo make install > /dev/null
cd ..

############################
# Intel Parallel Studio XE #
############################
if [ -f ${EDGE_CURRENT_DIR}/parallel_studio_xe_*.tgz ]
then
  tar -xzf ${EDGE_CURRENT_DIR}/parallel_studio_xe_*.tgz
  cd parallel_studio_xe_*
  sed -i 's/ACCEPT_EULA=decline/ACCEPT_EULA=accept/' ./silent.cfg
  sed -i 's/ACTIVATION_TYPE=exist_lic/ACTIVATION_TYPE=license_file/' ./silent.cfg
  sed -i "s<#ACTIVATION_LICENSE_FILE=<ACTIVATION_LICENSE_FILE=$(ls ${EDGE_CURRENT_DIR}/*.lic)<" ./silent.cfg
  sed -i "s/ARCH_SELECTED=ALL/ARCH_SELECTED=INTEL64/" ./silent.cfg
  sudo ./install.sh -s silent.cfg
  cd ..
fi

###########
# Git LFS #
###########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.rpm.sh | sudo bash
  sudo yum install -y -q -e 0 git-lfs
fi

##################
# Python modules #
##################
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo pip install pip==18.0
  sudo pip install --upgrade setuptools
  sudo pip3 install pip==18.0
  sudo pip install --upgrade setuptools
fi
sudo pip3 -q install xmltodict matplotlib netCDF4 h5py


#########
# Scons #
#########
wget http://prdownloads.sourceforge.net/scons/scons-3.0.1.tar.gz -O scons.tar.gz
mkdir scons
tar -xzf scons.tar.gz -C scons --strip-components=1
cd scons
sudo python setup.py install
cd ..

###########
# Vagrant #
###########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  wget https://releases.hashicorp.com/vagrant/2.2.0/vagrant_2.2.0_x86_64.rpm -O vagrant.rpm
  sudo yum install -y -q -e 0 vagrant.rpm
fi

#################
# GCP specifics #
#################
if [[ $(curl metadata.google.internal -si | grep Google) ]]
then
  sudo python -m pip install google-api-python-client
fi

############
# Clean up #
############
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo yum clean all
  sudo rm -rf /var/cache/yum
fi

sudo rm -rf ${EDGE_TMP_DIR}
cd ${EDGE_CURRENT_DIR}

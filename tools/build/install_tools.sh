#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2019-2020, Alexander Breuer
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
  # install extra packages
  sudo dnf install -y -q -e 0 epel-release
  # recent build tools
  sudo dnf install -y -q -e 0 gcc
  sudo dnf install -y -q -e 0 libasan libubsan
  sudo dnf install -y -q -e 0 gcc-toolset-9
  #source /opt/rh/gcc-toolset-9/enable
  #echo "source /opt/rh/gcc-toolset-9/enable > /dev/null" | sudo tee --append /etc/bashrc
  # other
  sudo dnf install -y -q -e 0 hostname
  sudo dnf install -y -q -e 0 wget
  sudo dnf install -y -q -e 0 unzip
  sudo dnf install -y -q -e 0 bzip2
  sudo dnf install -y -q -e 0 m4
  sudo dnf install -y -q -e 0 autoconf
  sudo dnf install -y -q -e 0 make
  sudo dnf install -y -q -e 0 cmake3
  sudo dnf install -y -q -e 0 git
  sudo dnf install -y -q -e 0 git-lfs
  sudo dnf install -y -q -e 0 python3-libxml2
  sudo dnf install -y -q -e 0 python2 python3 python2-devel python3-devel python2-setuptools python3-setuptools python2-pip python3-pip
  sudo ln -s /usr/bin/python3 /usr/bin/python
  sudo ln -s /usr/bin/pip3 /usr/bin/pip
  sudo dnf install -y -q -e 0 cppcheck
  sudo dnf install -y -q -e 0 valgrind
  sudo dnf install -y -q -e 0 glibc-static
  sudo dnf install -y -q -e 0 libhugetlbfs libhugetlbfs-devel libhugetlbfs-utils
  sudo dnf install -y -q -e 0 irqbalance
  sudo dnf install -y -q -e 0 mesa-libGLU libXcursor libXft libXinerama
  sudo dnf install -y -q -e 0 atop
  sudo dnf install -y -q -e 0 htop
  sudo dnf install -y -q -e 0 yum-plugin-copr
  sudo dnf copr enable -y -q -e 0 genericmappingtools/gmt
  sudo dnf install -y -q -e 0 gmt
  sudo dnf install -y -q -e 0 ghostscript
  # UCVM and EDGEv dependencies
  sudo dnf install -y -q -e 0 metis metis-devel
  sudo dnf install -y -q -e 0 proj-static proj-devel
  # EDGEcut dependencies
  sudo pip install meshio > /dev/null
  # Pipeline dependencies
  sudo dnf install -y -q -e 0 npm
fi

########
# GMSH #
########
wget http://gmsh.info/bin/Linux/gmsh-4.5.6-Linux64.tgz -O gmsh.tgz
sudo tar -xf gmsh.tgz -C /usr --strip-components=1

#########
# Clang #
#########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo dnf install -y -q -e 0 llvm-toolset
fi

###########
# OpenMPI #
###########
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz -O openmpi.tar.gz
mkdir openmpi
tar -xf openmpi.tar.gz -C openmpi --strip-components=1
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

##################
# Python modules #
##################
sudo pip -q install xmltodict matplotlib netCDF4 h5py obspy

#########
# Scons #
#########
wget http://prdownloads.sourceforge.net/scons/scons-3.1.2.tar.gz -O scons.tar.gz
mkdir scons
tar -xzf scons.tar.gz -C scons --strip-components=1
cd scons
sudo python setup.py install
cd ..

#############
# Buildkite #
#############
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  sudo sh -c 'echo -e "[buildkite-agent]\nname = Buildkite Pty Ltd\nbaseurl = https://yum.buildkite.com/buildkite-agent/stable/x86_64/\nenabled=1\ngpgcheck=0\npriority=1" > /etc/yum.repos.d/buildkite-agent.repo'
  sudo dnf install -y -q -e 0 buildkite-agent
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
  sudo dnf clean all
  sudo rm -rf /var/cache/dnf
fi

sudo rm -rf ${EDGE_TMP_DIR}
cd ${EDGE_CURRENT_DIR}

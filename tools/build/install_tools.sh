#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2021, Friedrich Schiller University Jena
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
# Installs the tools used in CI/CD (tested for Fedora 33 and Ubuntu Focal).
##
EDGE_CURRENT_DIR=$(pwd)
EDGE_TMP_DIR=$(mktemp -d)
EDGE_N_BUILD_PROC=$(grep -c ^processor /proc/cpuinfo)

EDGE_DIST=$(cat /etc/*-release | grep PRETTY_NAME)

cd ${EDGE_TMP_DIR}

###############
# Basic Tools #
###############
if [[ ${EDGE_DIST} == *"Fedora"* ]]
then
  # recent build tools
  sudo dnf install -y -q -e 0 gcc-c++
  sudo dnf install -y -q -e 0 gcc-gfortran
  sudo dnf install -y -q -e 0 libasan libubsan
  sudo dnf install -y -q -e 0 glibc-static

  # other
  sudo dnf install -y -q -e 0 environment-modules
  sudo dnf install -y -q -e 0 hostname
  sudo dnf install -y -q -e 0 wget
  sudo dnf install -y -q -e 0 unzip
  sudo dnf install -y -q -e 0 bzip2
  sudo dnf install -y -q -e 0 expect
  sudo dnf install -y -q -e 0 m4
  sudo dnf install -y -q -e 0 autoconf
  sudo dnf install -y -q -e 0 automake
  sudo dnf install -y -q -e 0 libtool
  sudo dnf install -y -q -e 0 make
  sudo dnf install -y -q -e 0 cmake
  sudo dnf install -y -q -e 0 git
  sudo dnf install -y -q -e 0 git-lfs
  sudo dnf install -y -q -e 0 patch
  sudo dnf install -y -q -e 0 python3-libxml2 python3-scons
  sudo dnf install -y -q -e 0 python2 python3 python2-devel python3-devel python2-setuptools python3-setuptools python3-pip
  sudo dnf install -y -q -e 0 cppcheck
  sudo dnf install -y -q -e 0 valgrind
  sudo dnf install -y -q -e 0 irqbalance
  sudo dnf install -y -q -e 0 mesa-libGLU libXcursor libXft libXinerama
  sudo dnf install -y -q -e 0 atop
  sudo dnf install -y -q -e 0 htop
  sudo dnf install -y -q -e 0 npm
  sudo dnf install -y -q -e 0 yum-plugin-copr
  sudo dnf copr enable -y -q -e 0 genericmappingtools/gmt
  sudo dnf install -y -q -e 0 gmt
  sudo dnf install -y -q -e 0 ghostscript
  # UCVM and EDGEv dependencies
  sudo dnf install -y -q -e 0 proj-static proj-devel
  # EDGEcut dependencies
  sudo pip install meshio > /dev/null
  # Pipeline dependencies
  sudo dnf install -y -q -e 0 npm
elif [[ ${EDGE_DIST} == *"Ubuntu"* ]]
then
  export DEBIAN_FRONTEND=noninteractive
  sudo apt-get -qq update
  sudo apt install -qq -y build-essential
  sudo apt install -qq -y gfortran
  sudo apt install -qq -y autoconf
  sudo apt install -qq -y libtool
  sudo apt install -qq -y git
  sudo apt install -qq -y cmake
  sudo apt install -qq -y scons
  sudo apt install -qq -y python3 python3-pip python3-dev
  sudo apt install -qq -y cppcheck
  sudo apt install -qq -y valgrind
  sudo ln -s /usr/bin/python3 /usr/bin/python
  sudo ln -s /usr/bin/pip3 /usr/bin/pip
  # Gmsh dependencies
  sudo apt install libblas-dev liblapack-dev
fi

#########
# Clang #
#########
if [[ ${EDGE_DIST} == *"Fedora"* ]]
then
  sudo dnf install -y -q -e 0 llvm-devel clang-devel
  sudo dnf install -y -q -e 0 libomp-devel
fi

###########
# OpenMPI #
###########
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz -O openmpi.tar.gz
mkdir openmpi
tar -xf openmpi.tar.gz -C openmpi --strip-components=1
cd openmpi
./configure --enable-static=yes --enable-mpi1-compatibility > /dev/null
make -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo make install > /dev/null
sudo ldconfig > /dev/null
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

#########
# ArmIE #
#########
if [ -f ${EDGE_CURRENT_DIR}/ARM-Instruction-Emulator_*.tar.gz ]
then
  tar -xzf ${EDGE_CURRENT_DIR}/ARM-Instruction-Emulator_*.tar.gz
  cd ARM-Instruction-Emulator*
  sudo ./arm-instruction-emulator-* -a
  udo bash -c 'echo export MODULEPATH=\$MODULEPATH:/opt/arm/modulefiles/ >> /etc/bashrc'
fi

##################
# Python modules #
##################
sudo pip -q install xmltodict matplotlib obspy

#############
# Buildkite #
#############
if [[ ${EDGE_DIST} == *"Fedora"* ]]
then
  if [[ $(uname -m) == "aarch64" ]]
  then
    sudo sh -c 'echo -e "[buildkite-agent]\nname = Buildkite Pty Ltd\nbaseurl = https://yum.buildkite.com/buildkite-agent/stable/aarch64/\nenabled=1\ngpgcheck=0\npriority=1" > /etc/yum.repos.d/buildkite-agent.repo'
  else
    sudo sh -c 'echo -e "[buildkite-agent]\nname = Buildkite Pty Ltd\nbaseurl = https://yum.buildkite.com/buildkite-agent/stable/x86_64/\nenabled=1\ngpgcheck=0\npriority=1" > /etc/yum.repos.d/buildkite-agent.repo'
  fi
  sudo dnf install -y -q -e 0 buildkite-agent
fi

############
# Clean up #
############
if [[ ${EDGE_DIST} == *"Fedora"* ]]
then
  sudo dnf clean all
  sudo rm -rf /var/cache/dnf
fi

sudo rm -rf ${EDGE_TMP_DIR}
cd ${EDGE_CURRENT_DIR}
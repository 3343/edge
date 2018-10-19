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
# Installs the tools used in CI/CD (tested for Debian 9.3, Ubuntu 16.04).
##
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
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y gcc
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y gfortran
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y libiomp5 libiomp-dev
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y software-properties-common

# install clang
sudo echo "deb http://apt.llvm.org/stretch/ llvm-toolchain-stretch-7 main"     >  /etc/apt/sources.list.d/clang.list
sudo echo "deb-src http://apt.llvm.org/stretch/ llvm-toolchain-stretch-7 main" >> /etc/apt/sources.list.d/clang.list
sudo wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key|sudo apt-key add -
sudo apt-get update
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y clang-7 lldb-7 lld-7
sudo ln -s /usr/bin/clang-7   /usr/bin/clang
sudo ln -s /usr/bin/clang++-7 /usr/bin/clang++

# install valgrind
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 libc6-dbg
wget http://mirrors.kernel.org/sourceware/valgrind/valgrind-3.13.0.tar.bz2
tar -xjf valgrind-3.13.0.tar.bz2
cd valgrind-3.13.0
./configure
make
sudo make install
cd ..

# install Git LFS
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 -y git-lfs

# install python modules
sudo pip -q install xmltodict matplotlib netCDF4
sudo pip3 -q install xmltodict matplotlib netCDF4

# install Vagrant
wget https://releases.hashicorp.com/vagrant/2.1.2/vagrant_2.1.2_x86_64.deb -O vagrant.deb
sudo dpkg -i vagrant.deb

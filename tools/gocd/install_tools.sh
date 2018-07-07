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
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 wget
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 unzip
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 m4
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 dh-autoreconf
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 make
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 cmake
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 scons
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 git
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 libxml2-utils
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 python-pip python3-pip
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 openmpi-bin
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 cppcheck
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 gmsh
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 gcc
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 gfortran
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 clang
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 valgrind
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
sudo apt-get install -qq -o=Dpkg::Use-Pty=0 git-lfs
sudo su go -c "git lfs install"
sudo pip -q install xmltodict matplotlib netCDF4
sudo pip3 -q install xmltodict matplotlib netCDF4
wget https://releases.hashicorp.com/vagrant/2.1.2/vagrant_2.1.2_x86_64.deb -O vagrant.deb
sudo dpkg -i vagrant.deb
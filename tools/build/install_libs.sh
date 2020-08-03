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
# Installs EDGE's dependencies system-wide (tested for Debian 9.3, CentOS 7).
##
EDGE_CURRENT_DIR=$(pwd)
EDGE_TMP_DIR=$(mktemp -d)
EDGE_N_BUILD_PROC=$(cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}')

# detect compilers
[[ $(type -P mpiicc)  ]] && export CC=mpiicc   || export CC=mpicc
[[ $(type -P mpiicpc) ]] && export CXX=mpiicpc || export CXX=mpiCC


cd ${EDGE_TMP_DIR}

########
# zlib #
########
wget http://zlib.net/zlib-1.2.11.tar.gz -O zlib.tar.gz > /dev/null
mkdir zlib
tar -xzf zlib.tar.gz -C zlib --strip-components=1
cd zlib
CFLAGS="-fPIC" ./configure --static --prefix=/usr/local > /dev/null
sudo make install -j ${EDGE_N_BUILD_PROC} > /dev/null
cd ..

########
# HDF5 #
########
wget https://www.hdfgroup.org/package/gzip/?wpdmdl=13048 -O hdf5.tar.gz > /dev/null
mkdir hdf5
tar -xzf hdf5.tar.gz -C hdf5 --strip-components=1
cd hdf5
./configure --enable-shared=no --prefix=/usr/local > /dev/null
make -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo make install > /dev/null
cd ..

####################
# LIBXSMM and MOAB #
####################
git clone https://github.com/3343/edge.git
cd edge
git checkout develop
git submodule init
git submodule update

# build libxsmm
cd submodules/libxsmm
sudo make FORTRAN=0 BLAS=0 PREFIX=/usr/local install -j ${EDGE_N_BUILD_PROC} > /dev/null
cd ../..

# build moab
cd submodules/moab
LANG=C autoreconf -fi
CXXFLAGS="-DEIGEN_DONT_VECTORIZE -fPIC" ./configure --disable-debug --disable-optimize --enable-shared=no --enable-static=yes --with-pic=yes --disable-fortran --enable-tools --disable-blaslapack --with-eigen3=$(pwd)/../eigen --with-hdf5=yes --with-netcdf=no --with-pnetcdf=no > /dev/null
make -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo make install > /dev/null
cd ../..

cd ..

########
# CGAL #
########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz -O cgal.tar.xz > /dev/null
  mkdir cgal
  tar -xf cgal.tar.xz -C cgal --strip-components=1
  cd cgal
  cmake . -DWITH_CGAL_ImageIO=OFF -DWITH_CGAL_Qt5=OFF -DCMAKE_BUILD_TYPE=Release > /dev/null
  make -j ${EDGE_N_BUILD_PROC} > /dev/null
  sudo make install > /dev/null
  cd ..
fi

###########
# EDGEcut #
###########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  cd edge/tools/edge_cut
  scons cgal_dir=/usr/local -j ${EDGE_N_BUILD_PROC}
  sudo mv ./build/edge_cut /usr/local/bin
  cd ../../..
fi

########
# UCVM #
########
git clone https://github.com/SCECcode/UCVMC.git
cd UCVMC
git checkout v19.4.2
cd largefiles
# CVM-S4.26 and CVM-H
bash -c "printf \"no\nno\nno\nno\nno\nyes\nno\nyes\n\" | ./get_large_files.py"
./check_largefiles_md5.py
./stage_large_files.py
cd ..
# CVM-S4.26 and CVM-H
sudo bash -c "printf \"/usr/local/ucvm\nyes\nyes\n\" | PATH=$PATH python2 ./ucvm_setup.py" > /dev/null
cd ..

##########
# EDGE-V #
##########
cd edge/tools/edge_v
scons moab=yes ucvm=/usr/local/ucvm zlib=yes hdf5=yes netcdf=no
sudo mv build/edge_v /usr/local/bin
cd ../..

############
# Clean up #
############
sudo rm -rf ${EDGE_TMP_DIR}

cd ${EDGE_CURRENT_DIR}
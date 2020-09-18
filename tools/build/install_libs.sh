#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2020, Friedrich Schiller University Jena
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

git clone https://github.com/3343/edge.git > /dev/null
cd edge
git checkout develop > /dev/null
git submodule init > /dev/null
git submodule update > /dev/null

#######################
# Run default scripts #
#######################
sudo ./tools/build/deps/zlib.sh -o /usr/local -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo ./tools/build/deps/hdf5.sh -z /usr/local -o /usr/local -j ${EDGE_N_BUILD_PROC} > /dev/null
sudo ./tools/build/deps/libxsmm.sh -i $(pwd)/submodules/libxsmm -o /usr/local -j ${EDGE_N_BUILD_PROC} > /dev/null

###########
# EDGEcut #
###########
if [[ ${EDGE_DIST} == *"CentOS"* ]]
then
  cd tools/edge_cut
  scons -j ${EDGE_N_BUILD_PROC}
  sudo mv ./build/edge_cut /usr/local/bin
  cd ../..
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
cd tools/edge_v
sudo bash -c "PATH=$PATH scons parallel=omp ucvm=/usr/local/ucvm install_dir=/usr/local -j ${EDGE_N_BUILD_PROC}" > /dev/null
cd ../..

############
# Clean up #
############
sudo rm -rf ${EDGE_TMP_DIR}

cd ${EDGE_CURRENT_DIR}
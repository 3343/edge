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
# Installation of HPC-related tools.
##
EDGE_CURRENT_DIR=$(pwd)
EDGE_TMP_DIR=$(mktemp -d)
EDGE_N_BUILD_PROC=$(cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}')

EDGE_DIST=$(cat /etc/*-release | grep PRETTY_NAME)

# detect compilers
[[ $(type -P mpiicc)  ]] && export CC=mpiicc   || export CC=mpicc
[[ $(type -P mpiicpc) ]] && export CXX=mpiicpc || export CXX=mpiCC
[[ $(type -P mpiifort) ]] && export F90=mpiifort || export CXX=mpifort

cd ${EDGE_TMP_DIR}

##########
# STREAM #
##########
wget http://www.cs.virginia.edu/stream/FTP/Code/stream.c
if [[ $(type -P icc) ]]
then
  icc -Ofast -xhost -qopenmp \
      -DNTIMES=1000 -DSTREAM_ARRAY_SIZE=$(cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}')000000 \
      -qopt-streaming-cache-evict=0 -qopt-streaming-stores always -qopt-prefetch-distance=64,8 stream.c -o stream-bench
else
  gcc -Ofast -march=native -fopenmp \
      -DNTIMES=1000 -DSTREAM_ARRAY_SIZE=$(cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}')000000 \
      stream.c -o stream-bench
fi
sudo mv ./stream-bench /usr/local/bin

###########
# OSU MPI #
###########
wget http://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-5.4.4.tar.gz -O osu.tar.gz
mkdir osu; tar -xzf osu.tar.gz -C osu --strip-components=1
cd osu
./configure
make -j ${EDGE_N_BUILD_PROC}
sudo make install
for osu_type in collective one-sided pt2pt startup
do
  echo "export PATH=/usr/local/libexec/osu-micro-benchmarks/mpi/${osu_type}:\${PATH}" | sudo tee --append /etc/bashrc
done
cd ..

###########
# SCORE-P #
###########
wget https://www.vi-hps.org/upload/packages/scorep/scorep-4.1.tar.gz -O scorep.tar.gz
mkdir scorep; tar -xzf scorep.tar.gz -C scorep --strip-components=1
cd scorep
if [[ $(type -P mpiicc) ]] && [[ $(type -P mpiicpc) ]] && [[ $(type -P mpiifort) ]]
then
  ./configure --disable-static --enable-shared --with-mpi=intel --with-nocross-compiler-suite=intel
else
  ./configure --disable-static --enable-shared
fi
make -j ${EDGE_N_BUILD_PROC}
sudo make install
echo "export PATH=/opt/scorep/bin:\${PATH}" | sudo tee --append /etc/bashrc
cd ..

############
# SCALASCA #
############
wget http://apps.fz-juelich.de/scalasca/releases/scalasca/2.4/dist/scalasca-2.4.tar.gz -O scalasca.tar.gz
mkdir scalasca; tar -xzf scalasca.tar.gz -C scalasca --strip-components=1
cd scalasca
if [[ $(type -P mpiicc) ]] && [[ $(type -P mpiicpc) ]] && [[ $(type -P mpiifort) ]]
then
  ./configure --disable-static --enable-shared --with-mpi=intel --with-nocross-compiler-suite=intel
else
  ./configure --disable-static --enable-shared
fi
make -j ${EDGE_N_BUILD_PROC}
sudo make install
echo "export PATH=/opt/scalasca/bin:\${PATH}" | sudo tee --append /etc/bashrc
cd ..

############
# Clean up #
############
sudo rm -rf ${EDGE_TMP_DIR}

cd ${EDGE_CURRENT_DIR}
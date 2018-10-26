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
./configure
make -j ${EDGE_N_BUILD_PROC}
sudo make install\
for osu_type in collective one-sided pt2pt startup
do
  echo "export PATH=/usr/local/libexec/osu-micro-benchmarks/mpi/${osu_type}:\${PATH}" | sudo tee --append /etc/bashrc
done

############
# Clean up #
############
sudo rm -rf ${EDGE_TMP_DIR}

cd ${EDGE_CURRENT_DIR}
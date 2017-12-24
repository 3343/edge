#!/bin/sh
#
# Copyright (c) 2017-2018, Intel Corporation
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Intel Corporation nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

# global build group variables
EDGE_CXX=mpiicpc
EDGE_PAR_COMPILE=16
EDGE_ARCH=avx512
EDGE_PARALLEL=mpi+omp
EDGE_ELEMENT=tet4
EDGE_EQUATION=elastic
EDGE_BIN_DIR=./bin/

# configs to build, $order_$precision_$cfr
EDGE_CONFIGS="2_64_1 2_32_1 2_64_8 2_32_16 3_64_1 3_32_1 3_64_8 3_32_16 4_64_1 4_32_1 4_64_8 4_32_16 5_64_1 5_32_1 5_64_8 5_32_16 6_64_1 6_32_1 6_64_8 6_32_16 7_64_1 7_32_1 7_64_8 7_32_16"

# save current location
PWD_JUMP_BACK=`pwd`

if [ -z ${EDGE_ROOT+x} ]
then 
  echo "EDGE_ROOT is not set, please set it!"
  exit -1
fi

# switch into EDGE root
cd ${EDGE_ROOT}

# remove bin dir
mkdir -p ${EDGE_BIN_DIR}

for c in ${EDGE_CONFIGS}
do
  # extract config
  EDGE_ORDER=`echo ${c} | awk -F"_" '{print $1}'`
  EDGE_PRECISION=`echo ${c} | awk -F"_" '{print $2}'`
  EDGE_CFR=`echo ${c} | awk -F"_" '{print $3}'`

  # cleanup
  rm -rf build/
  rm -rf .sconf_temp
  rm -rf .sconsign.dblite

  # finally let's build the code
  CXX=${EDGE_CXX} scons equations=${EDGE_EQUATION} order=${EDGE_ORDER} precision=${EDGE_PRECISION} cfr=${EDGE_CFR} element_type=${EDGE_ELEMENT} parallel=${EDGE_PARALLEL} arch=${EDGE_ARCH} xsmm=./libs zlib=./libs hdf5=./libs netcdf=./libs moab=./libs -j ${EDGE_PAR_COMPILE}

  # copy binary
  cp ./build/edge ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
  unlink ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
  ln -s edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR} ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
done

cd ${PWD_JUMP_BACK}

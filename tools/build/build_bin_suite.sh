#!/bin/sh
#
# Copyright (c) 2017-2018, Intel Corporation
# Copyright (c) 2017-2018, Regents of the University of California
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
if [[ -z $EDGE_DEPS ]]
then
  EDGE_DEPS=./deps
fi

if [[ -z $EDGE_CXX  ]]
then
  EDGE_CXX=mpiCC
fi

if [[ -z $EDGE_PAR_COMPILE ]]
then
  EDGE_PAR_COMPILE=16
fi

if [[ -z $EDGE_ARCH ]]
then
  EDGE_ARCH=avx512
fi

if [[ -z $EDGE_PARALLEL ]]
then
  EDGE_PARALLEL=mpi+omp
fi

if [[ -z $EDGE_ELEMENT ]]
then
  EDGE_ELEMENT=tet4
fi

if [[ -z $EDGE_EQUATION ]]
then
  EDGE_EQUATION=elastic
fi

if [[ -z $EDGE_BIN_DIR ]]
then
  EDGE_BIN_DIR=./bin/
fi

# configs to build, $order_$precision_$cfr
if [[ -z $EDGE_CONFIGS ]]
then
  if [[ "${EDGE_ARCH}" == "avx512" || "${EDGE_ARCH}" == "skx" || "${EDGE_ARCH}" == "knl" ]]
  then
    EDGE_CONFIGS="2_64_1 2_32_1 2_64_8 2_32_16 3_64_1 3_32_1 3_64_8 3_32_16 4_64_1 4_32_1 4_64_8 4_32_16 5_64_1 5_32_1 5_64_8 5_32_16 6_64_1 6_32_1 6_64_8 6_32_16 7_64_1 7_32_1 7_64_8 7_32_16"
  elif [[ "${EDGE_ARCH}" == "snb" || "${EDGE_ARCH}" == "hsw" ]]
  then
    EDGE_CONFIGS="2_64_1 2_32_1 2_64_4 2_32_8 3_64_1 3_32_1 3_64_4 3_32_8 4_64_1 4_32_1 4_64_4 4_32_8 5_64_1 5_32_1 5_64_4 5_32_8 6_64_1 6_32_1 6_64_4 6_32_8 7_64_1 7_32_1 7_64_4 7_32_8"
  else
    echo "unknown arch was specified!"
    exit -1
  fi
fi

# save current location
PWD_JUMP_BACK=`pwd`

if [[ -z ${EDGE_ROOT+x} ]]
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
  CXX=${EDGE_CXX} scons equations=${EDGE_EQUATION} order=${EDGE_ORDER} precision=${EDGE_PRECISION} cfr=${EDGE_CFR} element_type=${EDGE_ELEMENT} parallel=${EDGE_PARALLEL} arch=${EDGE_ARCH} xsmm=${EDGE_DEPS} zlib=${EDGE_DEPS} hdf5=${EDGE_DEPS} netcdf=${EDGE_DEPS} moab=${EDGE_DEPS} -j ${EDGE_PAR_COMPILE}

  # copy binary
  cp ./build/edge ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
  if [[ -f ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR} ]]
   then
    unlink ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
  fi
  ln -s edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR} ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
done

cd ${PWD_JUMP_BACK}

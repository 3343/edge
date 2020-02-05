#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (breuer AT mytum.de)
#
# @section LICENSE
# Copyright (c) 2019, Alexander Breuer
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
# Installs MOAB.
##
help() {
cat << EOF
Usage: ${0##*/} [-h -d] [-z ZLIB_INSTALL_DIR -5 HDF5_INSTALL_DIR -e EIGEN_DIR -i MOAB_DIR -o INSTALL_DIR -j N_BUILD_PROCS]
Installs MOAB.
     -h This help message.
     -d Enables debug build if set.
     -z ZLIB_INSTALL_DIR Absolute path of the zlib installation directory.
     -5 HDF5_INSTALL_DIR Absolute path of the hdf5 installation directory.
     -e EIGEN_DIR Absolute path of the Eigen directory.
     -i MOAB_DIR Absolute path of the MOAB source directory.
     -o INSTALL_DIR Absolute path of the installation directory, will be created if missing.
     -j N_BUILD_PROCS (optional) number of build processes, defaults to one.
EOF
}

DEBUG_BUILD=0
while getopts "hz:5:e:i:o:j:d" opt; do
  case "$opt" in
    h)
      help
      exit 0
      ;;
    z)
      ZLIB_INSTALL_DIR=$OPTARG
      ;;
    5)
      HDF5_INSTALL_DIR=$OPTARG
      ;;
    e)
      EIGEN_DIR=$OPTARG
      ;;
    i)
      MOAB_DIR=$OPTARG
      ;;
    o)
      INSTALL_DIR=$OPTARG
      ;;
    d)
      DEBUG_BUILD=1
      ;;
    j)
      N_BUILD_PROCS=$OPTARG
      ;;
    '?')
       help >&2
       exit 1
       ;;
    esac
done
shift "$((OPTIND-1))"

if [[ ${ZLIB_INSTALL_DIR:0:1} != "/" ]]
then
  echo "Error: ${ZLIB_INSTALL_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${HDF5_INSTALL_DIR:0:1} != "/" ]]
then
  echo "Error: ${HDF5_INSTALL_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${MOAB_DIR:0:1} != "/" ]]
then
  echo "Error: ${MOAB_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${EIGEN_DIR:0:1} != "/" ]]
then
  echo "Error: ${EIGEN_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${INSTALL_DIR:0:1} != "/" ]]
then
  echo "Error: ${INSTALL_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${N_BUILD_PROCS} == "" ]]
then
  N_BUILD_PROCS=1
fi

# move to temporary directory
TMP_DIR=$(mktemp -d)
cd ${TMP_DIR}

# copy over and move into dir
cp -r ${EIGEN_DIR} .
cp -r ${MOAB_DIR} .
cd $(basename ${MOAB_DIR})

autoreconf -fi

BUILD_FLAGS="--enable-shared=no              \
             --enable-static=yes             \
             --with-pic=yes                  \
             --disable-fortran               \
             --enable-tools                  \
             --with-zlib=${ZLIB_INSTALL_DIR} \
             --with-hdf5=${HDF5_INSTALL_DIR} \
             --with-netcdf=no                \
             --disable-blaslapack            \
             --with-eigen3=${EIGEN_DIR}      \
             --with-pnetcdf=no               \
             --with-metis=no                 \
            --prefix=${INSTALL_DIR}"

if [[ ${DEBUG_BUILD} == 0 ]]
then
  BUILD_FLAGS="--disable-debug   \
               --enable-optimize \
               ${BUILD_FLAGS}"
else
  BUILD_FLAGS="--enable-debug   \
               --disable-optimize \
               ${BUILD_FLAGS}"
fi

CXXFLAGS="-fPIC" LIBS="${ZLIB_INSTALL_DIR}/lib/libz.a" ./configure ${BUILD_FLAGS}
make -j ${N_BUILD_PROCS}
make install
rm ${INSTALL_DIR}/lib/*MOAB*.la
rm ${INSTALL_DIR}/lib/*iMesh*.la

rm -r ${TMP_DIR}

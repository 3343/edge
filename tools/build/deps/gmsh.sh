#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (breuer AT mytum.de)
#
# @section LICENSE
# Copyright (c) 2020, Alexander Breuer
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
# Installs Gmsh.
##
GMSH_LINK=https://gmsh.info/src/gmsh-4.6.0-source.tgz
GMSH_SHA256=0f2c55e50fb6c478ebc8977f6341c223754cbf3493b7b0d683b4395ae9f2ad1c

help() {
cat << EOF
Usage: ${0##*/} [-h] [-l BLAS_LAPACK_LIBRARIES -p EDGE_PARTITION_PLUGIN_DIR -o INSTALL_DIR -j N_BUILD_PROCS]
Installs Gmsh.
     -h This help message.
     -l BLAS_LAPACK_LIBRARIES (optional) string passed to Gmsh for the BLAS/LAPACK installation.
     -p EDGE_PARTITION_PLUGIN_DIR path to EDGE's partition plugin.
     -o INSTALL_DIR Absolute path of the installation directory, will be created if missing.
     -j N_BUILD_PROCS (optional) number of build processes, defaults to one.
EOF
}

while getopts "hl:p:o:j:" opt; do
  case "$opt" in
    h)
      help
      exit 0
      ;;
    l)
      BLAS_LAPACK_LIBRARIES=$OPTARG
      ;;
    p)
      EDGE_PARTITION_PLUGIN_DIR=$OPTARG
      ;;
    o)
      INSTALL_DIR=$OPTARG
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

if [[ ${INSTALL_DIR:0:1} != "/" ]]
then
  echo "Error: ${INSTALL_DIR} is not absolute"
  help >&2
  exit 1
fi

if [[ ${EDGE_PARTITION_PLUGIN_DIR} != "" ]]
then
  if [[ ${EDGE_PARTITION_PLUGIN_DIR:0:1} != "/" ]]
  then
    echo "Error: ${EDGE_PARTITION_PLUGIN_DIR} is not absolute"
    help >&2
    exit 1
  fi
fi

if [[ ${N_BUILD_PROCS} == "" ]]
then
  N_BUILD_PROCS=1
fi

# move to temporary directory
TMP_DIR=$(mktemp -d)
cd ${TMP_DIR}

# download Gmsh
wget ${GMSH_LINK} -O gmsh.tgz
if [[ $(sha256sum gmsh.tgz | cut -d' ' -f1) != ${GMSH_SHA256} ]]
then
  echo "Error: Checksum not matching"
  exit 1
fi

mkdir -p gmsh/build
tar -xf gmsh.tgz -C gmsh --strip-components=1
cd gmsh

# patch Gmsh with partition plugin if requested
if [[ ${EDGE_PARTITION_PLUGIN_DIR} != "" ]]
then
  patch -p1 < ${EDGE_PARTITION_PLUGIN_DIR}/gmsh.patch
  cp ${EDGE_PARTITION_PLUGIN_DIR}/EdgePartition.* Plugin
fi

cd build

# assemble cmake command
cmake_command="cmake"
if [[ ${BLAS_LAPACK_LIBRARIES} != "" ]]
then
  cmake_command="${cmake_command} -DBLAS_LAPACK_LIBRARIES=${BLAS_LAPACK_LIBRARIES}"
fi
cmake_command="${cmake_command} -DENABLE_BUILD_LIB=ON -DENABLE_OPENMP=OFF -DENABLE_MMG3D=OFF -DENABLE_FLTK=OFF -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} .."

# configure and build
$(${cmake_command})
make -j ${N_BUILD_PROCS}
make install

rm -r ${TMP_DIR}

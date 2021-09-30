#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (breuer AT mytum.de)
#
# @section LICENSE
# Copyright (c) 2021, Friedrich Schiller University Jena
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
# Installs HDF5.
##
set -euo pipefail

HDF5_LINK=https://s3.amazonaws.com/hdf-wordpress-1/wp-content/uploads/manual/HDF5/HDF5_1_10_5/source/hdf5-1.10.5.tar.bz2
HDF5_SHA256=68d6ea8843d2a106ec6a7828564c1689c7a85714a35d8efafa2fee20ca366f44

help() {
cat << EOF
Usage: ${0##*/} [-h] [-z ZLIB_INSTALL_DIR -o INSTALL_DIR -j N_BUILD_PROCS]
Installs HDF5.
     -h This help message.
     -z ZLIB_INSTALL_DIR path to the zlib installation directory.
     -o INSTALL_DIR Absolute path of the installation directory, will be created if missing.
     -j N_BUILD_PROCS (optional) number of build processes, defaults to one.
EOF
}

while getopts "hz:o:j:" opt; do
  case "$opt" in
    h)
      help
      exit 0
      ;;
    z)
      ZLIB_INSTALL_DIR=$OPTARG
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

if [[ ${ZLIB_INSTALL_DIR:0:1} != "/" ]]
then
  echo "Error: ${ZLIB_INSTALL_DIR} is not absolute"
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

# download hdf5
wget --no-check-certificate ${HDF5_LINK} -O hdf5.tar.bz2
if [[ $(sha256sum hdf5.tar.bz2 | cut -d' ' -f1) != ${HDF5_SHA256} ]]
then
  echo "Error: Checksum not matching"
  exit 1
fi

mkdir hdf5
tar -xf hdf5.tar.bz2 -C hdf5 --strip-components=1
cd hdf5

LDFLAGS="-fPIC" CFLAGS="-fPIC" ./configure --enable-shared=no --with-zlib=${ZLIB_INSTALL_DIR} --prefix=${INSTALL_DIR}
make -j ${N_BUILD_PROCS}
make install
rm ${INSTALL_DIR}/lib/*hdf5*.la

rm -r ${TMP_DIR}
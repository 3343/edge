#!/bin/bash
set -euo pipefail
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
# Installs zlib.
##
set -euo pipefail

ZLIB_LINK=https://www.zlib.net/zlib-1.2.11.tar.xz
ZLIB_SHA256=4ff941449631ace0d4d203e3483be9dbc9da454084111f97ea0a2114e19bf066

help() {
cat << EOF
Usage: ${0##*/} [-h] [-o INSTALL_DIR -j N_BUILD_PROCS]
Installs zlib.
     -h This help message.
     -o INSTALL_DIR Absolute path of the installation directory, will be created if missing.
     -j N_BUILD_PROCS (optional) number of build processes, defaults to one.
EOF
}

while getopts "ho:j:" opt; do
  case "$opt" in
    h)
      help
      exit 0
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

if [[ ${N_BUILD_PROCS} == "" ]]
then
  N_BUILD_PROCS=1
fi

# move to temporary directory
TMP_DIR=$(mktemp -d)
cd ${TMP_DIR}

# download zlib
wget --no-check-certificate ${ZLIB_LINK} -O zlib.tar.xz
if [[ $(sha256sum zlib.tar.xz | cut -d' ' -f1) != ${ZLIB_SHA256} ]]
then
  echo "Error: Checksum not matching"
  exit 1
fi

mkdir zlib
tar -xf zlib.tar.xz -C zlib --strip-components=1
cd zlib
LDFLAGS="-fPIC" CFLAGS="-fPIC" ./configure --static --prefix=${INSTALL_DIR}
make -j ${N_BUILD_PROCS}
make install
rm ${INSTALL_DIR}/lib/pkgconfig/zlib.pc

rm -r ${TMP_DIR}
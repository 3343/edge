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
# Installs METIS.
##
set -euo pipefail

METIS_LINK=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
METIS_SHA256=76faebe03f6c963127dbb73c13eab58c9a3faeae48779f049066a21c087c5db2

help() {
cat << EOF
Usage: ${0##*/} [-h] [-p EDGE_PATCH_PATH -o INSTALL_DIR -j N_BUILD_PROCS]
Installs METIS.
     -h This help message.
     -p EDGE_PATCH_PATH patch-file specific to EDGE.
     -o INSTALL_DIR Absolute path of the installation directory, will be created if missing.
     -j N_BUILD_PROCS (optional) number of build processes, defaults to one.
EOF
}

while getopts "hp:o:j:" opt; do
  case "$opt" in
    h)
      help
      exit 0
      ;;
    p)
      EDGE_PATCH_PATH=$OPTARG
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

if [[ ${EDGE_PATCH_PATH:0:1} != "/" ]]
then
  echo "Error: ${EDGE_PATCH_PATH} is not absolute"
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

# download METIS
wget --no-check-certificate ${METIS_LINK} -O metis.tar.gz
if [[ $(sha256sum metis.tar.gz | cut -d' ' -f1) != ${METIS_SHA256} ]]
then
  echo "Error: Checksum not matching"
  exit 1
fi

mkdir metis
tar -xf metis.tar.gz -C metis --strip-components=1
cd metis

# patch if requested
if [[ ${EDGE_PATCH_PATH} != "" ]]
then
  patch -p1 < ${EDGE_PATCH_PATH}
fi

# install
make config prefix=${INSTALL_DIR}
make -j ${N_BUILD_PROCS}
make install

rm -r ${TMP_DIR}
#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (alex.breuer AT uni-jena.de)
#
# @section LICENSE
# Copyright (c) 2020, Friedrich Schiller University Jena
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
# Creates a visit file for the given EDGE-output directory.
help() {
cat << EOF
Usage: ${0##*/} [-h] [-i INPUT_DIR -t TIME_STEP -o INSTALL_DIR]
Creates a visit file for the given EDGE-output directory.
     -h This help message.
     -i INPUT_DIR directory containing EDGE's wave field output.
     -t INITIAL_TIME initial time of the first output (default: 0).
     -s TIME_STEP optional time step (default: 1).
     -o OUTPUT_PATH output path to the visit-file.
EOF
}

TIME_STEP=1
INITIAL_TIME=0

while getopts "hi:t:s:o:" opt; do
  case "$opt" in
    h)
      help
      exit 0
      ;;
    i)
      INPUT_DIR=$OPTARG
      ;;
    t)
      INITIAL_TIME=$OPTARG
      ;;
    s)
      TIME_STEP=$OPTARG
      ;;
    o)
      OUTPUT_PATH=$OPTARG
      ;;
    '?')
       help >&2
       exit 1
       ;;
    esac
done
shift "$((OPTIND-1))"

if [[ ${INPUT_DIR} == "" ]]
then
  echo "Error input dir missing"
  help >&2
  exit 1
fi

if [[ ${OUTPUT_PATH} == "" ]]
then
  echo "Error output path missing"
  help >&2
  exit 1
fi

n_blocks=$(ls -d ${INPUT_DIR}/* | wc -l)
n_steps=$(ls ${INPUT_DIR/0/*} | wc -l)

echo "!NBLOCKS ${n_blocks}" > ${OUTPUT_PATH}
for time in $(seq 0 $(expr ${n_steps} - 1) )
do
  abs_time=$(echo "${INITIAL_TIME}+${time}*${TIME_STEP}" | bc)
  echo "!TIME ${abs_time}" >> ${OUTPUT_PATH}
  for part in $(ls -d ${INPUT_DIR}/*)
  do
    echo ${part}/wf_$(basename ${part})_${time}.vtk >> ${OUTPUT_PATH}
  done
done
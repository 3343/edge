##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018-2019, Regents of the University of California
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
# Starts the given number of GoCD agents.
##

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-n N_AGENTS -p ABS_PATH_TO_DIR]
Starts the given number of agents in the given directory.
     -h This help message.
     -n number of agents to start.
     -p absolute path to the directory, where the agents store there info (agent_i with i in 1 .. N_AGENTS).
EOF
}

# get arguments
while getopts "hn:p:" opt; do
  case "$opt" in
    h)
       show_help
       exit 0
       ;;
    n)
       N_AGENTS=$OPTARG
       ;;
    p)
       ABS_PATH_TO_DIR=$OPTARG
       ;;
    '?')
       show_help >&2
       exit 1
       ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [[ $OPTIND < 3 ]]
then
  show_help >&2
  exit 1
fi

echo "setting server certificate"
sudo mkdir -p /var/go
sudo bash -c "openssl s_client -showcerts -connect deb.dial3343.org:8154 </dev/null 2>/dev/null|openssl x509 -outform PEM > /var/go/root-cert.pem"

CUR_DIR=$(pwd)
TMP_DIR=$(mktemp -d)
cd ${TMP_DIR}

echo "downloading GoCD"
wget https://download.gocd.org/binaries/18.12.0-8222/generic/go-agent-18.12.0-8222.zip
unzip go-agent-18.12.0-8222.zip
ln -s go-agent-18.12.0 go-agent

cd ${CUR_DIR}
for agent in `seq 1 ${N_AGENTS}`
do
  echo "setting up dir for agent #${agent}"
  mkdir ${ABS_PATH_TO_DIR}/agent_${agent}

  echo "starting agent #${agent}"
  cd ${ABS_PATH_TO_DIR}/agent_${agent}
  AGENT_WORK_DIR=$(pwd) GO_SERVER_URL=https://deb.dial3343.org:8154/go nohup ${TMP_DIR}/go-agent/agent.sh &
done
cd ${CUR_DIR}

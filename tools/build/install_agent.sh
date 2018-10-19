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
# Installs the GoCD agent and adds the certificate of the given server.
##

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-s HOST:PORT]
Installs a GoCD server .
     -h This help message.
     -s HOST:PORT host and port of the GoCD server for which the certificate is added.
EOF
}

# get arguments
while getopts "hs:" opt; do
  case "$opt" in
    h)
       show_help
       exit 0
       ;;
    s)
       GOCD_SERVER=$OPTARG
       ;;
    '?')
       show_help >&2
       exit 1
       ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [[ $OPTIND < 2 ]]
then
  show_help >&2
  exit 1
fi

# add GoCD packages and key
sudo echo "deb https://download.gocd.org /" > /etc/apt/sources.list.d/gocd.list
curl -s https://download.gocd.org/GOCD-GPG-KEY.asc | sudo apt-key add -
sudo apt-get -qq update

# install Java 8
sudo add-apt-repository -y ppa:openjdk-r/ppa > /dev/null
sudo apt-get -qq update
sudo apt-get -qq install openjdk-8-jre > /dev/null

sudo apt-get -qq install go-agent
sudo su go -c "git lfs install"

# add server key
openssl s_client -showcerts -connect ${GOCD_SERVER} </dev/null 2>/dev/null|openssl x509 -outform PEM > /var/go/root-cert.pem

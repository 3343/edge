##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017, Regents of the University of California
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
# Required packages for continuous delivery through GoCD (http://linuxbrew.sh/).
# Tested on: Cent OS 7, Debian 8.7
#
# Install essentials first:
#   sudo apt-get install build-essential curl git python-setuptools ruby
#
# glibc-devel is most certainly required before installing anything through linuxbrew:
#   sudo apt-get install libc6-dbg  # Debian
#   sudo yum install glibc-devel
##

ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
PATH="$HOME/.linuxbrew/bin:$PATH"

# disable use of linux brew glibc (e.g., causes trouble with Intel compilers)
# consider using linuxbrew's bottle if your glibc >= 2.19 (see http://linuxbrew.sh/)
export HOMEBREW_BUILD_FROM_SOURCE=1
echo 'export HOMEBREW_BUILD_FROM_SOURCE=1' >>~/.bashrc

# disable auto-uipdate as it tends to break things
export HOMEBREW_NO_AUTO_UPDATE=1
echo 'export HOMEBREW_NO_AUTO_UPDATE=1' >>~/.bashrc


# disable analytics
brew analytics off

KNL=no
SYSTEM_GCC=no
LINUXBREW_PREFIX=~/.linuxbrew

if [ "SYSTEM_GCC" == "yes" ]
then
  GCC_VERSION="$(gcc -dumpversion | cut -f1-2 -d.)"
  ln -s $(which gcc)      $LINUXBREW_PREFIX/bin/gcc-$GCC_VERSION
  ln -s $(which g++)      $LINUXBREW_PREFIX/bin/g++-$GCC_VERSION
  ln -s $(which gfortran) $LINUXBREW_PREFIX/bin/gfortran-$GCC_VERSION
  export HOMEBREW_CXX=g++-$GCC_VERSION
  export HOMEBREW_CC=gcc-$GCC_VERSION
  brew tap homebrew/versions
else
  brew install gcc
fi

brew install binutils --without-gold # fails with gold linker

if [ "KNL" == "yes" ]
then
  # workaround for KNL
  # avoids "140250 Illegal instruction" in installation
  ln -s /usr/bin/gcc .linuxbrew/bin/gcc-4.8
  ln -s /usr/bin/g++ .linuxbrew/bin/g++-4.8
  HOMEBREW_CC=gcc-4.8 HOMEBREW_CXX=g++-4.8 brew install python
else
  brew install python
fi

pip install xmltodict
pip install matplotlib
#brew install makedepend # possible cmake dependency
brew install cmake
brew install cppcheck
brew install lapack # unresolved gmsh dependency
brew install homebrew/science/gmsh
brew install llvm --without-glibc
brew install scons
brew install valgrind
brew install r --without-x11

echo "Remark: Source the linuxbrew profile accordingly in /etc/default/go-agent*, e.g. through source /var/go/.bash_profile"

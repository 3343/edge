#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2016-2017, Regents of the University of California
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
# Replaces vars and functions with compatible notation.
##

#for l_file in SolverElastic2D.inc SolverElastic2DAn.inc
#for l_file in Trafo2D.inc TrafoInv2D.inc EntFixHarten2D.inc
for l_file in EntFixHarten3D.inc
#for l_file in Trafo3D.inc TrafoInv3D.inc FlMid3D.inc FrMid3D.inc
#for l_file in MiddleStateJumpL2D.inc MiddleStateJumpR2D.inc MiddleStateFlux2D.inc
#for l_file in MiddleStateJumpL3D.inc MiddleStateJumpR3D.inc MiddleStateFlux3D.inc
do
echo "processing ${l_file}"

# remove unwanted stuff
sed -i 's/"//g' ${l_file}
sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\n//g' ${l_file}

# replace math functions
sed -i 's/Power/std::pow/g' ${l_file}
sed -i 's/Sqrt/std::sqrt/g' ${l_file}

# replace variable names
if [[ $l_file == "MiddleStateFlux2D.inc" || $l_file == "MiddleStateFlux3D.inc" ]]
then
  sed -i 's/lam/i_lam/g' ${l_file}
else
  sed -i 's/laml/i_lamL/g' ${l_file}
  sed -i 's/lamr/i_lamR/g' ${l_file}
fi

if [[ $l_file == "MiddleStateFlux2D.inc" || $l_file == "MiddleStateFlux3D.inc" ]]
then
  sed -i 's/mu/i_mu/g' ${l_file}
else
  sed -i 's/mul/i_muL/g' ${l_file}
  sed -i 's/mur/i_muR/g' ${l_file}
fi

if [[ $l_file == "MiddleStateFlux2D.inc" || $l_file == "MiddleStateFlux3D.inc" ]]
then
  sed -i 's/rho/i_rho/g' ${l_file}
else
  sed -i 's/rhol/i_rhoL/g' ${l_file}
  sed -i 's/rhor/i_rhoR/g' ${l_file}
fi

sed -i 's/nx/i_nx/g' ${l_file}
sed -i 's/ny/i_ny/g' ${l_file}
sed -i 's/nz/i_nz/g' ${l_file}

sed -i 's/sx/i_sx/g' ${l_file}
sed -i 's/sy/i_sy/g' ${l_file}
sed -i 's/sz/i_sz/g' ${l_file}

sed -i 's/tx/i_tx/g' ${l_file}
sed -i 's/ty/i_ty/g' ${l_file}
sed -i 's/tz/i_tz/g' ${l_file}

sed -i 's/scale/i_scale/g' ${l_file}
sed -i 's/Max/std::max/g' ${l_file}

done

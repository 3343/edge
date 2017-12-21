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
# Converts array info to string.
##

##
# Converts 2d float input to a comma-separated (fastest dim) and newline-separated string.
#
# @param i_in 2d input.
# @return string representation.
##
def float2d( i_in ):
  # output string
  l_outStr = ''

  # iterate first dim
  for l_d1 in i_in:
    # iterate second dim
    for l_d2 in l_d1:
      l_outStr = l_outStr + str(float(l_d2)) + ','
    l_outStr = l_outStr + '\n'

  return l_outStr

##
# Converts 2d integal input to a comma-separated (fastest dim) and newline-separated string.
#
# @param i_in 2d input.
# @return string representation.
##
def int2d( i_in ):
  # output string
  l_outStr = ''

  # iterate over first dim
  for l_d1 in i_in:
    # iterate over second dim
    for l_d2 in l_d1:
      l_outStr = l_outStr + str(l_d2) + ','
    l_outStr = l_outStr + '\n'

  return l_outStr

##
# Converts 3d integral input to a comma-separated string (fastest dim), newline-separated (second fastest) and double-newline-separated string.
#
# @param i_in 3d input.
# @return string representation.
##
def int3d( i_in ):
  # output string
  l_outStr = ''

  # iterate over first dim
  for l_d1 in i_in:
    # iterate over second dim
    for l_d2 in l_d1:
      # iterate over third dim
      for l_d3 in l_d2:
        l_outStr = l_outStr + str(l_d3) + ','
      l_outStr = l_outStr + '\n'
    l_outStr = l_outStr + '\n'

  return l_outStr

##
# Converts 3d float input to a comma-separated string (fastest dim), newline-separated (second fastest) and double-newline-separated string. 
#
# @param i_in 3d input.
# @return string representation.
##
def float3d( i_in ):
  # output string
  l_outStr = ''

  # iterate over first dim
  for l_d1 in i_in:
    # iterate over second dim
    for l_d2 in l_d1:
      # iterate over third dim
      for l_d3 in l_d2:
        l_outStr = l_outStr + str(float(l_d3)) + ','
      l_outStr = l_outStr + '\n'
    l_outStr = l_outStr + '\n'

  return l_outStr

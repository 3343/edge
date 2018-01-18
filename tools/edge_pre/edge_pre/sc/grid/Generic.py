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
# Generic derivation of sub-grid information.
##

##
# Derives sub-cells adjacent to sub-cells (faces as bridge).
#
# @param i_scSfSvIn sub-vertices adjacent to inner sub-cells (faces as bridge)
# @param i_scSfSvSend sub-vertices adjacent to send sub-cells (faces as bridge)
# @param i_scSfSvRecv sub-vertices adjacent to recv sub-cells (faces as bridge)
# @return three lists (inner, send, recv sub-cells) containing the ids of adjacent sub-cells.
##
def scSfSc( i_scSfSvIn,
            i_scSfSvSend,
            i_scSfSvRecv ):
  # combine all sub-cells
  l_scSfSv = i_scSfSvIn+i_scSfSvSend+i_scSfSvRecv

  # resulting adjacency
  l_scSfSc = []

  # iterate over sub-cells
  for l_c1 in range(len(l_scSfSv)):
    l_scSfSc = l_scSfSc + [[]]

    # iterate over faces io the sub-cell
    for l_f1 in l_scSfSv[l_c1]:
      # fill up undefined (recveive sub-cell)
      if( l_f1 == [-1 for l_en in range(len(l_f1)) ] ):
        l_scSfSc[-1] = l_scSfSc[-1] + [-1]
      else:
        # iterate over all other sub-cells
        for l_c2 in range(len(l_scSfSv)):
          if( l_c1 != l_c2 ):
            # iterate over faces of other sub-cell
            for l_f2 in l_scSfSv[l_c2]:
              # check that this is neither and undefined face and set the vertices match
              if set(l_f1) == set(l_f2):
                l_scSfSc[-1] = l_scSfSc[-1] + [l_c2]

  # split by inner, send, receive
  return l_scSfSc[0:len(i_scSfSvIn)],\
         l_scSfSc[len(i_scSfSvIn):len(i_scSfSvIn) + len(i_scSfSvSend)],\
         l_scSfSc[len(i_scSfSvIn) + len(i_scSfSvSend) : ]

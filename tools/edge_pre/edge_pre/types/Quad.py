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
# Quadrilaterals.
##
class Quad:
  def __init__(self, i_deg):
    # number of dimensions
    self.n_dims = 2

    # number of basis functions
    self.n_basis = (i_deg+1)*(i_deg+1)

    # number of vertices
    self.n_ves = 4

    # number of faces
    self.n_fas = 4

    # number of sub-faces per face
    self.n_sfs = 2*i_deg+1

    # number of sub-cells per element
    self.n_scs = self.n_sfs**2

    # number of inner sub-cells per element
    self.n_scs_in = self.n_scs - 4*self.n_sfs + 4

    # vertices
    self.ves = [ [0,0], [1,0], [1,1], [0,1] ]

    # volume of the reference element
    self.vol = 1

  ##
  # Integration intervals for the element.
  #
  # @param i_syms symbols
  # @return integration intervals for the element.
  ##
  def intEl( self, i_sym ):
    return( [ (i_sym[0], self.ves[0][0], self.ves[1][0]),
              (i_sym[1], self.ves[0][1], self.ves[3][1]) ] )

  ##
  # Determines the local face coordinates of adjacent elements.
  #
  # @param i_symsFa symbols used for the face parametrization.
  # @return list containing one tuple with the new coordinates per face vertex orientation.
  ##
  def faToFa( self, i_symsFa ):
    assert( len(i_symsFa) == 1 )

    l_chi = i_symsFa[0]

    return [ (1-l_chi,) ]

  ##
  # Determines the element coordinates based on the faces coordinates.
  #
  # @param i_symsFa symbols used for the face parametrization.
  # @return list containing one tuple with the element coordinates per face.
  ##
  def faToEl( self, i_symsFa ):
    assert( len(i_symsFa) == 1 )

    l_chi = i_symsFa[0]

    return [ (l_chi,   0       ),
             (0,       l_chi   ),
             (1-l_chi, 1,      ),
             (0,       1-l_chi ) ]
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
# Projection operators.
##
import sympy

##
# Derives the scatter operator (DG -> sub-cell)
#
# @param i_basis basis functions.
# @param i_ints sub-cell integration intervals.
##
def scatter( i_basis,
             i_ints ):
  # scatter matrix
  l_mat = sympy.zeros( len(i_basis), len(i_ints) )

  # compute entries
  for l_co in range(len(i_ints)):
    l_scVol = sympy.integrate( 1, *i_ints[l_co] )
    for l_ro in range(len(i_basis)):
      # intergrate basis for sub-cell
      l_mat[l_ro, l_co] = sympy.integrate( i_basis[l_ro], *i_ints[l_co] )
      # devide by sub-cell volume (integration vs. average)
      l_mat[l_ro, l_co] = l_mat[l_ro, l_co] / l_scVol

  return l_mat

##
# Derives the gather operator (sub-cell -> DG)
#
# @param i_basis basis functions.
# @param i_ints sub-cell integration intervals.
##
def gather( i_basis,
            i_ints ):
  # left matrix Al
  l_aL = sympy.zeros( len(i_basis)+1, len(i_basis)+1 )

  # least squares part of Al
  for l_sc in range(len(i_ints)):
    l_scVol = sympy.integrate( 1, *i_ints[l_sc] )

    for l_ro in range(len(i_basis)):
      for l_co in range(len(i_basis)):
        # sum product over all sub-cells
        l_cont =   sympy.integrate( i_basis[l_ro], *i_ints[l_sc] ) \
                 * sympy.integrate( i_basis[l_co], *i_ints[l_sc] )
        # scale
        l_cont = l_cont * 2
        l_cont = l_cont / l_scVol**2

        # add sub-cell contribution
        l_aL[l_ro, l_co] = l_aL[l_ro, l_co] + l_cont

  # lagrange multiplier part of Al
  for l_rc in range(len(i_basis)):
    # integrate over entire DG element, by adding sub-intervals
    for l_sc in range(len(i_ints)):
      l_aL[ len(i_basis), l_rc ] =  l_aL[ len(i_basis), l_rc ] \
                                   -sympy.integrate( i_basis[l_rc], *i_ints[l_sc] )
    # symmetry for last column
    l_aL[ l_rc, len(i_basis) ] = l_aL[ len(i_basis), l_rc ]

  # right matrix Ar
  l_aR = sympy.zeros( len(i_basis)+1, len(i_ints) )

  # least squares part is scaled scatter matrix
  l_aR[0:len(i_basis), :] = scatter( i_basis, i_ints )
  l_aR = l_aR * 2

  # lagrange multiplier part of Ar
  for l_sc in range(len(i_ints)):
    l_aR[len(i_basis),l_sc] = -sympy.integrate( 1, *i_ints[l_sc] )

  # determine gather operator
  l_gather = l_aL.inv() * l_aR
  l_gather = l_gather[0:len(i_basis), :]

  # derivation is assuming DG DOFs and sub-cell DOFs as column-vectors -> transpose
  l_gather = l_gather.T

  return l_gather

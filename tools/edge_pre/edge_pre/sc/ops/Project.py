##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017-2018, Regents of the University of California
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
# @param i_syms volume symbols.
# @param i_basis basis functions.
# @param i_int reference integration interval.
# @param i_maps mappings from reference element to sub-cells.
# @param i_aDets absolute determinants of the Jacobians.
##
def scatter( i_syms,
             i_basis,
             i_int,
             i_maps,
             i_aDets ):
  # scatter matrix
  l_mat = sympy.zeros( len(i_basis), len(i_maps) )

  # volume of the reference element
  l_refVol = sympy.integrate( 1, *i_int )

  # compute entries
  for l_co in range(len(i_maps)):
    l_scVol = l_refVol * i_aDets[l_co]
    for l_ro in range(len(i_basis)):
      # assemble basis
      l_basis = i_basis[l_ro]

      # assemble substitutions
      l_subs = {}
      for l_di in range( len(i_int) ):
        l_subs[ i_syms[l_di] ] = i_maps[l_co][l_di]

      # substitute
      l_basis = l_basis.subs( l_subs, simultaneous=True )

      # intergrate basis for sub-cell
      l_mat[l_ro, l_co] = sympy.integrate( l_basis, *i_int )
      l_mat[l_ro, l_co] = l_mat[l_ro, l_co] * i_aDets[l_co]
      # devide by sub-cell volume (integration vs. average)
      l_mat[l_ro, l_co] = l_mat[l_ro, l_co] / l_scVol

  return l_mat

##
# Derives the gather operator (sub-cell -> DG)
#
# @param i_syms volume symbols.
# @param i_basis basis functions.
# @param i_int reference integration interval.
# @param i_maps mappings from reference element to sub-cells.
# @param i_aDets absolute determinants of the Jacobians.
##
def gather( i_syms,
            i_basis,
            i_int,
            i_maps,
            i_aDets ):
  # left matrix Al
  l_aL = sympy.zeros( len(i_basis)+1, len(i_basis)+1 )

  # assemble substitutions
  l_subs = []
  for l_sc in range(len(i_maps)):
    l_subs = l_subs + [ {} ]
    for l_di in range( len(i_syms) ):
      l_subs[-1][ i_syms[l_di] ] = i_maps[l_sc][l_di]

  # volume of the reference element
  l_volRef = sympy.integrate( 1, *i_int )

  # least squares part of Al
  for l_sc in range(len(i_maps)):
    l_scVol = l_volRef * i_aDets[l_sc]

    for l_ro in range(len(i_basis)):
      for l_co in range(len(i_basis)):
        # integrands
        l_int0 = i_basis[l_ro].subs( l_subs[l_sc], simultaneous=True )
        l_int1 = i_basis[l_co].subs( l_subs[l_sc], simultaneous=True )

        # sum product over all sub-cells
        l_cont =   sympy.integrate( l_int0, *i_int ) \
                 * sympy.integrate( l_int1, *i_int )
        # scale
        l_cont = l_cont * i_aDets[l_sc]**2
        l_cont = l_cont * 2
        l_cont = l_cont / l_scVol**2

        # add sub-cell contribution
        l_aL[l_ro, l_co] = l_aL[l_ro, l_co] + l_cont

  # lagrange multiplier part of Al
  for l_rc in range(len(i_basis)):
    # integrate over entire DG element, by adding sub-intervals
    for l_sc in range(len(i_maps)):
      l_scInt = sympy.integrate( i_basis[l_rc].subs( l_subs[l_sc], simultaneous=True ), *i_int )
      l_scInt = l_scInt * i_aDets[l_sc]

      l_aL[ len(i_basis), l_rc ] =  l_aL[ len(i_basis), l_rc ] - l_scInt
    # symmetry for last column
    l_aL[ l_rc, len(i_basis) ] = l_aL[ len(i_basis), l_rc ]

  # right matrix Ar
  l_aR = sympy.zeros( len(i_basis)+1, len(i_maps) )

  # least squares part is scaled scatter matrix
  l_aR[0:len(i_basis), :] = scatter( i_syms,\
                                     i_basis,\
                                     i_int,\
                                     i_maps,\
                                     i_aDets )
  l_aR = l_aR * 2

  # lagrange multiplier part of Ar
  for l_sc in range(len(i_maps)):
    l_aR[len(i_basis),l_sc] = -l_volRef * i_aDets[l_sc]

  # determine gather operator
  l_gather = l_aL.inv() * l_aR
  l_gather = l_gather[0:len(i_basis), :]

  # derivation is assuming DG DOFs and sub-cell DOFs as column-vectors -> transpose
  l_gather = l_gather.T

  return l_gather

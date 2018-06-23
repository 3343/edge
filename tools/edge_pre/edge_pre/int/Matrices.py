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
# Derives matrices based on integrations.
##
import sympy
from . import Scalar

##
# Substitues the given symbols in the functions.
#
# @param i_subs list of tuples describing the substitutions.
# @param i_funs functions.
# @return matrix (#subs x #funs) with substituted functions.
##
def subs( i_subs,
          i_funs ):
  # subs matrix
  l_subs = sympy.zeros( len(i_subs), len(i_funs) )

  for l_ro in range(len(i_subs)):
    for l_co in range(len(i_funs)):
      # copy fun
      l_subs[l_ro, l_co] = i_funs[l_co]
      # substitute
      l_subs[l_ro, l_co] = l_subs[l_ro, l_co].subs( i_subs[l_ro][0], i_subs[l_ro][1] )

  return l_subs

##
# Substitutes all symbols in the functions.
#
# @param i_subs list of tuples describing the substitutions.
# @param i_funs functions.
# @return list with substitued functions. All substitutions are applied to every function.
##
def subsAll( i_subs,
             i_funs ):
  l_subs = []

  for l_ba in i_funs:
    # add fun
    l_subs = l_subs + [ sympy.sympify(l_ba) ]

    # apply substitutions
    for l_su in i_subs:
      l_subs[-1] = l_subs[-1].subs( l_su[0], l_su[1] )

  return l_subs

##
# Integrates the functions for multiple integration intervals.
#
# @param i_ints list of lists of integration intervals.
# @param i_funs functions which are integrated.
# @return matrix (#ints x #funs) with integrated functions. Number of rows is the slowest dim in i_ints.
##
def intL( i_ints,
          i_funs ):
  l_intL = sympy.zeros( len(i_ints), len(i_funs) )

  for l_ro in range(len(i_ints)):
    for l_co in range(len(i_funs)):
      # integrate
      l_intL[l_ro, l_co] = Scalar.int( i_funs[l_co], i_ints[l_ro] )[0]

  return l_intL

##
# Computes the mass matrix.
#
# @param i_ints integration intervals.
# @param i_funs functions
##
def mass( i_ints,
          i_funs ):
  # mass matrix
  l_mass = sympy.zeros( len(i_funs), len(i_funs) )

  for l_ro in range(len(i_funs)):
    for l_co in range(len(i_funs)):
      l_prod = i_funs[l_ro] * i_funs[l_co]
      l_mass[l_ro, l_co] = Scalar.int( l_prod, i_ints )[0]

  return l_mass

##
# Computes the stiffness matrices.
#
# @param i_syms symbols, used for the derivative computation.
# @param i_basis element basis.
# @param i_int integration intervals for the element.
#
# @return stiffness matrices, one per given symbol.
##
def stiff( i_syms,
           i_basis,
           i_int ):
  l_stiff = []

  for l_sy in i_syms:
    # derivative computation
    l_der = lambda bf: sympy.diff( bf, l_sy)

    # integration
    l_int = lambda bf: Scalar.int( bf, i_int )[0]

    l_stiff = l_stiff + [[]]
    l_stiff[-1] = sympy.Matrix.transpose( sympy.Matrix( i_basis ) )
    l_stiff[-1] = l_stiff[-1].applyfunc( l_der )
    l_stiff[-1] = sympy.Matrix( i_basis ) * l_stiff[-1]
    l_stiff[-1] = l_stiff[-1].applyfunc( l_int )

  # return stiffness matrices
  return l_stiff

##
# Computes the flux matrices.
#   Remkark: For line elements standard flux matrices are written to the projections.
#
#
# @param i_basisFa face basis.
# @param i_basisEl element basis.
# @param i_symsFa symbols for the faces.
# @param i_symsEl symbols for the element. If only one, special handling for line elements is triggered.
# @param i_intFa integration intervals for the faces.
# @param i_faToFa mapping of face coordinates to adjacent element in dependency of the vertex orientation.
# @param i_faToEl mapping from face-coordinate to volune coordinates in dependency of the face.
#
# @return three lists, 1) matrices for the projection to/from face-basis, 2) face mass-matrix, 3) matrices for neighboring face-orientation and mirroring.
##
def flux( i_basisFa,
          i_basisEl,
          i_symsFa,
          i_symsEl,
          i_intFa,
          i_faToFa,
          i_faToEl ):
  # special handling for line elements (nothing to reduce)
  if( len(i_symsEl) == 1 ):
    assert( i_basisFa     == None )
    assert( i_symsFa      == None )
    assert( i_intFa       == None )
    assert( len(i_faToEl) == 2    )
    assert( i_faToFa      == None )

    # evaluation of basis function ats
    l_eval = [[],[]]
    for l_sd in range(2):
      for l_ba in i_basisEl:
        l_eval[l_sd] = l_eval[l_sd] + [l_ba.subs(i_symsEl[0], i_faToEl[l_sd][0])]
      # convert to sympy matrix
      l_eval[l_sd] = sympy.Matrix( l_eval[l_sd] )

    # pseudo face mass matrix and orientations
    l_massFa  = sympy.eye( 1 )
    l_ori     = [sympy.eye( 1 )]

    # return this as only result
    return l_eval, l_massFa, l_ori

  # derive number of vertices
  l_nVes = len(i_faToFa)

  # derive number of faces
  l_nFas = len(i_faToEl)

  # assemble integrands for projection
  l_proj = []
  for l_fa in range(l_nFas):
    # add new matrix
    l_proj = l_proj + [[]]

    for l_ro in i_basisEl:
      # add a new row
      l_proj[-1] = l_proj[-1] + [[]]

      l_in = l_ro
      # replace symbols
      for l_sy in range( len(i_symsEl) ):
        l_in = l_in.subs( i_symsEl[l_sy], i_faToEl[l_fa][l_sy] )

      for l_co in i_basisFa:
        # add entry
        l_proj[-1][-1] = l_proj[-1][-1] + [l_co*l_in]

  # do the integrations
  for l_fa in range(l_nFas):
    for l_ro in range( len(l_proj[l_fa]) ):
      l_proj[l_fa][l_ro] = list( intL( [i_intFa], l_proj[l_fa][l_ro] ) )
    # create sympy matrix
    l_proj[l_fa] = sympy.Matrix( l_proj[l_fa] )

  # derive face mass matrix
  l_massFa = mass( i_intFa, i_basisFa )

  # multiply with inverse mass matrix
  for l_fa in range(l_nFas):
    l_proj[l_fa] = l_proj[l_fa] * l_massFa.inv()

  # compute matrices accounting for vertex orientation (and mirroring)
  l_ori = []
  for l_ve in range(l_nVes):
    # add new matrix
    l_ori = l_ori + [[]]

    for l_ro in i_basisFa:
      # add row
      l_ori[-1] = l_ori[-1] + [[]]

      # init integrand
      l_in = l_ro

      # replace coords with temporary symbols to avoid replacing conflicts
      for l_sy in range( len(i_symsFa) ):
        l_in = l_in.subs( i_symsFa[l_sy], sympy.symbols('tmp_'+str(l_sy)) )

      for l_sy in range( len(i_symsFa) ):
        l_in = l_in.subs( sympy.symbols('tmp_'+str(l_sy)),  i_faToFa[l_ve][l_sy]  )

      for l_co in i_basisFa:
        # add entry
        l_ori[-1][-1] = l_ori[-1][-1] + [l_in * l_co]

  # perform the integrations
  for l_ve in range(l_nVes):
    for l_ro in range( len(l_ori[l_ve]) ):
      l_ori[l_ve][l_ro] = list( intL( [i_intFa], l_ori[l_ve][l_ro] ) )
    # create sympy matrix
    l_ori[l_ve] = sympy.Matrix( l_ori[l_ve] )

  return l_proj, l_massFa, l_ori
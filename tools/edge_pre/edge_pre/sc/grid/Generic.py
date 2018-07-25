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
# Generic derivation of sub-grid information.
##
import sympy

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

##
# Integration intervals for the sub-cells
#
# @param i_deg polynomial degree.
# @param i_ty element type.
# @param i_syms symbols.
# @param i_shape shape functions.
# @param i_svs sub-vertex coords.
# @param i_scSv sub-vertex ids, adjacent to sub-cells (list of 2: inner, send).
# @param i_scSfScRecv sub-cells adjacent to sub-cells (face as bridge) for recv-elements.
# @return 1) mappings, 2) absolute values of Jacobi determinant.
##
def intSc( i_deg, i_ty, i_syms, i_shape, i_svs, i_scSv, i_scSfScRecv ):
  # get sub-cell coords
  for l_ty in range(2):
    for l_sc in range(len(i_scSv[l_ty])):
      for l_sv in range(len(i_scSv[l_ty][l_sc])):
        i_scSv[l_ty][l_sc][l_sv] = i_svs[ i_scSv[l_ty][l_sc][l_sv] ]

  # derive mappings from reference to sub-cell
  l_maps = [[],[],[]]
  for l_ty in range(2):
    for l_sc in i_scSv[l_ty]:
      l_maps[l_ty] = l_maps[l_ty] + [ refToPhy( l_sc, i_shape ) ]

  # compute determinants of jacobian
  l_aDets = [[],[],[]]
  for l_ty in range(2):
    for l_ma in range(len(l_maps[l_ty])):
      l_aDets[l_ty] = l_aDets[l_ty] + [ absJacDet( i_syms, l_maps[l_ty][l_ma] ) ]

  l_mapsFlat   = l_maps[0] + l_maps[1]
  l_aDetsFlat = l_aDets[0] + l_aDets[1]

  for l_fa in range(i_ty.n_fas):
    # add empty lists for this face
    l_maps[2] = l_maps[2] + [[]]
    l_aDets[2] = l_aDets[2] + [[]]

    for l_sc in range(i_ty.n_sfs):
      l_scId = i_scSfScRecv[l_fa*i_ty.n_sfs + l_sc][0]
      l_maps[2][-1]  = l_maps[2][-1]  + [ l_mapsFlat[ l_scId ]  ]
      l_aDets[2][-1] = l_aDets[2][-1] + [ l_aDetsFlat[ l_scId ] ]

  # simplify mappings
  for l_ty in range( len(l_maps) ):
    for l_ma in range( len(l_maps[l_ty]) ):
      for l_di in range( len(l_maps[l_ty][l_ma]) ):
        l_maps[l_ty][l_ma][l_di] = sympy.simplify( l_maps[l_ty][l_ma][l_di] )

  return l_maps, l_aDets

##
# Derives the mapping from reference coordinates to physical coordinates.
#
# @param i_ves vertices.
# @param i_shapes shape functions.
# @return respective mapping.
##
def refToPhy( i_ves,
              i_shapes ):
  assert( len(i_ves) == len(i_shapes) )

  # create empty mapping
  l_map = [ 0 for _ in range( len(i_ves[0]) ) ]

  for l_di in range(len(l_map)):
    for l_ve in range(len(i_shapes)):
      l_map[l_di] = l_map[l_di] + i_ves[l_ve][l_di] * i_shapes[l_ve]

  return l_map

##
# Computes the absolute value of the Jacobi determinant for the given mapping.
#
# @param i_syms symbols.
# @param i_map mapping.
# @return absolute jacobi determinant.
##
def absJacDet( i_syms,
               i_map ):
  assert( len(i_syms) == len(i_map) )

  # derive jacobian
  l_jac = []

  for l_ma in i_map:
    l_jac = l_jac + [[]]
    for l_sy in i_syms:
      l_jac[-1] = l_jac[-1] + [ sympy.diff( l_ma, l_sy ) ]

  # compute absolute determinant
  l_jac = sympy.Matrix( l_jac )
  l_aDet = sympy.Abs( l_jac.det() )

  return l_aDet
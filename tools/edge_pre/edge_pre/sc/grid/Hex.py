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
# Sub-grids for hexes.
##
import fractions
import edge_pre.types.Hex
from . import Generic
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.backends.backend_pdf
import math

##
# Generates vertices of the sub-grid for the given polynomial degree.
#
# @param i_deg polynomial degree.
# @return list containing the coordinates of the vertices. [*][]: vertex id, [][*]: dimension
##
def svs( i_deg ):
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  # inc from one vertex to the next
  l_inc = [ l_ty.ves[1][0] - l_ty.ves[0][0],
            l_ty.ves[2][1] - l_ty.ves[0][1],
            l_ty.ves[4][2] - l_ty.ves[3][2] ]
  l_inc = [ fractions.Fraction( l_in, l_ty.n_ses ) for l_in in l_inc ]

  # vertices
  l_svs = []
  for l_e3 in range(l_ty.n_ses+1):
    for l_e2 in range(l_ty.n_ses+1):
      for l_e1 in range(l_ty.n_ses+1):
        l_svs = l_svs + [ [ l_ty.ves[0][0]+ l_e1*l_inc[0],
                            l_ty.ves[0][1]+ l_e2*l_inc[1],
                            l_ty.ves[0][2]+ l_e3*l_inc[2] ] ]

  return l_svs

##
# Derives sub-vertices adjacent to the sub-cells.
#
#   The order of the ghost sub-cells follows our convention for DG-faces:
#   Note: The illustration's node-ordering is different from the sub-vertex ordering,
#         strictly following the dimensions.
#   face 0: 0-3-2-1
#   face 1: 0-1-5-4
#   face 2: 1-2-6-5
#   face 3: 3-7-6-2
#   face 4: 0-4-7-3
#   face 5: 4-5-6-7
#
#
#           7 x*******************x 6
#            **                  **
#           * *                 * *
#          *  *                *  *
#         *   *               *   *
#        *    *              *    *
#     4 x*******************x 5   *
#       *     *             *     *
#       *   3 x************ * ****x 2
#       *    *              *    *
#       *   *               *   *
#       |  /                *  *
#  zeta | / eta             * *
#       |/                  **
#       x---****************x
#     0   xi                 1
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the vertex ids of the sub-cells.
##
def scSv( i_deg ):
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  ##
  # Derives four sub-vertices adjacent to a subcell.
  # Order is dimension-wise.
  #
  # @param i_f1 id in first dim.
  # @param i_f2 id in second dim.
  # @param i_f3 id in third dim.
  ##
  def svIds( i_f1, i_f2, i_f3 ):
    l_svIds = []
    l_svIds = l_svIds + [  i_f3    * (l_ty.n_ses+1)**2 +  i_f2    * (l_ty.n_ses+1) + i_f1    ]
    l_svIds = l_svIds + [  i_f3    * (l_ty.n_ses+1)**2 +  i_f2    * (l_ty.n_ses+1) + i_f1 + 1]
    l_svIds = l_svIds + [  i_f3    * (l_ty.n_ses+1)**2 + (i_f2+1) * (l_ty.n_ses+1) + i_f1    ]
    l_svIds = l_svIds + [  i_f3    * (l_ty.n_ses+1)**2 + (i_f2+1) * (l_ty.n_ses+1) + i_f1 + 1]
    l_svIds = l_svIds + [ (i_f3+1) * (l_ty.n_ses+1)**2 +  i_f2    * (l_ty.n_ses+1) + i_f1    ]
    l_svIds = l_svIds + [ (i_f3+1) * (l_ty.n_ses+1)**2 +  i_f2    * (l_ty.n_ses+1) + i_f1 + 1]
    l_svIds = l_svIds + [ (i_f3+1) * (l_ty.n_ses+1)**2 + (i_f2+1) * (l_ty.n_ses+1) + i_f1    ]
    l_svIds = l_svIds + [ (i_f3+1) * (l_ty.n_ses+1)**2 + (i_f2+1) * (l_ty.n_ses+1) + i_f1 + 1]

    # left ghost
    if( i_f1 == -1 ):
      l_svIds[0] = l_svIds[2] = l_svIds[4] = l_svIds[6] = -1
    # right ghost
    if( i_f1 == l_ty.n_ses ):
      l_svIds[1] = l_svIds[3] = l_svIds[5] = l_svIds[7] = -1
    # front ghost
    if( i_f2 == -1 ):
      l_svIds[0] = l_svIds[1] = l_svIds[4] = l_svIds[5] = -1
    # back ghost
    if( i_f2 == l_ty.n_ses ):
      l_svIds[2] = l_svIds[3] = l_svIds[6] = l_svIds[7] = -1
    # bottom ghost
    if( i_f3 == -1 ):
      l_svIds[0] = l_svIds[1] = l_svIds[2] = l_svIds[3] = -1
    # top ghost
    if( i_f3 == l_ty.n_ses ):
      l_svIds[4] = l_svIds[5] = l_svIds[6] = l_svIds[7] = -1

    # move the ghost entries to the end
    l_svEl    = []
    l_svGhost = []
    for l_id in l_svIds:
      if l_id != -1:
        l_svEl =    l_svEl + [l_id]
      else:
        l_svGhost = l_svGhost + [l_id]

    return l_svEl+l_svGhost

  # inner sub-cells
  l_scSvIn = []
  for l_e3 in range( 1, l_ty.n_ses-1):
    for l_e2 in range(1, l_ty.n_ses-1):
      for l_e1 in range(1, l_ty.n_ses-1):
        l_scSvIn = l_scSvIn + [ svIds( l_e1, l_e2, l_e3 ) ]
  
  # send sub-cells
  l_scSvSend = []

  # bottom DG face
  for l_e2 in range(0, l_ty.n_ses):
    for l_e1 in range(0, l_ty.n_ses):
      l_scSvSend = l_scSvSend + [ svIds( l_e1, l_e2, 0 ) ]

  # front DG face
  for l_e3 in range(1, l_ty.n_ses):
    for l_e1 in range(0, l_ty.n_ses):
      l_scSvSend = l_scSvSend + [ svIds( l_e1, 0, l_e3 ) ]

  # right DG face
  for l_e3 in range(1, l_ty.n_ses):
    for l_e2 in range(1, l_ty.n_ses):
      l_scSvSend = l_scSvSend + [ svIds( l_ty.n_ses-1, l_e2, l_e3 ) ]

  # back DG face
  for l_e3 in range(1, l_ty.n_ses):
    for l_e1 in range(0, l_ty.n_ses-1):
      l_scSvSend = l_scSvSend + [ svIds( l_e1, l_ty.n_ses-1, l_e3 ) ]

  # left DG face
  for l_e3 in range(1, l_ty.n_ses):
    for l_e2 in range(1, l_ty.n_ses-1):
      l_scSvSend = l_scSvSend + [ svIds( 0, l_e2, l_e3 ) ]

  # top DG face
  for l_e2 in range(1, l_ty.n_ses-1):
    for l_e1 in range(1, l_ty.n_ses-1):
      l_scSvSend = l_scSvSend + [ svIds( l_e1, l_e2, l_ty.n_ses-1 ) ]

  # receiver sub-cells
  l_scSvRecv = []

  # bottom DG face
  for l_e1 in range(0, l_ty.n_ses):
    for l_e2 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( l_e1, l_e2, -1 ) ]

  # front DG face
  for l_e3 in range(0, l_ty.n_ses):
    for l_e1 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( l_e1, -1, l_e3 ) ]

  # right DG face
  for l_e3 in range(0, l_ty.n_ses):
    for l_e2 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( l_ty.n_ses, l_e2, l_e3 ) ]

  # back DG face
  for l_e1 in range(0, l_ty.n_ses):
    for l_e3 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( l_e1, l_ty.n_ses, l_e3 ) ]

  # left DG face
  for l_e2 in range(0, l_ty.n_ses):
    for l_e3 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( -1, l_e2, l_e3 ) ]

  # top DG face
  for l_e2 in range(0, l_ty.n_ses):
    for l_e1 in range(0, l_ty.n_ses):
      l_scSvRecv = l_scSvRecv + [ svIds( l_e1, l_e2, l_ty.n_ses ) ]

  return l_scSvIn, l_scSvSend, l_scSvRecv

##
# Derives sub-vertices adjacent to sub-cells (faces as brige).
#
# @param i_deg polynomial degree.
##
def scSfSv( i_deg ):
  # determine sub-vertices adjacent to sub-cells
  l_scSv = scSv( i_deg )

  l_scSfSv = [ [], [], [] ]

  # set sub-faces
  for l_ty in range(3):
    for l_sc in l_scSv[l_ty]:
      l_scSfSv[l_ty] = l_scSfSv[l_ty] + [[ [l_sc[0], l_sc[2], l_sc[3], l_sc[1]],
                                           [l_sc[0], l_sc[1], l_sc[5], l_sc[4]],
                                           [l_sc[1], l_sc[3], l_sc[7], l_sc[5]],
                                           [l_sc[2], l_sc[6], l_sc[7], l_sc[3]],
                                           [l_sc[0], l_sc[4], l_sc[6], l_sc[2]],
                                           [l_sc[4], l_sc[5], l_sc[7], l_sc[6]]
                                        ]]

  # reset and sort vertices of ghost faces
  for l_sc in range( len(l_scSfSv[2]) ):
    for l_fa in range(6):
      for l_ve1 in range(4):
        if l_scSfSv[2][l_sc][l_fa][l_ve1] == -1:
          for l_ve2 in range(4):
            l_scSfSv[2][l_sc][l_fa][l_ve2] = -1
      # sort
      l_scSfSv[2][l_sc][l_fa].sort()

  return l_scSfSv[0], l_scSfSv[1], l_scSfSv[2]

##
# Derives sub-cells adjacent to sub-cells (faces as bridge).
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the ids of adjacent sub-cells.
##
def scSfSc( i_deg ):
  # determine sub-vertices adjacent to sub-cells (faces as bridge)
  l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = scSfSv( i_deg )

  # determine sub-cells adjacent to sub-cells (faces as bridge)
  l_scSfScIn, l_scSfScSend, l_scSfScRecv = Generic.scSfSc( l_scSfSvIn,
                                                           l_scSfSvSend,
                                                           l_scSfSvRecv )

  return l_scSfScIn, l_scSfScSend, l_scSfScRecv

##
# Integration intervals for the sub-cells
#
# @param i_deg polynomial degree.
# @param i_syms symbols.
# @return 1) mappings, 2) absolute values of Jacobi determinant.
##
def intSc( i_deg, i_syms ):
  assert( len(i_syms) == 3 )

  # hex shape functions
  l_shape = [ (1 - i_syms[0]) * (1 - i_syms[1]) * (1 - i_syms[2]),
               i_syms[0]      * (1 - i_syms[1]) * (1 - i_syms[2]),
               i_syms[0]      *  i_syms[1]      * (1 - i_syms[2]),
              (1 - i_syms[0]) *  i_syms[1]      * (1 - i_syms[2]),
              (1 - i_syms[0]) * (1 - i_syms[1]) *  i_syms[2],
               i_syms[0]      * (1 - i_syms[1]) *  i_syms[2],
               i_syms[0]      *  i_syms[1]      *  i_syms[2],
              (1 - i_syms[0]) *  i_syms[1]      *  i_syms[2] ]

  # get sub-cell coords
  l_svs = svs( i_deg )
  l_scSv = scSv( i_deg )[0:2]

  # assemble mappings and dets for surface sub-cells
  l_scSfSc = scSfSc( i_deg )

  l_ty = edge_pre.types.Hex.Hex( i_deg )

  return Generic.intSc( i_deg, l_ty, i_syms, l_shape, l_svs, l_scSv[0:2], l_scSfSc[2] )

##
# Derives the integration intervals for the sub-faces at the DG-faces.
#
# @param i_deg degree of the polynomial basis.
# @param i_symsS symbols of the surface integration.
# @param i_symsV symbols of the volume integration.
##
def intSfDg( i_deg, i_symsS, i_symsV ):
  assert( len(i_symsV) == 3)
  assert( len(i_symsS) == 2)


  # hex type
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  # substitutes for the volume coordinates
  l_subs = ( # bottom
             ( ( i_symsV[0], i_symsS[1]     ),
               ( i_symsV[1], i_symsS[0]     ),
               ( i_symsV[2], l_ty.ves[0][2] ) ),
             # front
             ( ( i_symsV[0], i_symsS[0]     ),
               ( i_symsV[1], l_ty.ves[0][1] ),
               ( i_symsV[2], i_symsS[1]     ) ),
             # right
             ( ( i_symsV[0], l_ty.ves[1][0] ),
               ( i_symsV[1], i_symsS[0]     ),
               ( i_symsV[2], i_symsS[1]     ) ),
             # back
             ( ( i_symsV[0], i_symsS[1]     ),
               ( i_symsV[1], l_ty.ves[2][1] ),
               ( i_symsV[2], i_symsS[0]     ) ),
             # left
             ( ( i_symsV[0], l_ty.ves[0][0] ),
               ( i_symsV[1], i_symsS[1]     ),
               ( i_symsV[2], i_symsS[0]     ) ),
             # top
             ( ( i_symsV[0], i_symsS[0]     ),
               ( i_symsV[1], i_symsS[1]     ),
               ( i_symsV[2], l_ty.ves[4][2] ) )
           )

  # increase from one DG-vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  assert(  l_ty.ves[3][1] - l_ty.ves[0][1] == l_inc )
  assert(  l_ty.ves[4][2] - l_ty.ves[0][2] == l_inc )

  # sub-edge length
  l_inc = fractions.Fraction( l_inc, l_ty.n_ses )

  # sub-face intgration intervals
  l_intSf = [ [], [], [], [], [], [] ]

  # set intervals
  assert(  l_ty.ves[0][0] == l_ty.ves[0][1] == l_ty.ves[0][2] == 0 )

  for l_e2 in range(l_ty.n_ses):
    for l_e1 in range(l_ty.n_ses):
      for l_fa in range(6):
        l_intSf[l_fa] = l_intSf[l_fa] + [ [ (i_symsS[0], l_e1*l_inc, (l_e1+1)*l_inc ),
                                            (i_symsS[1], l_e2*l_inc, (l_e2+1)*l_inc ) ] ]

  return l_subs, l_intSf

##
# Derives the types of the sub-cells' faces.
# Orientation is given relative to the face-normal of the DG-element, pointing for third vertex.
#  * Left: Sub-cell normal, pointing to third point, is in the same direction.
#  * Right: only "left" for hexes
#
#  0: face is partial DG-face 0, left
#  1: face is partial DG-face 1, left
#  2: face is partial DG-face 2, left
#  3: face is partial DG-face 3, left
#  4: face is partial DG-face 4, left
#  5: face is partial DG-face 5, left
#
#  6: face is partial DG-face 0, right
#  7: face is partial DG-face 1, right
#  8: face is partial DG-face 2, right
#  9: face is partial DG-face 3, right
# 10: face is partial DG-face 4, right
# 11: face is partial DG-face 5, right
#
# 12: face is inside sub-grid, equivalent to DG-face 0, left
# 13: face is inside sub-grid, equivalent to DG-face 1, left
# 14: face is inside sub-grid, equivalent to DG-face 2, left
# 15: face is inside sub-grid, equivalent to DG-face 3, left
# 16: face is inside sub-grid, equivalent to DG-face 4, left
# 17: face is inside sub-grid, equivalent to DG-face 5, left
#
# 18: face is inside sub-grid, equivalent to DG-face 0, right
# 19: face is inside sub-grid, equivalent to DG-face 1, right
# 20: face is inside sub-grid, equivalent to DG-face 2, right
# 21: face is inside sub-grid, equivalent to DG-face 3, right
# 22: face is inside sub-grid, equivalent to DG-face 4, right
# 23: face is inside sub-grid, equivalent to DG-face 5, right
#
# @return two (inner, send) lists of lists. [*][]: sub-cell, [][*]: face of sub-cell.
##
def scTySf( i_deg ):
  # special handling for degree 0 elements
  if( i_deg == 0 ):
    return [], [ [0,1,2,3,4,5] ]

  # get type
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  # set inner sub-cell types
  l_scTySfIn = [ [12, 13, 14, 15, 16, 17 ] for _ in range(l_ty.n_scs_in) ]

  l_scSfSc = scSfSc(i_deg)[1]

  l_scTySfSend = []

  # iterate over send sub-cells and set types based on adjacency info
  for l_sc1 in l_scSfSc:
    # init with DG
    l_scTySfSend = l_scTySfSend + [ [ 0, 1, 2, 3, 4, 5 ] ]

    # set inner if not a DG sub-face
    for l_sc2 in range(len(l_sc1)):
      if( l_sc1[l_sc2] < l_ty.n_scs ):
        l_scTySfSend[-1][l_sc2] = l_scTySfSend[-1][l_sc2] + 12

  return l_scTySfIn, l_scTySfSend

##
# Define sub-cell reordering based on vertex-combinations, given two DG-faces with adjacent sub-cells.
#
# The reference tetrahedron with counter-clockwise face vertices (w.r.t. opposite face) is given by:
#
#   face 0: 0-3-2-1
#   face 1: 0-1-5-4
#   face 2: 1-2-6-5
#   face 3: 3-7-6-2
#   face 4: 0-4-7-3
#   face 5: 4-5-6-7
#
#
#           7 x*******************x 6
#            **                  **
#           * *                 * *
#          *  *                *  *
#         *   *               *   *
#        *    *              *    *
#     4 x*******************x 5   *
#       *     *             *     *
#       *   3 x************ * ****x 2
#       *    *              *    *
#       *   *               *   *
#       |  /                *  *
#  zeta | / eta             * *
#       |/                  **
#       x---****************x
#     0   xi                 1
#
# The reference storage for a face, in terms of sub-cells, is given by (deg=1):
#
#       .                       .
#      .                       .
#     .                       .
#    #-----------------------*
#    |       |       |       |
#    |   6   |   7   |   8   |
#    |       |       |       |
#    |-------|-------|-------|
#    |       |       |       |
#    |   3   |   4   |   5   |
#    |       |       |       |
#    |-----------------------|
#    |       |       |       |  .
#    |   0   |   1   |   2   | .
#    |       |       |       |.
#    x-----------------------o
#
# Thus, for the first DG-face (0-3-2-1), the vertices are given as:
#   x = 0, o = 3, * = 2, # = 1
#
# Now, if two hexes are adjacent to each other, the returned reordering reflects
#   1) The position of the first vertex (x)
#   2) The change in orientation.
#
# For example, if x is equivalent to o of the second DG-face, the deg=1-reordering reads as:
#
#       .                       .
#      .                       .
#     .                       .
#    *-----------------------#
#    |       |       |       |
#    |   8   |   7   |   6   |
#    |       |       |       |
#    |-------|-------|-------|
#    |       |       |       |
#    |   5   |   4   |   3   |
#    |       |       |       |
#    |-----------------------|
#    |       |       |       |  .
#    |   2   |   1   |   0   | .
#    |       |       |       |.
#    o-----------------------x
#
#    2-1-0-5-4-3-8-7-6
#
# @param i_deg degree.
# @return required reordering, as seen from the adjacent DG-element.
##
def scDgAd( i_deg ):
  # get type
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  # assemble initial sub-cells at the DG-face
  l_scFa = []
  for l_d1 in range( l_ty.n_ses ):
    l_scFa = l_scFa + [[]]
    for l_d0 in range( l_ty.n_ses ):
      l_scFa[-1] = l_scFa[-1] + [ l_d1*l_ty.n_ses + l_d0 ]

  # derive solutions
  l_scDgAd = []

  # transpose based on vertex positions
  for l_ve in range(4):
    if l_ve % 2 == 1:
      l_scDgAd = l_scDgAd + [ l_scFa ]
    else:
      l_scDgAd = l_scDgAd + [ list(map(list, zip(*l_scFa))) ]

  # perform reverse in inner dimension
  for l_ve in [1, 2]:
    l_scDgAd[l_ve] = [ l_in[::-1] for l_in in l_scDgAd[l_ve] ]

  # perfrom reverse in outer dimension
  for l_ve in [2, 3]:
    l_scDgAd[l_ve] = l_scDgAd[l_ve][::-1]

  # zip inner dimension
  for l_ve in range(4):
    l_scDgAd[l_ve] = [ l_en for l_in in l_scDgAd[l_ve] for l_en in l_in ]

  return l_scDgAd

##
# Plots the sub-grid.
#
# @param i_out path to output file.
# @param i_deg degree.
##
def plot( i_out, i_deg ):
  ##
  # Defines a path for the given sub-cell.
  #
  # @param i_ves vertices of the sub-grid.
  # @param i_scSv vertices adjacent to the sub-cell.
  ##
  def path( i_ves, i_scSv ):
    l_code = [ matplotlib.path.Path.MOVETO,
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.CLOSEPOLY ]

    l_path = matplotlib.path.Path( [ i_ves[ i_scSv[0] ][0:2],
                                     i_ves[ i_scSv[1] ][0:2],
                                     i_ves[ i_scSv[3] ][0:2],
                                     i_ves[ i_scSv[2] ][0:2],
                                     i_ves[ i_scSv[0] ][0:2] ],
                                     l_code )
    return l_path

  ##
  # Annotates the plot with given string at the center sub-cell.
  #
  # @param i_str string which will be added to the plot.
  # @param i_svs sub-vertex coordinates.
  # @param i_scSv vertices adjacent to the element.
  # @param i_col color of the string.
  # @param io_plot plot to which the annotation is added.
  ##
  def annotate( i_str, i_svs, i_scSv, i_col, io_plot ):
    l_center = [0, 0]
    for l_sv in i_scSv:
      for l_di in range(2):
        l_center[l_di] = l_center[l_di] + i_svs[l_sv][l_di] / len(i_scSv)

    io_plot.annotate( i_str,
                      xy=l_center,
                      ha='center',
                      va='center',
                      color=i_col )

  ##
  # Checks if the given vertices are part of a single x-y slice.
  #
  # @param i_z id of the slice, 0 is the bottom of the DG-element.
  # @param i_ves ids of sub-vertices.
  # @return true if all sub-vertices are in the slice.
  ##
  def zSlice( i_z, i_svs ):
    # all vertices part of the slice
    l_range = [ i_z * (l_ty.n_ses+1)**2, (i_z+1) * (l_ty.n_ses+1)**2 ]

    # check if any vertex violates the range
    for l_sv in i_svs:
      if l_sv < l_range[0] or l_sv > l_range[1]:
        return False

    return True

  # get info
  l_svs    = svs( i_deg )
  l_scSv   = scSv( i_deg )

  # generate plot
  l_size = max( int( math.sqrt( len( l_scSv[0])+len(l_scSv[1]) ) ) - 2, 3 )

  # get type
  l_ty = edge_pre.types.Hex.Hex( i_deg )

  # output pdf
  l_pdf = matplotlib.backends.backend_pdf.PdfPages( i_out )

  # iterate over z-dimension
  for l_z in range(l_ty.n_ses):
    l_fig = matplotlib.pyplot.figure( figsize=(l_size, l_size) )
    if( i_deg > 1 ):
      matplotlib.pyplot.suptitle( 'Hex is sliced in z-direction. Shown are ids of sub-cells above the slice.' )
    matplotlib.pyplot.title( 'z-slice #'+str(l_z) )

    l_plt = l_fig.add_subplot(111)
    l_plt.set_aspect(1)

    l_id = 0
    for l_in in l_scSv[0]:
      # only plot sub-cell, which are part of the slice
      if zSlice( l_z, l_in[0:4] ):
        l_patch = matplotlib.patches.PathPatch( path(l_svs, l_in), lw=2, facecolor='white' )
        l_plt.add_patch( l_patch )
        annotate( str(l_id), l_svs, l_in, 'black', l_plt )
      l_id = l_id + 1
    for l_se in l_scSv[1]:
      # only plot sub-cell, which are part of the slice
      if zSlice( l_z, l_se[0:4] ):
        l_patch = matplotlib.patches.PathPatch( path(l_svs, l_se), lw=2, facecolor='black', edgecolor='white' )
        l_plt.add_patch( l_patch )
        annotate( str(l_id), l_svs, l_se, 'white', l_plt )
      l_id = l_id + 1

    # save plot
    l_pdf.savefig( l_fig )
  l_pdf.close()
  matplotlib.pyplot.close()

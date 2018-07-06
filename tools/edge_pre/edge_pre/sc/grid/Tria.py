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
# Sub-grids for triangles.
##
import fractions
import sympy
import edge_pre.types.Tria
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from . import Generic
import math

##
# Generates vertices of the sub-grid for the given polynomial degree.
#
# @param i_deg polynomial degree.
# @return list containing the coordinates of the vertices. [*][]: vertex id, [][*]: dimension
##
def svs( i_deg ):
  l_ty = edge_pre.types.Tria.Tria( i_deg )

  # inc from one vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  assert( l_inc == l_ty.ves[2][1] - l_ty.ves[0][0] )
  l_inc = fractions.Fraction( l_inc, l_ty.n_sfs )

  # vertices
  l_svs = []
  for l_f2 in range(l_ty.n_sfs+1):
    for l_f1 in range(l_ty.n_sfs+1 - l_f2):
      l_svs = l_svs + [ [ l_ty.ves[0][0] + l_f1*l_inc,
                          l_ty.ves[0][1] + l_f2*l_inc ] ]

  return l_svs

##
# Derives sub-vertices adjacent to the sub-cells.
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the vertex ids of the sub-cells.
##
def scSv( i_deg ):
  l_ty = edge_pre.types.Tria.Tria( i_deg )

  ##
  # Derives the vertex id of the quad in the triangular quad-grid (counter-clockwise).
  # If the upper right corner is outside the sub-grid, -1 is is returned.
  #
  # @param i_f1 sub-face id in the first dimension (bottom DG-face)
  # @param i_f2 sub-face id in the second dimensions (left DG-face).
  #
  # @return vertices adjacent to the quad; -1 for the upper-right corner if outside the sub-grid
  ##
  def svIdsQuad( i_f1, i_f2 ):
    l_svIds = [i_f1,-1,-1,-1]

    # bottom-left vertex
    for l_d2 in range(i_f2):
      l_svIds[0] = l_svIds[0] + l_ty.n_sfs+1-l_d2

    # bottom-right
    l_svIds[1] = l_svIds[0]+1

    # top-left
    l_svIds[3] = l_svIds[0] +  l_ty.n_sfs-i_f2+1

    # top-right
    if( i_f1 < l_ty.n_sfs-i_f2-1 ): l_svIds[2] = l_svIds[3]+1

    return l_svIds


  # derive inner sub-triangles
  l_scSvIn = []
  for l_f2 in range(0, l_ty.n_sfs):
    for l_f1 in range(0, l_ty.n_sfs-l_f2):
      l_svIdsQuad = svIdsQuad( l_f1, l_f2 )

      # add lower-left triangle if not at at DG-faces
      if( l_f1 > 0 and l_f2 > 0 and l_f1 + l_f2 < l_ty.n_sfs-1 ):
        l_scSvIn = l_scSvIn + [ [ l_svIdsQuad[0], l_svIdsQuad[1], l_svIdsQuad[3] ] ]

      # add upper-right trinagle if it exists
      if( l_svIdsQuad[2] != -1 ):
        l_scSvIn = l_scSvIn + [ [ l_svIdsQuad[1], l_svIdsQuad[2], l_svIdsQuad[3] ] ]


  # derive send sub-triangles
  l_scSvSend = []
  # bottom
  for l_f1 in range(0, l_ty.n_sfs):
    l_svIdsQuad = svIdsQuad( l_f1, 0 )
    l_scSvSend = l_scSvSend + [ [ l_svIdsQuad[0], l_svIdsQuad[1], l_svIdsQuad[3] ] ]

  # right
  for l_f2 in range(1, l_ty.n_sfs):
    for l_f1 in range(0, l_ty.n_sfs-l_f2):
      l_svIdsQuad = svIdsQuad( l_f1, l_f2 )
      if( l_f1 + l_f2 == l_ty.n_sfs-1 ):
        l_scSvSend = l_scSvSend + [ [ l_svIdsQuad[0], l_svIdsQuad[1], l_svIdsQuad[3] ] ]

  # left
  for l_f2 in range(l_ty.n_sfs-2, 0, -1):
    l_svIdsQuad = svIdsQuad( 0, l_f2 )
    l_scSvSend = l_scSvSend + [ [ l_svIdsQuad[0], l_svIdsQuad[1], l_svIdsQuad[3] ] ]


  # derive receive sub-triangles
  l_scSvRecv = []
  # bottom
  for l_f1 in range(0, l_ty.n_sfs):
    l_svIdsQuad = svIdsQuad( l_f1, 0 )
    l_scSvRecv = l_scSvRecv + [ [ l_svIdsQuad[0], l_svIdsQuad[1], -1 ] ]

  # right
  for l_f2 in range(0, l_ty.n_sfs):
    for l_f1 in range(0, l_ty.n_sfs-l_f2):
      l_svIdsQuad = svIdsQuad( l_f1, l_f2 )
      if( l_f1 + l_f2 == l_ty.n_sfs-1 ):
        l_scSvRecv = l_scSvRecv + [ [ l_svIdsQuad[1], l_svIdsQuad[3], -1 ] ]

  # left
  for l_f2 in range(l_ty.n_sfs-1, -1, -1):
    l_svIdsQuad = svIdsQuad( 0, l_f2 )
    l_scSvRecv = l_scSvRecv + [ [ l_svIdsQuad[3], l_svIdsQuad[0], -1 ] ]

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
      l_scSfSv[l_ty] = l_scSfSv[l_ty] + [ [ [l_sc[0], l_sc[1]],
                                            [l_sc[1], l_sc[2]],
                                            [l_sc[2], l_sc[0]] ] ]

  # reset receive-vertices
  for l_sc in range(len(l_scSfSv[2])):
    for l_fa in range(3):
      if( l_scSfSv[2][l_sc][l_fa][0] == -1 or l_scSfSv[2][l_sc][l_fa][1] == -1 ):
        l_scSfSv[2][l_sc][l_fa][0] = -1
        l_scSfSv[2][l_sc][l_fa][1] = -1

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
  assert( len(i_syms) == 2 )

  # tria shape functions
  l_shape = [ 1 - i_syms[0] - i_syms[1],
              i_syms[0],
              i_syms[1] ]

  # get sub-cell coords
  l_svs = svs( i_deg )
  l_scSv = scSv( i_deg )[0:2]

  # assemble mappings and dets for surface sub-cells
  l_scSfSc = scSfSc( i_deg )

  l_ty = edge_pre.types.Tria.Tria( i_deg )

  return Generic.intSc( i_deg, l_ty, i_syms, l_shape, l_svs, l_scSv[0:2], l_scSfSc[2] )

##
# Derives the integration intervals for the sub-faces at the DG-faces.
#
# @param i_deg degree of the polynomial basis.
# @param i_symsS symbols of the surface integration.
# @param i_symsV symbols of the volume integration.
##
def intSfDg( i_deg, i_symsS, i_symsV ):
  assert( len(i_symsV) == 2)
  assert( len(i_symsS) == 1)

  # tria type
  l_ty = edge_pre.types.Tria.Tria( i_deg )

  # substitutes for the volume coordinates
  l_subs = ( # bottom
             ( ( i_symsV[0], i_symsS[0]                ),
               ( i_symsV[1], l_ty.ves[0][1]            ) ),
             # right
             ( ( i_symsV[0], l_ty.ves[1][0]-i_symsS[0] ),
               ( i_symsV[1], i_symsS[0]                ) ),
             # left
             ( ( i_symsV[0], l_ty.ves[0][0]            ),
               ( i_symsV[1], i_symsS[0]                ) )
           )

  # sub-face intgration intervals
  l_intSf = [ [], [], [] ]

  # inc from one vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  assert( l_inc == l_ty.ves[2][1] - l_ty.ves[0][0] )
  l_inc = fractions.Fraction( l_inc, l_ty.n_sfs )

  # set intervals
  for l_sf in range(l_ty.n_sfs):
    # bottom
    l_intSf[0] = l_intSf[0] + [[ ( i_symsS[0],
                                   l_ty.ves[0][0] +  l_sf   *l_inc,
                                   l_ty.ves[0][0] + (l_sf+1)*l_inc ) ]]

    # right
    l_intSf[1] = l_intSf[1] + [[ ( i_symsS[0],
                                   l_ty.ves[0][0] +  l_sf   *l_inc,
                                   l_ty.ves[0][0] + (l_sf+1)*l_inc ) ]]

    # left
    l_intSf[2] = l_intSf[2] + [[ ( i_symsS[0],
                                   l_ty.ves[2][1] - (l_sf+1)*l_inc,
                                   l_ty.ves[2][1] -  l_sf   *l_inc ) ]]

  return l_subs, l_intSf
  
##
# Derives the types of the sub-cells' faces.
# Orientation is given relative to the face-normal of the DG-element, pointing for third vertex.
#  * Left: Sub-cell normal, pointing to third point, is in the same direction.
#  * Right: Sub-cell normal, pointing to third point, is in the opposite direction.
#
#  0: face is partial DG-face 0, left
#  1: face is partial DG-face 1, left
#  2: face is partial DG-face 2, left
#
#  3: face is partial DG-face 0, right
#  4: face is partial DG-face 1, right
#  5: face is partial DG-face 2, right
#
#  6: face is inside sub-grid, parallel to DG-face 0, left
#  7: face is inside sub-grid, parallel to DG-face 1, left
#  8: face is inside sub-grid, parallel to DG-face 2, left
#
#  9: face is inside sub-grid, parallel to DG-face 0, right
# 10: face is inside sub-grid, parallel to DG-face 1, right
# 11: face is inside sub-grid, parallel to DG-face 2, right
#
# @return two (inner, send) lists of lists. [*][]: sub-cell, [][*]: face of sub-cell.
##
def scTySf( i_deg ):
  # special handling for degree 0 elements
  if( i_deg == 0 ):
    return [], [ [0,1,2] ]

  # get vertices adjacent to the sub-cells
  l_scSvIn, l_scSvSend, l_scSvRecv = scSv( i_deg )

  # set inner sub-cell types
  l_scTySfIn = []
  for l_scSv in l_scSvIn:
    # ascending same orientation as DG-element
    if l_scSv[1] < l_scSv[2]:
      l_scTySfIn = l_scTySfIn + [[6, 7, 8]]
    # descending is "mirrored" at second DG-face
    else:
      l_scTySfIn = l_scTySfIn + [[11, 9, 10]]

  # set send sub-cell types
  l_ty = edge_pre.types.Tria.Tria( i_deg )
  l_scTySfSend = []

  # bottom-left corner
  l_scTySfSend = l_scTySfSend + [[0, 7, 2]]

  # bottom, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[0, 7, 8]]

  # bottom right corner
  l_scTySfSend = l_scTySfSend + [[0, 1, 8]]

  # "right", no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[6, 1, 8]]

  # top-left corner
  l_scTySfSend = l_scTySfSend + [[6, 1, 2]]

  # left, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[6, 7, 2]]

  return l_scTySfIn, l_scTySfSend

##
# Define sub-cell reordering based on vertex-combinations, given two DG-faces with adjacent sub-cells.
#   For two-dimensional elements, there is only one possible combination.
#   The order is simply reversed.
#
# @param i_deg degree.
# @return required reordering, as seen from the adjacent DG-element.
##
def scDgAd( i_deg ):
  # get type
  l_ty = edge_pre.types.Tria.Tria( i_deg )

  l_scDgAd = list( range(l_ty.n_sfs) )
  l_scDgAd.reverse()

  return [ l_scDgAd ]

##
# Plots the sub-grid.
#
# @param i_out path to output file.
# @param i_svs vertices of the sub-grid.
# @param i_scSvIn sub-vertices adjancent to inner sub-cells.
# @param i_scSvSend sub-vertices adjacent to send sub-cells.
##
def plot( i_out, i_svs, i_scSvIn, i_scSvSend ):
  l_size = max( int( math.sqrt( len( i_scSvIn)+len(i_scSvSend) ) ) - 1, 3 )
  l_fig = matplotlib.pyplot.figure( figsize=(l_size, l_size) )

  l_plt = l_fig.add_subplot(111)
  l_plt.set_aspect(1)

  ##
  # Defines a path for the given sub-cell.
  #
  # @param i_svs vertices of the sub-grid.
  # @param i_scSv vertices adjacent to the sub-cell.
  ##
  def path( i_ves, i_scSv ):
    l_code = [ matplotlib.path.Path.MOVETO,
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.CLOSEPOLY ]

    l_path = matplotlib.path.Path( [ i_svs[ i_scSv[0] ],
                                     i_svs[ i_scSv[1] ],
                                     i_svs[ i_scSv[2] ],
                                     i_svs[ i_scSv[0] ] ],
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

  l_id = 0
  for l_in in i_scSvIn:
    l_patch = matplotlib.patches.PathPatch( path(i_svs, l_in), lw=2, facecolor='white' )
    l_plt.add_patch( l_patch )
    annotate( str(l_id), i_svs, l_in, 'black', l_plt )
    l_id = l_id + 1
  for l_se in i_scSvSend:
    l_patch = matplotlib.patches.PathPatch( path(i_svs, l_se), lw=2, facecolor='black', edgecolor='white' )
    l_plt.add_patch( l_patch )
    annotate( str(l_id), i_svs, l_se, 'white', l_plt )
    l_id = l_id + 1

  # save plot
  l_fig.savefig( i_out )
  matplotlib.pyplot.close()

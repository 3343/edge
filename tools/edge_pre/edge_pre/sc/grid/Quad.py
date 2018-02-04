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
# Sub-grids for quads.
##
import fractions
import edge_pre.types.Quad
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
  l_ty = edge_pre.types.Quad.Quad( i_deg )

  # inc from one vertex to the next
  l_inc = [ l_ty.ves[1][0] - l_ty.ves[0][0],
            l_ty.ves[3][1] - l_ty.ves[0][1] ]
  l_inc = [ fractions.Fraction( l_in, l_ty.n_sfs ) for l_in in l_inc ]

  # vertices
  l_svs = []
  for l_f2 in range(l_ty.n_sfs+1):
    for l_f1 in range(l_ty.n_sfs+1):
      l_svs = l_svs + [ [ l_ty.ves[0][0]+ l_f1*l_inc[0],
                          l_ty.ves[0][1]+ l_f2*l_inc[1] ] ]

  return l_svs

##
# Derives sub-vertices adjacent to the sub-cells.
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the vertex ids of the sub-cells.
##
def scSv( i_deg ):
  l_ty = edge_pre.types.Quad.Quad( i_deg )

  ##
  # Derives four sub-vertices adjacent to a subcell.
  # Order is counter-clockwise
  #
  # @param i_f1 id in first dim.
  # @param i_f2 id in second dim.
  ##
  def svIds( i_f1, i_f2 ):
    l_svIds = []
    l_svIds = l_svIds + [  i_f2    * (l_ty.n_sfs+1) + i_f1    ]
    l_svIds = l_svIds + [  i_f2    * (l_ty.n_sfs+1) + i_f1 + 1]
    l_svIds = l_svIds + [ (i_f2+1) * (l_ty.n_sfs+1) + i_f1 + 1]
    l_svIds = l_svIds + [ (i_f2+1) * (l_ty.n_sfs+1) + i_f1    ]

    return l_svIds

  # inner sub-cells
  l_scSvIn = []
  for l_f2 in range(1, l_ty.n_sfs-1):
    for l_f1 in range(1, l_ty.n_sfs-1):
      l_scSvIn = l_scSvIn + [ svIds( l_f1, l_f2 ) ]

  # send sub-cells
  l_scSvSend = []

  # first DG face
  for l_f1 in range(0, l_ty.n_sfs):
    l_scSvSend = l_scSvSend +  [ svIds( l_f1, 0 ) ]

  #second DG face
  for l_f2 in range(1, l_ty.n_sfs):
    l_scSvSend = l_scSvSend +  [ svIds( l_ty.n_sfs-1, l_f2 ) ]

  # third DG face
  for l_f1 in range(l_ty.n_sfs-2, -1, -1):
    l_scSvSend = l_scSvSend +  [ svIds( l_f1, l_ty.n_sfs-1 ) ]

  # fourth DG face
  for l_f2 in range(l_ty.n_sfs-2,  0, -1):
    l_scSvSend = l_scSvSend +  [ svIds( 0, l_f2 ) ]


  # receive sub-cells
  l_scSvRecv = []

  # first DG face
  for l_f1 in range(0, l_ty.n_sfs):
    l_scSvRecv = l_scSvRecv + [ [l_f1, l_f1+1, -1, -1] ]

  # second DG face
  for l_f2 in range(0, l_ty.n_sfs):
    l_scSvRecv = l_scSvRecv + [ [(l_f2+1)*(l_ty.n_sfs+1) - 1,
                                 (l_f2+2)*(l_ty.n_sfs+1) - 1,
                                  -1, -1 ] ]

  # third DG face
  for l_f1 in range(l_ty.n_sfs-1, -1, -1):
    l_scSvRecv = l_scSvRecv + [ [ l_f1+1 + (l_ty.n_sfs+1)*l_ty.n_sfs,
                                  l_f1   + (l_ty.n_sfs+1)*l_ty.n_sfs,
                                 -1, -1] ]

  # fourth DG face
  for l_f2 in range(l_ty.n_sfs-1,  -1, -1):
    l_scSvRecv = l_scSvRecv + [ [(l_f2+1)*(l_ty.n_sfs+1),
                                 (l_f2  )*(l_ty.n_sfs+1),
                                  -1, -1 ] ]

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
                                            [l_sc[2], l_sc[3]],
                                            [l_sc[3], l_sc[0]] ] ]

  # reset receive-vertices
  for l_sc in range(len(l_scSfSv[2])):
    for l_fa in range(4):
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
# @param i_syms symbols
# @return integration intervals for inner sub-cells, send sub-cells, and per DG-face of surface sub-cells.
##
def intSc( i_deg, i_syms ):
  assert( len(i_syms) == 2 )

  # get vertex coordinates
  l_svs = svs( i_deg )

  # get vertices adjacent to sub-cells
  l_scSv = scSv( i_deg )

  l_intSc = [[],[]]

  # iterate over inner and send
  for l_ty in range(2):
    # iterate over sub-cells
    for l_sc in l_scSv[l_ty]:
      l_intSc[l_ty] = l_intSc[l_ty] + [[]]

      l_intSc[l_ty][-1] = l_intSc[l_ty][-1] + [ ( i_syms[0],
                                                  l_svs[ l_sc[0] ][0],
                                                  l_svs[ l_sc[1] ][0], ) ]
      l_intSc[l_ty][-1] = l_intSc[l_ty][-1] + [ ( i_syms[1],
                                                  l_svs[ l_sc[0] ][1],
                                                  l_svs[ l_sc[3] ][1], ) ]

  # determine integration intervals for sub-cells adjacent to the DG surface
  l_ty = edge_pre.types.Quad.Quad( i_deg )

  # inc from one vertex to the next
  l_inc = [ l_ty.ves[1][0] - l_ty.ves[0][0],
            l_ty.ves[3][1] - l_ty.ves[0][1] ]
  l_inc = [ fractions.Fraction( l_in, l_ty.n_sfs ) for l_in in l_inc ]

  # sub-cells adjacent ot the surface
  l_intSurf = [[], [], [], []]

  # bottom
  for l_sc in range(l_ty.n_sfs):
    l_intSurf[0] = l_intSurf[0] + [ [ ( i_syms[0], l_sc*l_inc[0], (l_sc+1)*l_inc[0] ),
                                      ( i_syms[1], 0,             l_inc[1]          ) ] ]

  # right
  for l_sc in range(l_ty.n_sfs):
    l_intSurf[1] = l_intSurf[1] + [ [ ( i_syms[0], (l_ty.n_sfs-1)*l_inc[0],  l_ty.n_sfs*l_inc[0] ),
                                      ( i_syms[1],  l_sc         *l_inc[1], (l_sc+1)   *l_inc[1] ) ] ]

  # top
  for l_sc in range(l_ty.n_sfs-1, -1, -1):
    l_intSurf[2] = l_intSurf[2] + [ [ ( i_syms[0],  l_sc         *l_inc[0], (l_sc+1)   *l_inc[0] ),
                                      ( i_syms[1], (l_ty.n_sfs-1)*l_inc[1],  l_ty.n_sfs*l_inc[1] ) ] ]

  # left
  for l_sc in range(l_ty.n_sfs-1, -1, -1):
    l_intSurf[3] = l_intSurf[3] + [ [ ( i_syms[0],  0,                      l_inc[0]             ),
                                      ( i_syms[1],  l_sc         *l_inc[1], (l_sc+1)   *l_inc[1] ) ] ]
                                    

  return l_intSc[0], l_intSc[1], l_intSurf

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

  # quad type
  l_ty = edge_pre.types.Quad.Quad( i_deg )

  # substitutes for the volume coordinates
  l_subs = ( # bottom
             ( ( i_symsV[0], i_symsS[0]     ),
               ( i_symsV[1], l_ty.ves[0][1] ) ),
             # right
             ( ( i_symsV[0], l_ty.ves[1][0] ),
               ( i_symsV[1], i_symsS[0]     ) ),
             # top
             ( ( i_symsV[0], i_symsS[0]     ),
               ( i_symsV[1], l_ty.ves[3][1] ) ),
             # left
             ( ( i_symsV[0], l_ty.ves[0][0] ),
               ( i_symsV[1], i_symsS[0]     ) )
           )

  # sub-face intgration intervals
  l_intSf = [ [], [], [], [] ]

  # inc from one vertex to the next
  l_inc = [ l_ty.ves[1][0] - l_ty.ves[0][0],
            l_ty.ves[3][1] - l_ty.ves[0][1] ]
  l_inc = [ fractions.Fraction( l_in, l_ty.n_sfs ) for l_in in l_inc ]

  # set intervals
  for l_sf in range(l_ty.n_sfs):
    # bottom
    l_intSf[0] = l_intSf[0] + [[ ( i_symsS[0],
                                   l_ty.ves[0][0] +  l_sf   *l_inc[0],
                                   l_ty.ves[0][0] + (l_sf+1)*l_inc[0] ) ]]

    # right
    l_intSf[1] = l_intSf[1] + [[ ( i_symsS[0],
                                   l_ty.ves[1][1] +  l_sf   *l_inc[1],
                                   l_ty.ves[1][1] + (l_sf+1)*l_inc[1] ) ]]

    # top
    l_intSf[2] = l_intSf[2] + [[ ( i_symsS[0],
                                   l_ty.ves[2][0] - (l_sf+1)*l_inc[0],
                                   l_ty.ves[2][0] -  l_sf   *l_inc[0] ) ]]

    # left
    l_intSf[3] = l_intSf[3] + [[ ( i_symsS[0],
                                   l_ty.ves[3][1] - (l_sf+1)*l_inc[1],
                                   l_ty.ves[3][1] -  l_sf   *l_inc[1] ) ]]

  return l_subs, l_intSf
  

##
# Derives the types of the sub-cells' faces.
# Orientation is given relative to the face-normal of the DG-element, pointing for third vertex.
#  * Left: Sub-cell normal, pointing to third point, is in the same direction.
#  * Right: only "left" for quads
#
#  0: face is partial DG-face 0, left
#  1: face is partial DG-face 1, left
#  2: face is partial DG-face 2, left
#  3: face is partial DG-face 3, left
#
#  4: face is partial DG-face 0, right
#  5: face is partial DG-face 1, right
#  6: face is partial DG-face 2, right
#  7: face is partial DG-face 3, right
#
#  8: face is inside sub-grid, equivalent to DG-face 0, left
#  9: face is inside sub-grid, equivalent to DG-face 1, left
# 10: face is inside sub-grid, equivalent to DG-face 2, left
# 11: face is inside sub-grid, equivalent to DG-face 3, left
#
# 12: face is inside sub-grid, equivalent to DG-face 0, right
# 13: face is inside sub-grid, equivalent to DG-face 1, right
# 14: face is inside sub-grid, equivalent to DG-face 2, right
# 15: face is inside sub-grid, equivalent to DG-face 3, right
#
# @return two (inner, send) lists of lists. [*][]: sub-cell, [][*]: face of sub-cell.
##
def scTySf( i_deg ):
  # get type
  l_ty = edge_pre.types.Quad.Quad( i_deg )

  # set inner sub-cell types
  l_scTySfIn = [ [8, 9, 10, 11] for _ in range(l_ty.n_scs_in) ]

  # set send sub-cell types
  l_scTySfSend = []

  # bottom-left corner
  l_scTySfSend = l_scTySfSend + [[0, 9, 10, 3]]

  # bottom, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[0, 9, 10, 11]]

  # bottom right corner
  l_scTySfSend = l_scTySfSend + [[0, 1, 10, 11]]

  # right, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[8, 1, 10, 11]]

  # top-right corer
  l_scTySfSend = l_scTySfSend + [[8, 1, 2, 11]]

  # top, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[8, 9, 2, 11]]

  # top-left coner
  l_scTySfSend = l_scTySfSend + [[8, 9, 2, 3]]

  # left, no corner
  for l_sc in range(1, l_ty.n_sfs-1):
    l_scTySfSend = l_scTySfSend + [[8, 9, 10, 3]]

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
  l_ty = edge_pre.types.Quad.Quad( i_deg )

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
  l_size = max( int( math.sqrt( len( i_scSvIn)+len(i_scSvSend) ) ) - 2, 3 )
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
               matplotlib.path.Path.LINETO,
               matplotlib.path.Path.CLOSEPOLY ]

    l_path = matplotlib.path.Path( [ i_svs[ i_scSv[0] ],
                                     i_svs[ i_scSv[1] ],
                                     i_svs[ i_scSv[2] ],
                                     i_svs[ i_scSv[3] ],
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

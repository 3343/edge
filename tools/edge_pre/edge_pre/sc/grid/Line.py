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
# Sub-grids for line elements.
##
import fractions
import edge_pre.types.Line
from . import Generic
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

##
# Generates vertices of the sub-grid for the given polynomial degree.
#
# @param i_deg polynomial degree.
# @return list containing the coordinates of the vertices. [*][]: vertex id, [][*]: dimension
##
def svs( i_deg ):
  l_ty = edge_pre.types.Line.Line( i_deg )

  # inc from one vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  l_inc = fractions.Fraction( l_inc, l_ty.n_scs )

  # vertices
  l_svs = [ [ l_ty.ves[0][0] ] ]
  while( l_svs[-1][0] != l_ty.ves[1][0] ):
    l_svs = l_svs + [ [ l_svs[-1][0] + l_inc ] ]

  return l_svs

##
# Derives sub-vertices adjacent to the sub-cells.
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the vertex ids of the sub-cells.
##
def scSv( i_deg ):
  # special handling for degree 0
  if( i_deg == 0 ):
    return [], [ [0,1] ], [ [-1,0], [1, -1] ]

  l_ty = edge_pre.types.Line.Line( i_deg )

  # inner-sub-cells
  l_scSv = []
  for l_sc in range(1, l_ty.n_scs-1):
    l_scSv = l_scSv + [ [l_sc, l_sc+1] ]

  # send sub-cells
  l_scSvSend = [ [0,1], [l_ty.n_scs-1, l_ty.n_scs] ]

  # recv sub-cells
  l_scSvRecv = [ [0,-1], [l_ty.n_scs, -1] ]

  return l_scSv, l_scSvSend, l_scSvRecv

##
# Derives sub-cells adjacent to sub-cells (faces as bridge).
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the ids of adjacent sub-cells.
##
def scSfSc( i_deg ):
  # special handling for deg 0
  if( i_deg == 0 ):
    return [], [ [1, 2] ], [ [0, -1], [0, -1] ]

  l_ty = edge_pre.types.Line.Line( i_deg )

  l_scSfScIn = []
  # inner-sub-cells
  for l_sc in range(0, l_ty.n_scs-2):
    l_scSfScIn = l_scSfScIn + [ [l_sc-1, l_sc+1] ]

  # adjust inner-send adjacency
  l_scSfScIn[0][0]  = l_ty.n_scs-2
  l_scSfScIn[-1][1] = l_ty.n_scs-1

  # send-sub-cells
  l_scSfScSend = [ # 1st remote
                   [ l_ty.n_scs, 0 ],
                   # 2nd remote
                   [ l_ty.n_scs-3, l_ty.n_scs+1 ]
                 ]

  # receive-sub-cells
  l_scSfScRecv = [ # 1st local
                   [ l_ty.n_scs-2, -1],
                   # 2nd local
                   [ l_ty.n_scs-1, -1]
                 ]

  return l_scSfScIn, l_scSfScSend, l_scSfScRecv

##
# Derives the types of the sub-cells' faces.
# Orientation is given relative to the face-normal of the DG-element, pointing for third vertex.
#  * Left: Sub-cell normal, pointing to third point, is in the same direction.
#  * Right: only "left" for line elements
#
#  0: face is partial DG-face 0, left
#  1: face is partial DG-face 1, left
#
#  2: face is partial DG-face 0, right
#  3: face is partial DG-face 1, right
#
#  4: face is inside sub-grid, parallel to DG-face 0, left
#  5: face is inside sub-grid, parallel to DG-face 1, left
#
#  6: face is inside sub-grid, parallel to DG-face 0, right
#  7: face is inside sub-grid, parallel to DG-face 1, right
#
# @return two (inner, send) lists of lists. [*][]: sub-cell, [][*]: face of sub-cell.
##
def scTySf( i_deg ):
  # special handling for degree 0 elements
  if( i_deg == 0 ):
    return [], [ [0, 1] ]

  # get type
  l_ty = edge_pre.types.Line.Line( i_deg )

  # inner
  l_scTySfIn = [ [4,5] for _ in range(l_ty.n_scs-2) ]

  # send
  l_scTySfSend = [ [0, 5], [4, 1] ]

  return l_scTySfIn, l_scTySfSend

##
# Integration intervals for the sub-cells
#
# @param i_deg polynomial degree.
# @param i_syms symbols.
# @return 1) mappings, 2) absolute values of Jacobi determinant.
##
def intSc( i_deg, i_syms ):
  assert( len(i_syms) == 1 )

  # line shape functions
  l_shape = [ 1 - i_syms[0],
              i_syms[0] ]

  # get sub-cell coords
  l_svs = svs( i_deg )
  l_scSv = scSv( i_deg )[0:2]

  # assemble mappings and dets for surface sub-cells
  l_scSfSc = scSfSc( i_deg )

  l_ty = edge_pre.types.Line.Line( i_deg )

  return Generic.intSc( i_deg, l_ty, i_syms, l_shape, l_svs, l_scSv[0:2], l_scSfSc[2] )

##
# Plots the sub-grid.
#
# @param i_out path to output file.
# @param i_svs vertices of the sub-grid.
# @param i_scSvIn sub-vertices adjancent to inner sub-cells.
# @param i_scSvSend sub-vertices adjacent to send sub-cells.
##
def plot( i_out, i_svs, i_scSvIn, i_scSvSend ):
  l_fig = matplotlib.pyplot.figure()
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

    l_path = matplotlib.path.Path( [ i_svs[ i_scSv[0] ] + [0],
                                     i_svs[ i_scSv[1] ] + [0],
                                     i_svs[ i_scSv[1] ] + [0.1],
                                     i_svs[ i_scSv[0] ] + [0.1],
                                     i_svs[ i_scSv[0] ] + [0] ],
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
    l_center = i_svs[ i_scSv[0] ][0] + i_svs[ i_scSv[1] ][0]
    l_center = l_center * 0.5

    io_plot.annotate( i_str,
                      xy=[l_center, 0.05],
                      ha='center',
                      va='center',
                      rotation=90,
                      color=i_col )

  l_id = 0
  for l_in in i_scSvIn:
    l_patch = matplotlib.patches.PathPatch( path(i_svs, l_in), lw=2, facecolor='white' )
    l_plt.add_patch( l_patch )
    annotate( str(l_id), i_svs, l_in, 'black', l_plt )
    l_id = l_id + 1
  for l_se in i_scSvSend:
    l_patch = matplotlib.patches.PathPatch( path(i_svs, l_se), lw=2, facecolor='black' )
    l_plt.add_patch( l_patch )
    annotate( str(l_id), i_svs, l_se, 'white', l_plt )
    l_id = l_id + 1

  # save plot
  l_fig.savefig( i_out )
  matplotlib.pyplot.close()

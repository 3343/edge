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
# Sub-grids for tets.
##
from fractions import Fraction as Fra
import sympy
import edge_pre.types.Tet
from . import Generic

##
# Generates vertices of the sub-grid for the given polynomial degree.
#
# @param i_deg polynomial degree.
# @return list containing the coordinates of the vertices. [*][]: vertex id, [][*]: dimension
##
def svs( i_deg ):
  l_ty = edge_pre.types.Tet.Tet( i_deg )

  # inc from one vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  assert( l_inc == l_ty.ves[2][1] - l_ty.ves[0][1] )
  assert( l_inc == l_ty.ves[3][2] - l_ty.ves[0][2] )
  l_inc = Fra( l_inc, l_ty.n_ses )

  # vertices
  l_svs = []
  for l_e3 in range(l_ty.n_ses+1):
    for l_e2 in range(l_ty.n_ses+1 - l_e3):
      for l_e1 in range(l_ty.n_ses+1 - l_e3 - l_e2):
        l_svs = l_svs + [ [ l_ty.ves[0][0] + l_e1*l_inc,
                            l_ty.ves[0][1] + l_e2*l_inc,
                            l_ty.ves[0][2] + l_e3*l_inc ] ]

  return l_svs

##
# Derives the sub-tets in an element, separated by their tet-type.
#
# @param i_deg polynomial degree.
# @return 1) vertices of sub-cells per hex, sorted dimension-wise. Per entry by (up, down, H1, H2, H3, H4), empty list if not part of the element.
#         2) (i,j,k) -> mono indices.
##
def subTets( i_deg ):
  l_ty = edge_pre.types.Tet.Tet( i_deg )

  # get sub-vertices
  l_svs = svs( i_deg )

  # assemble mono-indices
  l_mono = {}
  l_id = 0
  for l_e3 in range(l_ty.n_ses+1):
    for l_e2 in range(l_ty.n_ses+1 - l_e3):
      for l_e1 in range(l_ty.n_ses+1 - l_e3 - l_e2):
        l_mono[ (l_e1, l_e2, l_e3) ] = l_id
        l_id = l_id + 1


  ##
  # Upward-pointing tetrahedron.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /|        /|
  #   / |       / |
  #  /__.______/  |
  #  3  | _____|__|
  #  |  2      |  /
  #  | /       | /
  #  0_________1/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetUp( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k,   i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k,   i_l+1) ] ]  ]
    return l_ves

  ##
  # Downward-pointing tetrahedron.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /3        /1
  #   / |       / |
  #  /__.______0  |
  #  |  | _____|__2
  #  |  /      |  /
  #  | /       | /
  #  |/________|/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetDo( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k+1, i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l+1) ] ]  ]
    return l_ves

  ##
  # First tet in octahedron holes.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /3        /|
  #   / |       / |
  #  /__.______/  |
  #  |  2 _____|__1
  #  |  /      |  /
  #  | /       | /
  #  |/________0/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetH1( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l+1) ] ]  ]
    return l_ves

  ##
  # Second tet in octahedron holes.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /2        /|
  #   / |       / |
  #  /__.______/  |
  #  3  | _____1__|
  #  |  /      |  /
  #  | /       | /
  #  |/________0/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetH2( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k,   i_l+1) ] ]  ]
    return l_ves

  ##
  # Third tet in octahedron holes.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /1        /|
  #   / |       / |
  #  /__.______/  |
  #  2  0 _____|__|
  #  |  /      |  /
  #  | /       | /
  #  |/________3/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetH3( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k,   i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l  ) ] ]  ]
    return l_ves

  ##
  # Fourth tet in octahedron holes.
  #   Vertex ordering is counterclockwise for faces if watching from outside-to-inside.
  #   Faces: 0-2-1, 0-1-3, 0-3-2, 1-2-3
  #     __________
  #    /2        /|
  #   / |       / |
  #  /__.______/  |
  #  |  | _____3__1
  #  |  /      |  /
  #  | /       | /
  #  |/________0/
  #
  # Reference:
  #   Dumbser, M., and Raphael L.
  #   "A simple robust and accurate a posteriori sub-cell finite volume limiter for the discontinuous Galerkin method on unstructured meshes."
  #   Journal of Computational Physics 319 (2016): 163-199.
  #
  # @param i_j index j.
  # @param i_k index k.
  # @param i_l index l.
  # @return vertices of the sub-cell.
  ##
  def tetH4( i_j,
             i_k,
             i_l ):
    l_ves = []
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k+1, i_l  ) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j,   i_k+1, i_l+1) ] ]  ]
    l_ves = l_ves + [  l_svs[ l_mono[ (i_j+1, i_k,   i_l+1) ] ]  ]
    return l_ves

  # generate sub-cell tets
  l_subTets = []

  for l_e3 in range(l_ty.n_ses):
    for l_e2 in range(l_ty.n_ses - l_e3):
      for l_e1 in range(l_ty.n_ses - l_e3 - l_e2):
        l_subTets = l_subTets + [ [ tetUp(l_e1, l_e2, l_e3) ] ]
    
        if(l_e1+1, l_e2+1, l_e3+1) in l_mono.keys():
          l_subTets[-1] = l_subTets[-1] + [ tetDo(l_e1, l_e2, l_e3) ]
        else:
          l_subTets[-1] = l_subTets[-1] + [ [] ]
        if(l_e1+1, l_e2+1, l_e3) in l_mono.keys():
          l_subTets[-1] = l_subTets[-1] + [ tetH1(l_e1, l_e2, l_e3) ]
          l_subTets[-1] = l_subTets[-1] + [ tetH2(l_e1, l_e2, l_e3) ]
          l_subTets[-1] = l_subTets[-1] + [ tetH3(l_e1, l_e2, l_e3) ]
          l_subTets[-1] = l_subTets[-1] + [ tetH4(l_e1, l_e2, l_e3) ]
        else:
          l_subTets[-1] = l_subTets[-1] + [ [], [], [], [] ]
  return l_subTets, l_mono

##
# Derives sub-vertices adjacent to the sub-cells.
#
# @param i_deg polynomial degree.
# @return three lists (inner, send, recv sub-cells) containing the vertex ids of the sub-cells.
##
def scSv( i_deg ):
  # get sub-vertices
  l_svs = svs( i_deg )

  # get sub-tets
  l_subTets = subTets( i_deg )[0]

  # determine vertices at the DG surface
  l_svsBnd = [ [], [], [], [] ]
  for l_sv in l_svs:
    # fist face
    if( l_sv[2] == 0 ):
      l_svsBnd[0] = l_svsBnd[0] + [l_sv]
    # second face
    if( l_sv[1] == 0 ):
      l_svsBnd[1] = l_svsBnd[1] + [l_sv]
    # third face
    if( l_sv[0] == 0 ):
      l_svsBnd[2] = l_svsBnd[2] + [l_sv]
    # fourth face
    if( l_sv[0] == sympy.sympify(1)-l_sv[1]-l_sv[2]):
      l_svsBnd[3] = l_svsBnd[3] + [l_sv]

  # make sure that the sizes match
  assert( len(l_svsBnd[0]) == len(l_svsBnd[1]) )
  assert( len(l_svsBnd[0]) == len(l_svsBnd[2]) )
  assert( len(l_svsBnd[0]) == len(l_svsBnd[3]) )

  ##
  # Determines the DG-faces of a sub-cell
  #
  # @return list of indices of DG-faces, empty if inner-sub-cell
  ##
  def dgFaces( i_svs ):
    l_dgFaces = []

    # iterate over DG-faces
    for l_dg in range(4):
      # count number of sub-vertices, which are part of the DG-face
      l_count = 0
      for l_sv in i_svs:
        if l_sv in l_svsBnd[l_dg]:
          l_count = l_count+1
      # three vertices per DG-face at most
      assert( l_count <= 3 )

      # add DG-face, if shared with sub-cell
      if( l_count == 3):
        l_dgFaces = l_dgFaces + [ l_dg ]
    return l_dgFaces

  ##
  # Gets the vertex ids of the given sub-vertices.
  #
  # @param i_svCoords coordinates of the vertices.
  # @return ids of the vertices.
  ##
  def getSvIds( i_svCoords ):
    l_ids = []
    for l_sv in i_svCoords:
      assert( l_sv in l_svs )
      l_ids = l_ids + [ l_svs.index( l_sv ) ]
    return l_ids

  # define empty inner, send- and receive sub-cells
  l_scSvIn, l_scSvSendFa = [], [ [], [], [], [] ]

  # determine inner and send
  for l_he in l_subTets:
    for l_sc in l_he:
      # ignore non-existent tets
      if l_sc == []: continue

      l_dgFaces = dgFaces( l_sc )

      # if empty: inner
      if len(l_dgFaces) == 0:
        l_scSvIn = l_scSvIn + [ getSvIds(l_sc) ]
        getSvIds( l_sc )
      # otherwise send
      else:
        for l_fa in l_dgFaces:
          l_scSvSendFa[l_fa] = l_scSvSendFa[l_fa] + [ getSvIds(l_sc) ]

  # check sizes
  assert( i_deg == 0 or
          len(l_scSvIn) ==   (2*i_deg+1)**3 -
                           4*(2*i_deg+1)**2 + # DG-faces
                           6*(2*i_deg+1)    - # edges
                              4 )             # vertices
  for l_fa in range(4):
    assert( len(l_scSvSendFa[l_fa]) == (2*i_deg+1)**2 )

  # unify send sub-cells (removes duplicates at edges and vertices)
  l_scSvSend = []
  for l_fa in l_scSvSendFa:
    for l_sc in l_fa:
      if l_sc not in l_scSvSend:
        l_scSvSend = l_scSvSend + [l_sc]

  # determine receive sub-cells
  l_scSvRecv = []
  for l_fa in range(4):
    for l_se in l_scSvSendFa[l_fa]:
      l_scSvRecv = l_scSvRecv + [[]]
      for l_sv in l_se:
        if( l_svs[l_sv] in l_svsBnd[l_fa] ):
          l_scSvRecv[-1] = l_scSvRecv[-1] + [l_sv]
      # check for all three vertices
      assert( len(l_scSvRecv[-1]) == 3)
      # add ghost entry for remaining vertex
      l_scSvRecv[-1] = l_scSvRecv[-1] + [-1]
      # exchange vertex 1 and 2 for DG-face 1 and 3
      # 1) we loose the tet's vertex 2 (0-2-1 storage for sub-face instead of 0-1-3)
      # 3) we loose the tet's vertex 0 (0-2-1 storage for sub-face instead of 1-2-3)
      if l_fa == 1 or l_fa == 3:
        l_tmp = l_scSvRecv[-1][1]
        l_scSvRecv[-1][1] = l_scSvRecv[-1][2]
        l_scSvRecv[-1][2] = l_tmp

  # perform transpose of hex-order for first and third face, following the element's face coords
  l_ty = edge_pre.types.Tet.Tet( i_deg )

  for l_fa in [0,2]:
    l_trans = []
    for l_e1 in range(l_ty.n_ses):
      # initial offset given by the preceeding faces
      l_off = (l_fa * l_ty.n_sfs)
      for l_e2 in range(l_ty.n_ses-l_e1):
        l_trans = l_trans + [ l_scSvRecv[l_e1*2+l_off] ]
        if( l_e2 != l_ty.n_ses-l_e1-1):
          l_trans = l_trans + [ l_scSvRecv[l_e1*2+l_off+1] ]

        # increase offset according to the transpose op
        l_off = (l_off + (l_ty.n_ses-l_e2)*2 - 1)
      
    # assign transpose
    l_scSvRecv[ (l_fa * l_ty.n_sfs):((l_fa+1) * l_ty.n_sfs) ] = l_trans

  # perform normalization for fourth face, currently storing in h4-up-h4-up-up instead of up-h4-up-h4-up
  # initial offset given by the preceeding faces
  l_off = (3 * l_ty.n_sfs)
  for l_e1 in range(l_ty.n_ses):
    for l_e2 in range( l_ty.n_ses-l_e1-1 ):
      l_tmp = l_scSvRecv[l_off]
      l_scSvRecv[l_off] = l_scSvRecv[l_off+1]
      l_scSvRecv[l_off+1] = l_tmp
      l_off = l_off+2
    l_off = l_off+1

  return l_scSvIn, l_scSvSend, l_scSvRecv

##
# Derives sub-vertices adjacent to sub-cells (faces as bridge).
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
      l_scSfSv[l_ty] = l_scSfSv[l_ty] + [ [ [l_sc[0], l_sc[2], l_sc[1]],
                                            [l_sc[0], l_sc[1], l_sc[3]],
                                            [l_sc[0], l_sc[3], l_sc[2]],
                                            [l_sc[1], l_sc[2], l_sc[3]] ] ]

  # reset receive-vertices
  for l_sc in range(len(l_scSfSv[2])):
    for l_fa in range(4):
      if -1 in l_scSfSv[2][l_sc][l_fa]:
        l_scSfSv[2][l_sc][l_fa][0] = -1
        l_scSfSv[2][l_sc][l_fa][1] = -1
        l_scSfSv[2][l_sc][l_fa][2] = -1

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

  # Tetrahedral shape functions
  l_shape = [ 1 - i_syms[0] - i_syms[1] - i_syms[2],
              i_syms[0],
              i_syms[1],
              i_syms[2] ]
  
  # get sub-cell coords
  l_svs = svs( i_deg )
  l_scSv = scSv( i_deg )[0:2]

  # assemble mappings and dets for surface sub-cells
  l_scSfSc = scSfSc( i_deg )

  l_ty = edge_pre.types.Tet.Tet( i_deg )

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

  # tet type
  l_ty = edge_pre.types.Tet.Tet( i_deg )

  # substitutes for the volume coordinates
  l_subs = ( ( ( i_symsV[0], i_symsS[1]              ),
               ( i_symsV[1], i_symsS[0]              ),
               ( i_symsV[2], 0                       ) ),
             ( ( i_symsV[0], i_symsS[0]              ),
               ( i_symsV[1], 0                       ),
               ( i_symsV[2], i_symsS[1]              ) ),
             ( ( i_symsV[0], 0                       ),
               ( i_symsV[1], i_symsS[1]              ),
               ( i_symsV[2], i_symsS[0]              ) ),
             ( ( i_symsV[0], 1-i_symsS[0]-i_symsS[1] ),
               ( i_symsV[1], i_symsS[0]              ),
               ( i_symsV[2], i_symsS[1]              ) ),
           )

  # sub-face intgration intervals
  l_intSf = []

  # inc from one vertex to the next
  l_inc = l_ty.ves[1][0] - l_ty.ves[0][0]
  assert( l_inc == l_ty.ves[2][1] - l_ty.ves[0][1] )
  assert( l_inc == l_ty.ves[3][2] - l_ty.ves[0][2] )
  l_inc = Fra( l_inc, l_ty.n_ses )

  # set intervals
  for l_e2 in range(l_ty.n_ses):
    for l_e1 in range(l_ty.n_ses-l_e2):
      l_intSf = l_intSf + [ [ ( i_symsS[0],
                                l_e1*l_inc,
                                (l_e1+1)*l_inc - (i_symsS[1]-l_e2*l_inc) ),
                              ( i_symsS[1],
                                l_e2*l_inc,
                                (l_e2+1)*l_inc ) ] ]

      # add second triangle only if not the last quad
      if( l_e1 != l_ty.n_ses-l_e2-1 ):
        l_intSf = l_intSf + [ [ ( i_symsS[0],
                                  (l_e1+1)*l_inc - (i_symsS[1]-l_e2*l_inc) ,
                                  (l_e1+1)*l_inc ),
                                ( i_symsS[1],
                                  l_e2*l_inc,
                                  (l_e2+1)*l_inc ) ] ]

  # check that we got all sub-triangles
  assert( len(l_intSf) == l_ty.n_sfs )

  return l_subs, [ l_intSf, l_intSf, l_intSf, l_intSf ]
  
##
# Derives the types of the sub-cells' faces.
#
# DG-face 0 has normal 0,0,1
# DG-face 1 has normal 0,1,0
# DG-face 2 has normal 1,0,0
# DG-Face 3 has normal sqrt(1/3), sqrt(1/3), sqrt(1/3)
#
# Sub-tet types H1-H4 introduce additional normals.
# In this case, "left" is the side which points to (0,0,0).
# H1:
#   (DG)  0,0,1               - left
#   (DG)  0,1,0               - right
#   (new) sqrt(2), sqrt(2), 0 - right
#   (new) sqrt(2), 0, sqrt(2) - left
#
# H2:
#   (DG)  0,0,1               - right
#   (DG)  0,1,0               - left
#   (new) sqrt(2), sqrt(2), 0 - left
#   (new) sqrt(2), 0, sqrt(2) - right
#
# H3:
#   (DG)  1,0,0                           - left
#   (DG)  sqrt(1/3), sqrt(1/3), sqrt(1/3) - right
#   (new) sqrt(1/2),sqrt(1/2),0           - left
#   (new) sqrt(2), 0, sqrt(2)             - left
#
# H4:
#   (DG)  1,0,0                           - right
#   (DG)  sqrt(1/3), sqrt(1/3), sqrt(1/3) - left
#   (new) sqrt(1/2), sqrt(1/2), 0         - right
#   (new) sqrt(1/2), 0, sqrt(2)           - right
#
# Additional face 0 (add-face 0): ( sqrt(2), sqrt(2), 0       )
# Additional face 1 (add-face 1): ( sqrt(2), 0,       sqrt(2) )
#
# Orientation is given relative to the face-normal of the DG-element, pointing to fourth vertex.
#  * Left: Sub-cell normal, pointing to fourth point, is in the same direction.
#  * Right: Sub-cell normal, pointing to fourth point, is in the opposite direction.
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
#  8: face is inside sub-grid, parallel to DG-face  0, left
#  9: face is inside sub-grid, parallel to DG-face  1, left
# 10: face is inside sub-grid, parallel to DG-face  2, left
# 11: face is inside sub-grid, parallel to DG-face  3, left
# 12: face is inside sub-grid, parallel to add-face 0, left.
# 13: face is inside sub-grid, parallel to add-face 1, left.
#
# 14: face is inside sub-grid, parallel to DG-face  0, right
# 15: face is inside sub-grid, parallel to DG-face  1, right
# 16: face is inside sub-grid, parallel to DG-face  2, right
# 17: face is inside sub-grid, parallel to DG-face  3, right
# 18: face is inside sub-grid, parallel to add-face 0, right
# 19: face is inside sub-grid, parallel to add-face 1, right
# #
# # @return two (inner, send) lists of lists. [*][]: sub-cell, [][*]: face of sub-cell.
# ##
def scTySf( i_deg ):
  # special handling for degree 0 elements
  if( i_deg == 0 ):
    return [], [ [0,1,2,3] ]

  # get sub-vertices
  l_svs = svs( i_deg )

  # get sub-faces
  l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = scSfSv( i_deg )

  # define normals
  l_normals = []
  l_normals = l_normals + [ sympy.Matrix( [0,0,1] )                ]
  l_normals = l_normals + [ sympy.Matrix( [0,1,0] )                ]
  l_normals = l_normals + [ sympy.Matrix( [1,0,0] )                ]
  l_normals = l_normals + [ sympy.Matrix( [-sympy.sqrt(Fra(1,3)),
                                           -sympy.sqrt(Fra(1,3)),
                                           -sympy.sqrt(Fra(1,3))] ) ]
  l_normals = l_normals + [ sympy.Matrix( [-sympy.sqrt(Fra(1,2)),
                                           -sympy.sqrt(Fra(1,2)),
                                            0] )                   ]
  l_normals = l_normals + [ sympy.Matrix( [-sympy.sqrt(Fra(1,2)),
                                            0,
                                           -sympy.sqrt(Fra(1,2))] ) ]

  # define empty sub-cell types
  l_scTySf = [[], []]

  # iterate over sub-faces
  for l_ty in range(2):
    for l_sc in [l_scSfSvIn, l_scSfSvSend][l_ty]:
      # add empty list
      l_scTySf[l_ty] = l_scTySf[l_ty] + [[]]

      for l_sf in l_sc:
        # assemble three points
        l_or = sympy.Matrix( l_svs[ l_sf[0] ] )
        l_p0 = sympy.Matrix( l_svs[ l_sf[1] ] )
        l_p1 = sympy.Matrix( l_svs[ l_sf[2] ] )

        # assemble edges
        l_e0 = l_p0 - l_or
        l_e1 = l_p1 - l_or

        # compute normalized cross product
        l_cr = l_e0.cross( l_e1 ).normalized()

        # determine type of normal
        for l_no in range(len(l_normals)):
          if l_cr   == l_normals[l_no]:
            l_ri = 1
            break
          elif l_cr == -l_normals[l_no]:
            l_ri = 0
            break
          # make sure we found the correct normal
          assert( l_no != len(l_normals)-1 )

        # determine if this a partial DG face
        l_pd = 0

        if l_ty == 1:
          # bottom
          if(   l_or[2] == 0 and l_p0[2] == 0 and l_p1[2] == 0 ):
            assert( l_no == 0 )
            l_pd = 1
          elif( l_or[1] == 0 and l_p0[1] == 0 and l_p1[1] == 0 ):
            assert( l_no == 1 )
            l_pd = 1
          elif( l_or[0] == 0 and l_p0[0] == 0 and l_p1[0] == 0 ):
            assert( l_no == 2)
            l_pd = 1
          elif( l_or[0] == sympy.sympify(1)-l_or[1]-l_or[2] and
                l_p0[0] == sympy.sympify(1)-l_p0[1]-l_p0[2] and
                l_p1[0] == sympy.sympify(1)-l_p1[1]-l_p1[2] ):
            assert( l_no == 3 )
            l_pd = 1

        # add sub-cell type
        if l_pd == 0:
          l_scTySf[l_ty][-1] = l_scTySf[l_ty][-1] + [ 8 + l_ri * 6 +  l_no ]
        else:
          assert( l_no < 4 )
          l_scTySf[l_ty][-1] = l_scTySf[l_ty][-1] + [ 0 + l_ri * 4 +  l_no ]

  return l_scTySf[0], l_scTySf[1]

##
# Define sub-cell reordering based on vertex-combinations, given two DG-faces with adjacent sub-cells.
#
# @param i_deg degree.
# @return required reordering, as seen from the adjacent DG-element.
##
def scDgAd( i_deg ):
  # get type
  l_ty = edge_pre.types.Tet.Tet( i_deg )

  l_qgrid = []

  # assemble pseudo quad grid, holding the default face ordering
  l_id = 0
  for l_e2 in range(l_ty.n_ses):
    # add quad-row to grid
    l_qgrid = l_qgrid + [[]]

    for l_e1 in range(l_ty.n_ses):
      # assemble this quad, empty list if no triangles present
      l_quad = []

      # add upward-pointing tria
      if l_e1 < l_ty.n_ses - l_e2:
        l_quad = l_quad + [l_id]
        l_id = l_id+1
      # add downward-pointing tria
      if l_e1 < l_ty.n_ses - l_e2 - 1:
        l_quad = l_quad + [l_id]
        l_id = l_id+1

      # add quad to grid
      l_qgrid[-1] = l_qgrid[-1] + [ l_quad ]
  
  # assemble the different reorderings, based on the faces' vertices lying on top of each other
  l_scDgAd = [ [], [], [] ]

  for l_e1 in range(l_ty.n_ses):
    for l_e2 in range(l_ty.n_ses):
      l_scDgAd[0] = l_scDgAd[0] + l_qgrid[l_e2][l_e1]

  for l_e2 in range(l_ty.n_ses):
    for l_e1 in range(l_ty.n_ses):
      l_scDgAd[1] = l_scDgAd[1] + list( reversed( l_qgrid[l_e2][l_ty.n_ses - l_e1 - 1] ) )

  # create a flat list for the last reordering
  l_fgrid = []
  for l_e2 in range(l_ty.n_ses):
    l_fgrid = l_fgrid + [ [ l_en for l_ro in l_qgrid[l_e2] for l_en in l_ro ] ]
  l_fgrid.reverse()

  # continue until all values have been added
  while l_fgrid[-1] != []:
    for l_e1 in range(l_ty.n_ses):
      l_newVals = []
      # add at most two values per column
      for _ in range(2):
        if l_fgrid[l_e1] != []:
          l_newVals = l_newVals + [ l_fgrid[l_e1][-1] ]
          l_fgrid[l_e1] = l_fgrid[l_e1][0:-1]
      
      # add reversed
      l_newVals.reverse()
      l_scDgAd[2] = l_scDgAd[2] + l_newVals

  return l_scDgAd
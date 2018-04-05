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
# Unit tests for tet sub-grids.
##
import unittest
from fractions import Fraction as Fra
from . import Tet
import sympy

# Reference grid:
#
# Elements are sorted by:
#   1) inner, send, receive
#   2) face if send or receive (bottom, front, right, back, left, top)
#   3) dimension x, y, z
#
# 2nd order
#
#   slice 3 ........x 
#                   * *    .      
#             layer *    *      .
#             2     *       *         x    
#   slice 2 ........*----------*    .   .
#                   *          |  *       .
#             layer *          |.   *       .
#             1     *         .|       *      .
#   slice 1 ........*------.------------- *     .
#                   *    .     |          |  *    .
#             layer *  .       |          |     *   .
#             0     *.         |          |        *
#   slice 0 ........x---------------------------------x
#
#
#     slice 0 
#
#          9
#          | *
#          |   *
#          |     *
#          |       *
#          |         *
#          7-----------8
#          | *         | *
#          |   *       |   *
#          |     *     |     *
#          |       *   |       *
#          |         * |         *
#          4-----------5-----------6
#          | *         | *         | *
#          |   *       |   *       |   *
#          |     *     |     *     |     *
#          |       *   |       *   |       *
#          |         * |         * |         *
#          0-----------1-----------2-----------3
#
#
#     slice 1
#
#          15
#          | *
#          |   *
#          |     *
#          |       *
#          |         *
#          13---------14
#          | *         | *
#          |   *       |   *
#          |     *     |     *
#          |       *   |       *
#          |         * |         *
#          10---------11----------12
#
#
#     slice 2
#
#          18
#          | *
#          |   *
#          |     *
#          |       *
#          |         *
#          16---------17
#
#
#    slice 3 (single point)
#
#          19
#
#
class TestGridTet( unittest.TestCase ):
  ##
  # Tests generation of vertices.
  ##
  def test_svs(self):
    # FV
    l_svs = Tet.svs( 0 )
    self.assertEqual( l_svs, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ] )

    # 2nd order
    l_svs = Tet.svs( 1 )
    self.assertEqual( l_svs, [ [ 0,        0,        0        ],
                               [ Fra(1,3), 0,        0        ],
                               [ Fra(2,3), 0,        0        ],
                               [ Fra(3,3), 0,        0        ],

                               [ 0,        Fra(1,3), 0        ],
                               [ Fra(1,3), Fra(1,3), 0        ],
                               [ Fra(2,3), Fra(1,3), 0        ],

                               [ 0,        Fra(2,3), 0        ],
                               [ Fra(1,3), Fra(2,3), 0        ],

                               [ 0,        Fra(3,3), 0        ],


                               [ 0,        Fra(0,3), Fra(1,3) ],
                               [ Fra(1,3), Fra(0,3), Fra(1,3) ],
                               [ Fra(2,3), Fra(0,3), Fra(1,3) ],

                               [ 0,        Fra(1,3), Fra(1,3) ],
                               [ Fra(1,3), Fra(1,3), Fra(1,3) ],

                               [ 0,        Fra(2,3), Fra(1,3) ],


                               [ 0,        Fra(0,3), Fra(2,3) ],
                               [ Fra(1,3), Fra(0,3), Fra(2,3) ],
                               [ Fra(0,3), Fra(1,3), Fra(2,3) ],

                               [ Fra(0,3), Fra(0,3), Fra(3,3) ] ] )

  #l_tmp = Tet.svs(2) 
  #l_tmp = [[float(x), float(y), float(z)] for x, y, z in l_tmp]
  #pprint.pprint( l_tmp )

  ##
  # Tests scSv adjacency.
  ##
  def test_scSv(self):
    # l_tmp = Tet.scSv( 2 )
    # l_tmp = l_tmp[0] + l_tmp[1]
    # pprint.pprint( l_tmp )

    # first order sub-grid
    l_scSvIn, l_scSvSend, l_scSvRecv = Tet.scSv( 0 )

    self.assertEqual( l_scSvIn, [] )

    self.assertEqual( l_scSvSend, [ [0, 1, 2, 3] ] )

    self.assertEqual( l_scSvRecv, [ [0, 1, 2, -1],
                                    [0, 3, 1, -1],
                                    [0, 2, 3, -1],
                                    [1, 3, 2, -1] ] )

    # second order sub-grid
    l_scSvIn, l_scSvSend, l_scSvRecv = Tet.scSv( 1 )

    self.assertEqual( l_scSvIn, [ [11, 14, 5, 13],
                                  [1, 5, 13, 11],
                                  [5, 14, 11, 2],
                                  [5, 14, 15, 13],
                                  [11, 14, 13, 18] ] )

    self.assertEqual( l_scSvSend, [ [0, 1, 4, 10],
                                    [1, 5, 4, 13],
                                    [1, 2, 5, 11],
                                    [2, 6, 5, 14],
                                    [2, 3, 6, 12],
                                    [4, 5, 7, 13],
                                    [5, 8, 7, 15],
                                    [5, 6, 8, 14],
                                    [7, 8, 9, 15],

                                    [1,  11, 13, 10],
                                    [2,  12, 14, 11],
                                    [10, 11, 13, 16],
                                    [11, 17, 18, 16],
                                    [11, 12, 14, 17],
                                    [16, 17, 18, 19],

                                    [4, 13, 10, 1],
                                    [7, 15, 13, 5],
                                    [13, 18, 16, 11],
                                    [13, 14, 15, 18],

                                    [2, 6, 14, 12],
                                    [5, 8, 15, 14],
                                    [11, 14, 18, 17] ] )

    # Imposed reordering for faces 0 and 2:
    #   0 - 0
    #   1 - 1
    #   2 - 5
    #   3 - 6
    #   4 - 8
    #   5 - 2
    #   6 - 3
    #   7 - 7
    #   8 - 4
    #
    # Imposed reordering for face 3:
    #   0 - 1
    #   1 - 0
    #   2 - 3
    #   3 - 2
    #   4 - 4
    #   5 - 6
    #   6 - 5
    #   7 - 7
    #   8 - 8
    self.assertEqual( l_scSvRecv, [ [0, 1, 4, -1], # 0
                                    [1, 5, 4, -1], # 1
                                    [4, 5, 7, -1], # 5
                                    [5, 8, 7, -1], # 6
                                    [7, 8, 9, -1], # 8
                                    [1, 2, 5, -1], # 2
                                    [2, 6, 5, -1], # 3
                                    [5, 6, 8, -1], # 7
                                    [2, 3, 6, -1], # 4

                                    [0, 10, 1,  -1],
                                    [1, 10, 11, -1],
                                    [1, 11,  2, -1],
                                    [2, 11, 12, -1],
                                    [2, 12,  3, -1],

                                    [10, 16, 11, -1],
                                    [11, 16, 17, -1],
                                    [11, 17, 12, -1],
                                    [16, 19, 17, -1],

                                    [0, 4, 10, -1],   # 0
                                    [4, 13, 10, -1],  # 1
                                    [10, 13, 16, -1], # 5
                                    [13, 18, 16, -1], # 6
                                    [16, 18, 19, -1], # 8
                                    [4, 7, 13, -1],   # 2
                                    [7, 15, 13, -1],  # 3
                                    [13, 15, 18, -1], # 7
                                    [7, 9, 15, -1],   # 4

                                    [3, 12,  6, -1],     # 1
                                    [6, 12, 14, -1],     # 0
                                    [6, 14,  8, -1],     # 3
                                    [8, 14, 15, -1],     # 2
                                    [8,  15,  9, -1],    # 4
                                    [12, 17, 14, -1],    # 6
                                    [14, 17, 18, -1],    # 5
                                    [14, 18, 15, -1],    # 7
                                    [17, 19, 18, -1] ] ) # 8

  ##
  # Tests scSfSv adjacency
  ##
  def test_scSfSv(self):
    # second order
    l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = Tet.scSfSv( 1 )

    # check inner sub-cells
    self.assertEqual( l_scSfSvIn, [ [ [11,  5, 14], [11, 14, 13], [11, 13,  5], [14,  5, 13] ],
                                    [ [ 1, 13,  5], [ 1,  5, 11], [ 1, 11, 13], [ 5, 13, 11] ],
                                    [ [ 5, 11, 14], [ 5, 14,  2], [ 5,  2, 11], [14, 11,  2] ],
                                    [ [ 5, 15, 14], [ 5, 14, 13], [ 5, 13, 15], [14, 15, 13] ],
                                    [ [11, 13, 14], [11, 14, 18], [11, 18, 13], [14, 13, 18] ] ] )

    # check a few send sub-cells
    self.assertEqual( l_scSfSvSend[8], [ [7, 9, 8], [7, 8, 15], [7, 15, 9], [8, 9, 15] ] )
    self.assertEqual( l_scSfSvSend[15], [ [4, 10, 13], [4, 13, 1], [4, 1, 10], [13, 10, 1] ] )
    self.assertEqual( l_scSfSvSend[17], [ [13, 16, 18], [13, 18, 11], [13, 11, 16], [18, 16, 11] ] )

    # check receive sub-cells
    #
    # Imposed reordering for faces 0 and 2:
    #   0 - 0
    #   1 - 1
    #   2 - 5
    #   3 - 6
    #   4 - 8
    #   5 - 2
    #   6 - 3
    #   7 - 7
    #   8 - 4
    #
    # Imposed reordering for face 3:
    #   0 - 1
    #   1 - 0
    #   2 - 3
    #   3 - 2
    #   4 - 4
    #   5 - 6
    #   6 - 5
    #   7 - 7
    #   8 - 8
    self.assertEqual( l_scSfSvRecv, [ [ [ 0,  4,  1], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 0
                                      [ [ 1,  4,  5], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 1
                                      [ [ 4,  7,  5], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 5
                                      [ [ 5,  7,  8], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 6
                                      [ [ 7,  9,  8], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 8
                                      [ [ 1,  5,  2], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 2
                                      [ [ 2,  5,  6], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 3
                                      [ [ 5,  8,  6], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 7
                                      [ [ 2,  6,  3], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 4

                                      [ [ 0, 1,  10], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [ 1, 11, 10], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [ 1,  2, 11], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [ 2, 12, 11], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [ 2,  3, 12], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],

                                      [ [10, 11, 16], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [11, 17, 16], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [11, 12, 17], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],
                                      [ [16, 17, 19], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],

                                      [ [ 0, 10,  4], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 0
                                      [ [ 4, 10, 13], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 1
                                      [ [10, 16, 13], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 5
                                      [ [13, 16, 18], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 6
                                      [ [16, 19, 18], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 8
                                      [ [ 4, 13,  7], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 2
                                      [ [ 7, 13, 15], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 3
                                      [ [13, 18, 15], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 7
                                      [ [ 7, 15,  9], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ], # 4

                                      [ [ 3, 6,  12], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 1
                                      [ [ 6, 14, 12], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 0
                                      [ [ 6,  8, 14], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 3
                                      [ [ 8, 15, 14], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 2
                                      [ [ 8,  9, 15], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 4
                                      [ [12, 14, 17], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 6
                                      [ [14, 18, 17], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 5
                                      [ [14, 15, 18], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ],    # 7
                                      [ [17, 18, 19], [-1, -1, -1], [-1, -1, -1], [-1, -1, -1] ] ] ) # 8


  ##
  # Tests scSfSc adjacency.
  ##
  def test_scSfSc(self):
    # second order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Tet.scSfSc( 1 )

    # check first inner sub-cell
    self.assertEqual( l_scSfScIn[0], [2, 4, 1, 3] )

    self.assertEqual( l_scSfScSend[-1], [4, 18, 17, 27+3*9+6] )

    # check a few receive sub-cells
    # Imposed reordering for faces 0 and 2:
    #   0 - 0
    #   1 - 1
    #   2 - 5
    #   3 - 6
    #   4 - 8
    #   5 - 2
    #   6 - 3
    #   7 - 7
    #   8 - 4
    self.assertEqual( l_scSfScRecv[0:18], [ [ 5, -1, -1, -1],  # 0
                                            [ 6, -1, -1, -1 ], # 1
                                            [ 10, -1, -1, -1 ], # 5
                                            [ 11, -1, -1, -1 ], # 6
                                            [ 13, -1, -1, -1 ], # 8
                                            [ 7, -1, -1, -1 ], # 2
                                            [ 8, -1, -1, -1 ], # 3
                                            [ 12, -1, -1, -1 ], # 7
                                            [ 9, -1, -1, -1 ], # 4

                                            [  5, -1, -1, -1 ],
                                            [ 14, -1, -1, -1 ],
                                            [  7, -1, -1, -1 ],
                                            [ 15, -1, -1, -1 ],
                                            [ 9, -1, -1, -1 ],
                                            [ 16, -1, -1, -1 ],
                                            [ 17, -1, -1, -1 ],
                                            [ 18, -1, -1, -1 ],
                                            [ 19, -1, -1, -1 ] ] )

  ##
  # Tests integration intervals for sub-cells.
  ##
  def test_intSc(self):
    l_xi0 = sympy.symbols('xi_0')
    l_xi1 = sympy.symbols('xi_1')
    l_xi2 = sympy.symbols('xi_2')

    #
    # check that we obtain the same integration results as for the entire tet
    #
    for l_de in range(0, 5):
      l_intRef = [ (l_xi0, 0, 1-l_xi1-l_xi2), (l_xi1, 0, 1-l_xi2), (l_xi2, 0, 1)  ]
      l_maps, l_absDets = Tet.intSc( 1, [l_xi0, l_xi1, l_xi2] )
      l_maps    = l_maps[0]    + l_maps[1]
      l_absDets = l_absDets[0] + l_absDets[1]

      l_nScs = len(l_maps)
      l_sumSc = [0, 0, 0, 0, 0]

      for l_sc in range(l_nScs):
        # assemble map dict
        l_td = { l_xi0: l_maps[l_sc][0],
                 l_xi1: l_maps[l_sc][1],
                 l_xi2: l_maps[l_sc][2] }

        l_sumSc[0] = l_sumSc[0] + sympy.integrate( sympy.sympify(1),                       *l_intRef ) * l_absDets[l_sc]
        l_sumSc[1] = l_sumSc[1] + sympy.integrate( l_xi0.subs( l_td, simultaneous=True ),  *l_intRef ) * l_absDets[l_sc]
        l_sumSc[2] = l_sumSc[2] + sympy.integrate( l_xi1.subs( l_td, simultaneous=True ),  *l_intRef ) * l_absDets[l_sc]
        l_sumSc[3] = l_sumSc[3] + sympy.integrate( l_xi2.subs( l_td, simultaneous=True  ), *l_intRef ) * l_absDets[l_sc]
        l_sumSc[4] = l_sumSc[4] + sympy.integrate( l_xi0.subs( l_td, simultaneous=True )*\
                                                   l_xi1.subs( l_td, simultaneous=True )*\
                                                   l_xi2.subs( l_td, simultaneous=True ),  *l_intRef ) * l_absDets[l_sc]


      l_sumEl = [0, 0, 0, 0, 0]
      l_sumEl[0] = sympy.integrate( sympy.sympify(1),  *l_intRef  )
      l_sumEl[1] = sympy.integrate( l_xi0,             *l_intRef  )
      l_sumEl[2] = sympy.integrate( l_xi1,             *l_intRef  )
      l_sumEl[3] = sympy.integrate( l_xi2,             *l_intRef  )
      l_sumEl[4] = sympy.integrate( l_xi0*l_xi1*l_xi2, *l_intRef  )

      self.assertEqual( l_sumSc, l_sumEl )

  ##
  # Tests integration intervals for DG sub-faces.
  ##
  def test_intSfDg(self):
    l_xi0  = sympy.symbols('xi_0')
    l_xi1  = sympy.symbols('xi_1')
    l_xi2  = sympy.symbols('xi_2')
    l_chi0 = sympy.symbols('chi_0')
    l_chi1 = sympy.symbols('chi_1')

    # second order
    l_subs, l_intSfDg = Tet.intSfDg( 1, [l_chi0, l_chi1], [l_xi0, l_xi1, l_xi2] )

    # check substitutions
    self.assertEqual( l_subs, ( ( ( l_xi0, l_chi1          ),
                                  ( l_xi1, l_chi0          ),
                                  ( l_xi2, 0               ) ),

                                ( ( l_xi0, l_chi0          ),
                                  ( l_xi1, 0               ),
                                  ( l_xi2, l_chi1          ) ),

                                ( ( l_xi0, 0               ),
                                  ( l_xi1, l_chi1          ),
                                  ( l_xi2, l_chi0          ) ),

                                ( ( l_xi0, 1-l_chi0-l_chi1 ),
                                  ( l_xi1, l_chi0          ),
                                  ( l_xi2, l_chi1          ) ) ) )

    for l_fa in l_intSfDg:
      self.assertEqual( l_fa, [ # first quad
                                [ ( l_chi0, Fra(0,3), Fra(1,3)-l_chi1 ),
                                  ( l_chi1, Fra(0,3), Fra(1,3)        ) ],
                                [ ( l_chi0, Fra(1,3)-l_chi1, Fra(1,3) ),
                                  ( l_chi1, Fra(0,3), Fra(1,3)        ) ],
                                # second quad
                                [ ( l_chi0, Fra(1,3), Fra(2,3)-l_chi1 ),
                                  ( l_chi1, Fra(0,3), Fra(1,3)        ) ],
                                [ ( l_chi0, Fra(2,3)-l_chi1, Fra(2,3) ),
                                  ( l_chi1, Fra(0,3), Fra(1,3)        ) ],
                                # third quad
                                [ ( l_chi0, Fra(2,3), Fra(3,3)-l_chi1 ),
                                  ( l_chi1, Fra(0,3), Fra(1,3)        ) ],
                                # fourth quad
                                [ ( l_chi0, Fra(0,3), Fra(1,3)-(l_chi1-Fra(1,3)) ),
                                  ( l_chi1, Fra(1,3), Fra(2,3)                   ) ],
                                [ ( l_chi0, Fra(1,3)-(l_chi1-Fra(1,3)), Fra(1,3) ),
                                  ( l_chi1, Fra(1,3), Fra(2,3)                   ) ],
                                # fifth quad
                                [ ( l_chi0, Fra(1,3), Fra(2,3)-(l_chi1-Fra(1,3)) ),
                                  ( l_chi1, Fra(1,3), Fra(2,3)                   ) ],
                                # sixth quad
                                [ ( l_chi0, Fra(0,3), Fra(1,3)-(l_chi1-Fra(2,3)) ),
                                  ( l_chi1, Fra(2,3), Fra(3,3)                   ) ],
                              ] )


  ##
  # Tests the derivation of sub-face types.
  ##
  def test_scTySf(self):
    # get sub-face types for order 2
    l_scTySfIn, l_scTySfSend = Tet.scTySf( 1 )

    # Remark: Clockwise storage if looking at triangle in normal direction

    # check inner types
    self.assertEqual( l_scTySfIn, [ [16, 14, 17, 15],
                                    [19, 16, 18, 11],
                                    [10, 12, 17, 13],
                                    [12,  9, 19, 14],
                                    [ 8, 13, 18, 15] ] )

    # check send first few types
    self.assertEqual( l_scTySfSend[0:5], [ [ 0,  1,  2, 11],
                                           [ 0, 13, 18, 15],
                                           [ 0,  1, 10, 11],
                                           [ 0, 13, 18, 15],
                                           [ 0,  1, 10,  3] ] )

    # check two more send types
    self.assertEqual( l_scTySfSend[12], [12,  1, 19, 14] ),
    self.assertEqual( l_scTySfSend[18], [ 8,  9,  2,  3] ),

  ##
  # Tests the sub-cell reordering for sub-cells at face of adjacent elements
  ##
  def test_scDgAd(self):
    # get reordering for deg 2
    l_scDgAd = Tet.scDgAd( 2 )

    # check reordering
    self.assertEqual( l_scDgAd[0],
                      [0,1,9,10,16,17,21,22,24,
                       2,3,11,12,18,19,23,
                       4,5,13,14,20,
                       6,7,15,
                       8])

    self.assertEqual( l_scDgAd[1],
                      [8,7,6,5,4,3,2,1,0,
                       15,14,13,12,11,10,9,
                       20,19,18,17,16,
                       23,22,21,
                       24])

    self.assertEqual( l_scDgAd[2],
                      [24,22,23,19,20,14,15,7,8,
                       21,17,18,12,13,5,6,
                       16,10,11,3,4,
                       9,1,2,
                       0])
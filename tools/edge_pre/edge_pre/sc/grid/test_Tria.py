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
# Unit tests for the tria sub-grid.
##
import unittest
from fractions import Fraction as Fra
from . import Tria
import sympy

# Reference grids:
#
# 2nd order
#
#    9
#    | *
#    |   *   14
# 15 |     *
#    |   7   *
#    |         *
#    7-----------8
#    | *         | *
#    |   *    2  |   *   13
# 16 |     *     |     *
#    |   8   *   |   6   *
#    |         * |         *
#    4-----------5-----------6
#    | *         | *         | *
#    |   *    0  |   *    1  |   *   12
# 17 |     *     |     *     |     *
#    |   3   *   |   4   *   |  5    *
#    |         * |         * |         *
#    0-----------1-----------2-----------3
#          9           10          11
#
# 3rd order
#
#    20
#    | *
#    |   *   34
# 35 |     *
#    |  21   *
#    |         *
#    18---------19
#    | *         | *
#    |   *   12  |   *   33
# 36 |     *     |     *
#    |  22   *   |  20   *
#    |         * |         *
#    15---------16----------17
#    | *         | *         | *
#    |   *    9  |   *   11  |   *   32
# 37 |     *     |     *     |     *
#    |  23   *   |   10  *   |  19   *
#    |         * |         * |         *
#    11---------12----------13----------14
#    | *         | *         | *         | *
#    |   *    4  |   *   6   |   *   8   |   *   31
# 38 |     *     |     *     |     *     |     *
#    |  24   *   |   5   *   |   7   *   |  18   *
#    |         * |         * |         * |         *
#    6-----------7-----------8-----------9----------10
#    | *         | *         | *         | *         | *
#    |   *    0  |   *    1  |   *    2  |   *    3  |   *   30
# 39 |     *     |     *     |     *     |     *     |     * 
#    |  13   *   |   14  *   |  15   *   |  16   *   |  17   *
#    |         * |         * |         * |         * |         *
#    0-----------1-----------2-----------3-----------4-----------5
#         25           26          27         28          29
#
class TestGridTria( unittest.TestCase ):
  ##
  # Tests generation of vertices.
  ##
  def test_svs(self):
    # FV
    l_svs = Tria.svs( 0 )
    self.assertEqual( l_svs, [ [0,0], [1,0], [0,1] ] )

    # 2nd order
    l_svs = Tria.svs( 1 )
    self.assertEqual( l_svs, [ [0, 0       ], [Fra(1,3),        0], [Fra(2,3),        0], [1,        0],
                               [0, Fra(1,3)], [Fra(1,3), Fra(1,3)], [Fra(2,3), Fra(1,3)],
                               [0, Fra(2,3)], [Fra(1,3), Fra(2,3)],
                               [0,        1] ] )

    # 3rd order
    l_svs = Tria.svs( 2 )
    self.assertEqual( l_svs, [ [0,               0], [Fra(1,5),        0], [Fra(2,5),        0],
                               [Fra(3,5),        0], [Fra(4,5),        0], [1,               0],

                               [0,        Fra(1,5)], [Fra(1,5), Fra(1,5)], [Fra(2,5), Fra(1,5)],
                               [Fra(3,5), Fra(1,5)], [Fra(4,5), Fra(1,5)],

                               [0,        Fra(2,5)], [Fra(1,5), Fra(2,5)], [Fra(2,5), Fra(2,5)],
                               [Fra(3,5), Fra(2,5)],

                               [0,        Fra(3,5)], [Fra(1,5), Fra(3,5)], [Fra(2,5), Fra(3,5)],

                               [0,        Fra(4,5)], [Fra(1,5), Fra(4,5)],

                               [0,               1] ] )

  ##
  # Tests scSv adjacency.
  ##
  def test_scSv(self):
    # second order sub-grid
    l_scSvIn, l_scSvSend, l_scSvRecv = Tria.scSv( 1 )

    # check inner sub-cells
    self.assertEqual( l_scSvIn, [ [1, 5, 4],
                                  [2, 6, 5],
                                  [5, 8, 7] ] )

    # check send sub-cells
    self.assertEqual( l_scSvSend, [ [0, 1, 4],
                                    [1, 2, 5],
                                    [2, 3, 6],
                                    [5, 6, 8],
                                    [7, 8, 9],
                                    [4, 5, 7] ] )

    self.assertEqual( l_scSvRecv, [ [0, 1, -1],
                                    [1, 2, -1],
                                    [2, 3, -1],
                                    [3, 6, -1],
                                    [6, 8, -1],
                                    [8, 9, -1],
                                    [9, 7, -1],
                                    [7, 4, -1],
                                    [4, 0, -1] ] )

    # third order sub-grid
    l_scSvIn, l_scSvSend, l_scSvRecv = Tria.scSv( 2 )

    # check inner sub-cells
    self.assertEqual( l_scSvIn, [ [ 1,  7,  6],
                                  [ 2,  8,  7],
                                  [ 3,  9,  8],
                                  [ 4, 10,  9],

                                  [ 7, 12, 11],
                                  [ 7,  8, 12],
                                  [ 8, 13, 12],
                                  [ 8,  9, 13],
                                  [ 9, 14, 13],
                                  [12, 16, 15],
                                  [12, 13, 16],
                                  [13, 17, 16],
                                  [16, 19, 18] ] )

    # check send sub-cells
    self.assertEqual( l_scSvSend, [ [ 0,  1,  6],
                                    [ 1,  2,  7],
                                    [ 2,  3,  8],
                                    [ 3,  4,  9],
                                    [ 4,  5, 10],
                                    [ 9, 10, 14],
                                    [13, 14, 17],
                                    [16, 17, 19],
                                    [18, 19, 20],
                                    [15, 16, 18],
                                    [11, 12, 15],
                                    [ 6,  7, 11] ] )

    # check receive sub-cells
    self.assertEqual( l_scSvRecv, [ [ 0,  1, -1],
                                    [ 1,  2, -1],
                                    [ 2,  3, -1],
                                    [ 3,  4, -1],
                                    [ 4,  5, -1],

                                    [ 5, 10, -1],
                                    [10, 14, -1],
                                    [14, 17, -1],
                                    [17, 19, -1],
                                    [19, 20, -1],

                                    [20, 18, -1],
                                    [18, 15, -1],
                                    [15, 11, -1],
                                    [11,  6, -1],
                                    [ 6,  0, -1] ] )

  ##
  # Tests scSfSv adjacency
  ##
  def test_scSfSv(self):
    # second order
    l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = Tria.scSfSv( 1 )

    # check inner sub-cells
    self.assertEqual( l_scSfSvIn, [ [ [ 1,  5], [ 5,  4], [ 4,  1] ],
                                    [ [ 2,  6], [ 6,  5], [ 5,  2] ],
                                    [ [ 5,  8], [ 8,  7], [ 7,  5] ] ] )

    # check send sub-cells
    self.assertEqual( l_scSfSvSend, [ [ [ 0,  1], [ 1,  4], [ 4,  0] ],
                                      [ [ 1,  2], [ 2,  5], [ 5,  1] ],
                                      [ [ 2,  3], [ 3,  6], [ 6,  2] ],
                                      [ [ 5,  6], [ 6,  8], [ 8,  5] ],
                                      [ [ 7,  8], [ 8,  9], [ 9,  7] ],
                                      [ [ 4,  5], [ 5,  7], [ 7,  4] ] ] )

    # check receive sub-cell
    self.assertEqual( l_scSfSvRecv, [ [ [ 0,  1], [-1, -1], [-1, -1] ],
                                      [ [ 1,  2], [-1, -1], [-1, -1] ],
                                      [ [ 2,  3], [-1, -1], [-1, -1] ],
                                      [ [ 3,  6], [-1, -1], [-1, -1] ],
                                      [ [ 6,  8], [-1, -1], [-1, -1] ],
                                      [ [ 8,  9], [-1, -1], [-1, -1] ],
                                      [ [ 9,  7], [-1, -1], [-1, -1] ],
                                      [ [ 7,  4], [-1, -1], [-1, -1] ],
                                      [ [ 4,  0], [-1, -1], [-1, -1] ] ] )


    # third order
    l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = Tria.scSfSv( 2 )

    # check inner sub-cells
    self.assertEqual( l_scSfSvIn, [ [ [ 1,  7], [ 7,  6], [ 6,  1] ],
                                    [ [ 2,  8], [ 8,  7], [ 7,  2] ],
                                    [ [ 3,  9], [ 9,  8], [ 8,  3] ],
                                    [ [ 4, 10], [10,  9], [ 9,  4] ],

                                    [ [ 7, 12], [12, 11], [11,  7] ],
                                    [ [ 7,  8], [ 8, 12], [12,  7] ],
                                    [ [ 8, 13], [13, 12], [12,  8] ],
                                    [ [ 8,  9], [ 9, 13], [13,  8] ],
                                    [ [ 9, 14], [14, 13], [13,  9] ],

                                    [ [12, 16], [16, 15], [15, 12] ],
                                    [ [12, 13], [13, 16], [16, 12] ],
                                    [ [13, 17], [17, 16], [16, 13] ],

                                    [ [16, 19], [19, 18], [18, 16] ] ] )

    # check send sub-cells
    self.assertEqual( l_scSfSvSend, [ [ [ 0,  1], [ 1,  6], [ 6,  0] ],
                                      [ [ 1,  2], [ 2,  7], [ 7,  1] ],
                                      [ [ 2,  3], [ 3,  8], [ 8,  2] ],
                                      [ [ 3,  4], [ 4,  9], [ 9,  3] ],
                                      [ [ 4,  5], [ 5, 10], [10,  4] ],

                                      [ [ 9, 10], [10, 14], [14,  9] ],
                                      [ [13, 14], [14, 17], [17, 13] ],
                                      [ [16, 17], [17, 19], [19, 16] ],
                                      [ [18, 19], [19, 20], [20, 18] ],

                                      [ [15, 16], [16, 18], [18, 15] ],
                                      [ [11, 12], [12, 15], [15, 11] ],
                                      [ [ 6,  7], [ 7, 11], [11,  6] ] ] )

    # check receive sub-cells
    self.assertEqual( l_scSfSvRecv, [ [ [ 0,  1], [-1, -1], [-1, -1] ],
                                      [ [ 1,  2], [-1, -1], [-1, -1] ],
                                      [ [ 2,  3], [-1, -1], [-1, -1] ],
                                      [ [ 3,  4], [-1, -1], [-1, -1] ],
                                      [ [ 4,  5], [-1, -1], [-1, -1] ],

                                      [ [ 5, 10], [-1, -1], [-1, -1] ],
                                      [ [10, 14], [-1, -1], [-1, -1] ],
                                      [ [14, 17], [-1, -1], [-1, -1] ],
                                      [ [17, 19], [-1, -1], [-1, -1] ],
                                      [ [19, 20], [-1, -1], [-1, -1] ],

                                      [ [20, 18], [-1, -1], [-1, -1] ],
                                      [ [18, 15], [-1, -1], [-1, -1] ],
                                      [ [15, 11], [-1, -1], [-1, -1] ],
                                      [ [11,  6], [-1, -1], [-1, -1] ],
                                      [ [ 6,  0], [-1, -1], [-1, -1] ] ] )



  ##
  # Tests scSfSc adjacency.
  ##
  def test_scSfSc(self):
    # second order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Tria.scSfSc( 1 )

    # check inner sub-cells
    self.assertEqual( l_scSfScIn, [ [ 4, 8, 3 ],
                                    [ 5, 6, 4 ],
                                    [ 6, 7, 8 ] ] )

    # check send sub-cells
    self.assertEqual( l_scSfScSend, [ [  9,  0, 17 ],
                                      [ 10,  1,  0 ],
                                      [ 11, 12,  1 ],

                                      [  1, 13,  2 ],
                                      [  2, 14, 15 ],

                                      [  0,  2, 16] ] )

    # check receive sub-cells
    self.assertEqual( l_scSfScRecv, [ [ 3, -1, -1 ],
                                      [ 4, -1, -1 ],
                                      [ 5, -1, -1 ],

                                      [ 5, -1, -1 ],
                                      [ 6, -1, -1 ],
                                      [ 7, -1, -1 ],

                                      [ 7, -1, -1 ],
                                      [ 8, -1, -1 ],
                                      [ 3, -1, -1 ] ] )

    # third order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Tria.scSfSc( 2 )

    # check inner sub-cells
    self.assertEqual( l_scSfScIn, [ [ 14, 24, 13 ],
                                    [ 15,  5, 14 ],
                                    [ 16,  7, 15 ],
                                    [ 17, 18, 16 ],

                                    [  5, 23, 24 ],
                                    [  1,  6,  4 ],
                                    [  7, 10,  5 ],
                                    [  2,  8,  6 ],
                                    [ 18, 19,  7 ],

                                    [ 10, 22, 23 ],
                                    [  6, 11,  9 ],
                                    [ 19, 20, 10 ],

                                    [ 20, 21, 22 ] ] )

    # check send sub-cells
    self.assertEqual( l_scSfScSend, [ [ 25,  0, 39 ],
                                      [ 26,  1,  0 ],
                                      [ 27,  2,  1 ],
                                      [ 28,  3,  2 ],
                                      [ 29, 30,  3 ],

                                      [  3, 31,  8 ],
                                      [  8, 32, 11 ],
                                      [ 11, 33, 12 ],
                                      [ 12, 34, 35 ],

                                      [  9, 12, 36 ],
                                      [  4,  9, 37 ],
                                      [  0,  4, 38 ] ] )

    # check receive sub-cells
    self.assertEqual( l_scSfScRecv, [ [ 13, -1, -1 ],
                                      [ 14, -1, -1 ],
                                      [ 15, -1, -1 ],
                                      [ 16, -1, -1 ],
                                      [ 17, -1, -1 ],

                                      [ 17, -1, -1 ],
                                      [ 18, -1, -1 ],
                                      [ 19, -1, -1 ],
                                      [ 20, -1, -1 ],
                                      [ 21, -1, -1 ],

                                      [ 21, -1, -1 ],
                                      [ 22, -1, -1 ],
                                      [ 23, -1, -1 ],
                                      [ 24, -1, -1 ],
                                      [ 13, -1, -1 ] ] )

  ##
  # Tests integration intervals for sub-cells.
  ##
  def test_intSc(self):
    l_xi0 = sympy.symbols('xi_0')
    l_xi1 = sympy.symbols('xi_1')

    # second order
    l_intIn, l_intSend, l_intSurf = Tria.intSc( 1, [l_xi0, l_xi1] )

    # inner intervals
    self.assertEqual( l_intIn, [ [ ( l_xi0, Fra(1,3)-l_xi1,            Fra(1,3) ),
                                   ( l_xi1, Fra(0,3),                  Fra(1,3) ) ],

                                 [ ( l_xi0, Fra(2,3)-l_xi1,            Fra(2,3) ),
                                   ( l_xi1, Fra(0,3),                  Fra(1,3) ) ],

                                 [ ( l_xi0, Fra(1,3)-(l_xi1-Fra(1,3)), Fra(1,3) ),
                                   ( l_xi1, Fra(1,3),                  Fra(2,3) ) ] ] )

    # send intervals
    self.assertEqual( l_intSend, [ [ ( l_xi0, Fra(0,3), Fra(1,3)-l_xi1            ),
                                     ( l_xi1, Fra(0,3), Fra(1,3)                  ) ],

                                   [ ( l_xi0, Fra(1,3), Fra(2,3)-l_xi1            ),
                                     ( l_xi1, Fra(0,3), Fra(1,3)                  ) ],

                                   [ ( l_xi0, Fra(2,3), Fra(3,3)-l_xi1            ),
                                     ( l_xi1, Fra(0,3), Fra(1,3)                  ) ],

                                   [ ( l_xi0, Fra(1,3), Fra(2,3)-(l_xi1-Fra(1,3)) ),
                                     ( l_xi1, Fra(1,3), Fra(2,3)                  ) ],

                                   [ ( l_xi0, Fra(0,3), Fra(1,3)-(l_xi1-Fra(2,3)) ),
                                     ( l_xi1, Fra(2,3), Fra(3,3)                  ) ],

                                   [ ( l_xi0, Fra(0,3), Fra(1,3)-(l_xi1-Fra(1,3)) ),
                                     ( l_xi1, Fra(1,3), Fra(2,3)                  ) ] ] )

    # DG intervals
    self.assertEqual( l_intSurf[0], [ [ ( l_xi0, Fra(0,3), Fra(1,3)-l_xi1 ),
                                        ( l_xi1, Fra(0,3), Fra(1,3)       ) ],

                                      [ ( l_xi0, Fra(1,3), Fra(2,3)-l_xi1 ),
                                        ( l_xi1, Fra(0,3), Fra(1,3)       ) ],

                                      [ ( l_xi0, Fra(2,3), Fra(3,3)-l_xi1 ),
                                        ( l_xi1, Fra(0,3), Fra(1,3)       ) ] ] )

    self.assertEqual( l_intSurf[1], [ [ ( l_xi0, Fra(2,3), Fra(3,3)-l_xi1 ),
                                        ( l_xi1, Fra(0,3), Fra(1,3)       ) ],

                                      [ ( l_xi0, Fra(1,3), Fra(2,3)-(l_xi1-Fra(1,3)) ),
                                        ( l_xi1, Fra(1,3), Fra(2,3)                  ) ],

                                      [ ( l_xi0, Fra(0,3), Fra(1,3)-(l_xi1-Fra(2,3)) ),
                                        ( l_xi1, Fra(2,3), Fra(3,3)                  ) ] ] )

    self.assertEqual( l_intSurf[2], [ [ ( l_xi0, Fra(0,3), Fra(1,3)-(l_xi1-Fra(2,3)) ),
                                        ( l_xi1, Fra(2,3), Fra(3,3)                  ) ],

                                      [ ( l_xi0, Fra(0,3), Fra(1,3)-(l_xi1-Fra(1,3)) ),
                                        ( l_xi1, Fra(1,3), Fra(2,3)                  ) ],

                                      [ ( l_xi0, Fra(0,3), Fra(1,3)-l_xi1 ),
                                        ( l_xi1, Fra(0,3), Fra(1,3)       ) ] ] )

  ##
  # Tests integration intervals for DG sub-faces.
  ##
  def test_intSfDg(self):
    l_xi0  = sympy.symbols('xi_0')
    l_xi1  = sympy.symbols('xi_1')
    l_chi0 = sympy.symbols('chi_0')

    # second order
    l_subs, l_intSfDg = Tria.intSfDg( 1, [l_chi0], [l_xi0, l_xi1] )

    # check substitutions
    self.assertEqual( l_subs, (
                                ( ( l_xi0, l_chi0   ),
                                  ( l_xi1, 0        ) ),
                                ( ( l_xi0, 1-l_chi0 ),
                                  ( l_xi1, l_chi0   ) ),
                                ( ( l_xi0, 0        ),
                                  ( l_xi1, l_chi0   ) )
                              ) )

    # check intervals
    self.assertEqual( l_intSfDg, [ [ # DG face0
                                     [ ( l_chi0, Fra(0,3), Fra(1,3) ) ],
                                     [ ( l_chi0, Fra(1,3), Fra(2,3) ) ],
                                     [ ( l_chi0, Fra(2,3), Fra(3,3) ) ]
                                   ],
                                   [ # DG face1
                                     [ ( l_chi0, Fra(0,3), Fra(1,3) ) ],
                                     [ ( l_chi0, Fra(1,3), Fra(2,3) ) ],
                                     [ ( l_chi0, Fra(2,3), Fra(3,3) ) ]
                                   ],
                                   [ # DG face3
                                     [ ( l_chi0, Fra(2,3), Fra(3,3) ) ],
                                     [ ( l_chi0, Fra(1,3), Fra(2,3) ) ],
                                     [ ( l_chi0, Fra(0,3), Fra(1,3) ) ]
                                   ] ] )

    # third order
    l_subs, l_intSfDg = Tria.intSfDg( 2, [l_chi0], [l_xi0, l_xi1] )

    # check substitutions
    self.assertEqual( l_subs, (
                                ( ( l_xi0, l_chi0   ),
                                  ( l_xi1, 0        ) ),
                                ( ( l_xi0, 1-l_chi0 ),
                                  ( l_xi1, l_chi0   ) ),
                                ( ( l_xi0, 0        ),
                                  ( l_xi1, l_chi0   ) )
                              ) )

    self.assertEqual( l_intSfDg, [ [ # DG face0
                                     [ ( l_chi0, Fra(0,5), Fra(1,5) ) ],
                                     [ ( l_chi0, Fra(1,5), Fra(2,5) ) ],
                                     [ ( l_chi0, Fra(2,5), Fra(3,5) ) ],
                                     [ ( l_chi0, Fra(3,5), Fra(4,5) ) ],
                                     [ ( l_chi0, Fra(4,5), Fra(5,5) ) ]
                                   ],
                                   [ # DG face1
                                     [ ( l_chi0, Fra(0,5), Fra(1,5) ) ],
                                     [ ( l_chi0, Fra(1,5), Fra(2,5) ) ],
                                     [ ( l_chi0, Fra(2,5), Fra(3,5) ) ],
                                     [ ( l_chi0, Fra(3,5), Fra(4,5) ) ],
                                     [ ( l_chi0, Fra(4,5), Fra(5,5) ) ]
                                   ],
                                   [ # DG face3
                                     [ ( l_chi0, Fra(4,5), Fra(5,5) ) ],
                                     [ ( l_chi0, Fra(3,5), Fra(4,5) ) ],
                                     [ ( l_chi0, Fra(2,5), Fra(3,5) ) ],
                                     [ ( l_chi0, Fra(1,5), Fra(2,5) ) ],
                                     [ ( l_chi0, Fra(0,5), Fra(1,5) ) ]
                                   ] ] )

  ##
  # Tests the derivation of sub-face types.
  #
  # Sub-faces at the DG-surface have types 0-5.
  # Sub-faces not at the DG-surface have types 6-11.
  ##
  def test_scTySf(self):
    # get sub-face types for order 2
    l_scTySfIn, l_scTySfSend = Tria.scTySf( 1 )

    # check inner types
    self.assertEqual( l_scTySfIn, [ [11, 9, 10],
                                    [11, 9, 10],
                                    [11, 9, 10] ] )

    # check send types
    self.assertEqual( l_scTySfSend, [ [0, 7, 2],
                                      [0, 7, 8],
                                      [0, 1, 8],
                                      [6, 1, 8],
                                      [6, 1, 2],
                                      [6, 7, 2] ] )

    # get sub-face types for order 3
    l_scTySfIn, l_scTySfSend = Tria.scTySf( 2 )

    # check inner types
    self.assertEqual( l_scTySfIn, [ [11, 9, 10],
                                    [11, 9, 10],
                                    [11, 9, 10],
                                    [11, 9, 10],

                                    [11, 9, 10],
                                    [ 6, 7,  8],
                                    [11, 9, 10],
                                    [ 6, 7,  8],
                                    [11, 9, 10],

                                    [11, 9, 10],
                                    [ 6, 7,  8],
                                    [11, 9, 10],

                                    [11, 9, 10] ] )

    # check send types
    self.assertEqual( l_scTySfSend, [ [0, 7, 2],
                                      [0, 7, 8],
                                      [0, 7, 8],
                                      [0, 7, 8],
                                      [0, 1, 8],
                                      [6, 1, 8],
                                      [6, 1, 8],
                                      [6, 1, 8],
                                      [6, 1, 2],
                                      [6, 7, 2],
                                      [6, 7, 2],
                                      [6, 7, 2] ] )

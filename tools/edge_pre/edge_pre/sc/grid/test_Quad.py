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
# Unit tests for the quad sub-grid.
##
import unittest
from fractions import Fraction as Fra
from . import Quad

# Reference grids:
#
# 2nd order
#
#        17      16       15
#    12-----13------14------15
#    |       |       |       |
# 18 |   7   |   6   |   5   | 14
#    |       |       |       |
#    8-------9------10------11
#    |       |       |       |
# 19 |   8   |   0   |   4   | 13
#    |       |       |       |
#    4-------5-------6-------7
#    |       |       |       |
# 20 |   1   |   2   |   3   | 12
#    |       |       |       |
#    0-------1-------2-------3
#        9      10      11
#
#
# 3rd order
#       40      39      38      37       36
#    30-----31------32------33------34------35
#    |       |       |       |       |       |
# 41 |  21   |  20   |  19   |  18   |   17  | 35
#    |       |       |       |       |       |
#    24-----25------26------27------28------29
#    |       |       |       |       |       |
# 42 |  22   |   6   |   7   |   8   |   16  | 34
#    |       |       |       |       |       |
#    18-----19------20------21------22------23
#    |       |       |       |       |       |
# 43 |  23   |   3   |   4   |   5   |   15  | 33
#    |       |       |       |       |       |
#    12-----13------14------15------16------17
#    |       |       |       |       |       |
# 44 |  24   |   0   |   1   |   2   |   14  | 32
#    |       |       |       |       |       |
#    6-------7-------8-------9------10------11
#    |       |       |       |       |       |
# 45 |   9   |  10   |  11   |   12  |   13  | 31
#    |       |       |       |       |       |
#    0-------1-------2-------3-------4-------5
#       26      27      28       29      30
#
class TestGridQuad( unittest.TestCase ):
  ##
  # Tests generation of vertices.
  ##
  def test_svs(self):
    # FV
    l_svs = Quad.svs( 0 )
    self.assertEqual( l_svs, [ [0,0], [1,0],
                               [0,1], [1,1] ] )

    # 2nd order
    l_svs = Quad.svs( 1 )
    self.assertEqual( l_svs, [ [0, 0       ], [Fra(1,3),        0], [Fra(2,3),        0], [1,        0],
                               [0, Fra(1,3)], [Fra(1,3), Fra(1,3)], [Fra(2,3), Fra(1,3)], [1, Fra(1,3)],
                               [0, Fra(2,3)], [Fra(1,3), Fra(2,3)], [Fra(2,3), Fra(2,3)], [1, Fra(2,3)],
                               [0,        1], [Fra(1,3),        1], [Fra(2,3),        1], [1,        1] ] )

    # 3rd order
    l_svs = Quad.svs( 2 )
    self.assertEqual( l_svs, [ [0,               0], [Fra(1,5),        0], [Fra(2,5),        0],
                               [Fra(3,5),        0], [Fra(4,5),        0], [1,               0],

                               [0,        Fra(1,5)], [Fra(1,5), Fra(1,5)], [Fra(2,5), Fra(1,5)],
                               [Fra(3,5), Fra(1,5)], [Fra(4,5), Fra(1,5)], [1,        Fra(1,5)],

                               [0,        Fra(2,5)], [Fra(1,5), Fra(2,5)], [Fra(2,5), Fra(2,5)],
                               [Fra(3,5), Fra(2,5)], [Fra(4,5), Fra(2,5)], [1,        Fra(2,5)],

                               [0,        Fra(3,5)], [Fra(1,5), Fra(3,5)], [Fra(2,5), Fra(3,5)],
                               [Fra(3,5), Fra(3,5)], [Fra(4,5), Fra(3,5)], [1,        Fra(3,5)],

                               [0,        Fra(4,5)], [Fra(1,5), Fra(4,5)], [Fra(2,5), Fra(4,5)],
                               [Fra(3,5), Fra(4,5)], [Fra(4,5), Fra(4,5)], [1,        Fra(4,5)],

                               [0,               1], [Fra(1,5),        1], [Fra(2,5),        1],
                               [Fra(3,5),        1], [Fra(4,5),        1], [1,               1] ] )

  ##
  # Tests scSv adjacency.
  ##
  def test_scSv(self):
    l_scSvIn, l_scSvSend, l_scSvRecv = Quad.scSv( 1 )

    # check inner sub-cell
    self.assertEqual( l_scSvIn, [ [5, 6, 10, 9] ] )

    # check send sub-cells
    self.assertEqual( l_scSvSend, [ [ 0,  1,  5,  4],
                                    [ 1,  2,  6,  5],
                                    [ 2,  3,  7,  6],
                                    [ 6,  7, 11, 10],
                                    [10, 11, 15, 14],
                                    [ 9, 10, 14, 13],
                                    [ 8,  9, 13, 12],
                                    [ 4,  5,  9,  8] ] )

    # check receive sub-cells
    self.assertEqual( l_scSvRecv, [ [ 0,  1, -1, -1],
                                    [ 1,  2, -1, -1],
                                    [ 2,  3, -1, -1],
                                    [ 3,  7, -1, -1],
                                    [ 7, 11, -1, -1],
                                    [11, 15, -1, -1],
                                    [15, 14, -1, -1],
                                    [14, 13, -1, -1],
                                    [13, 12, -1, -1],
                                    [12,  8, -1, -1],
                                    [ 8,  4, -1, -1],
                                    [ 4,  0, -1, -1] ] )

    l_scSvIn, l_scSvSend, l_scSvRecv = Quad.scSv( 2 )

    # check inner sub-cells
    self.assertEqual( l_scSvIn, [ [ 7,  8, 14, 13],
                                  [ 8,  9, 15, 14],
                                  [ 9, 10, 16, 15],

                                  [13, 14, 20, 19],
                                  [14, 15, 21, 20],
                                  [15, 16, 22, 21],

                                  [19, 20, 26, 25],
                                  [20, 21, 27, 26],
                                  [21, 22, 28, 27] ] )

    # check send sub-cells
    self.assertEqual( l_scSvSend, [ [ 0,  1,  7,  6],
                                    [ 1,  2,  8,  7],
                                    [ 2,  3,  9,  8],
                                    [ 3,  4, 10,  9],
                                    [ 4,  5, 11, 10],

                                    [10, 11, 17, 16],
                                    [16, 17, 23, 22],
                                    [22, 23, 29, 28],
                                    [28, 29, 35, 34],

                                    [27, 28, 34, 33],
                                    [26, 27, 33, 32],
                                    [25, 26, 32, 31],
                                    [24, 25, 31, 30],

                                    [18, 19, 25, 24],
                                    [12, 13, 19, 18],
                                    [ 6,  7, 13, 12] ] )

    # check receive sub-cells
    self.assertEqual( l_scSvRecv, [ [ 0,  1, -1, -1],
                                    [ 1,  2, -1, -1],
                                    [ 2,  3, -1, -1],
                                    [ 3,  4, -1, -1],
                                    [ 4,  5, -1, -1],

                                    [ 5, 11, -1, -1],
                                    [11, 17, -1, -1],
                                    [17, 23, -1, -1],
                                    [23, 29, -1, -1],
                                    [29, 35, -1, -1],

                                    [35, 34, -1, -1],
                                    [34, 33, -1, -1],
                                    [33, 32, -1, -1],
                                    [32, 31, -1, -1],
                                    [31, 30, -1, -1],

                                    [30, 24, -1, -1],
                                    [24, 18, -1, -1],
                                    [18, 12, -1, -1],
                                    [12,  6, -1, -1],
                                    [ 6,  0, -1, -1] ] )

  ##
  # Tests scSfSv adjacency
  ##
  def test_scSfSv(self):
    # second order
    l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = Quad.scSfSv( 1 )

    # check inner sub-cell
    self.assertEqual( l_scSfSvIn, [ [ [5,6], [6,10], [10,9], [9,5] ] ] )

    # check send sub-cell
    self.assertEqual( l_scSfSvSend, [ [ [ 0, 1], [ 1, 5], [ 5, 4], [ 4, 0] ],
                                      [ [ 1, 2], [ 2, 6], [ 6, 5], [ 5, 1] ],
                                      [ [ 2, 3], [ 3, 7], [ 7, 6], [ 6, 2] ],
                                      [ [ 6, 7], [ 7,11], [11,10], [10, 6] ],
                                      [ [10,11], [11,15], [15,14], [14,10] ],
                                      [ [ 9,10], [10,14], [14,13], [13, 9] ],
                                      [ [ 8, 9], [ 9,13], [13,12], [12, 8] ],
                                      [ [ 4, 5], [ 5, 9], [ 9, 8], [ 8, 4] ] ] )

    # check send sub-cell
    self.assertEqual( l_scSfSvRecv, [ [ [ 0, 1], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [ 1, 2], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [ 2, 3], [-1,-1], [-1,-1], [-1,-1] ],

                                      [ [ 3, 7], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [ 7,11], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [11,15], [-1,-1], [-1,-1], [-1,-1] ],

                                      [ [15,14], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [14,13], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [13,12], [-1,-1], [-1,-1], [-1,-1] ],

                                      [ [12, 8], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [ 8, 4], [-1,-1], [-1,-1], [-1,-1] ],
                                      [ [ 4, 0], [-1,-1], [-1,-1], [-1,-1] ] ] )

  ##
  # Tests scSfSc adjacency.
  ##
  def test_scSfSc(self):
    # second order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Quad.scSfSc( 1 )

    # check inner sub-cell
    self.assertEqual( l_scSfScIn, [ [ 2, 4, 6, 8 ] ] )

    # check send sub-cells
    self.assertEqual( l_scSfScSend, [ [ 9, 2, 8,20 ],
                                      [10, 3, 0, 1 ],
                                      [11,12, 4, 2 ],

                                      [ 3,13, 5, 0 ],
                                      [ 4,14,15, 6 ],

                                      [ 0, 5,16, 7 ],
                                      [ 8, 6,17,18 ],
                                      [ 1, 0, 7,19 ] ])

    # check receive sub-cells
    self.assertEqual( l_scSfScRecv, [ [ 1,-1,-1,-1 ],
                                      [ 2,-1,-1,-1 ],
                                      [ 3,-1,-1,-1 ],

                                      [ 3,-1,-1,-1 ],
                                      [ 4,-1,-1,-1 ],
                                      [ 5,-1,-1,-1 ],

                                      [ 5,-1,-1,-1 ],
                                      [ 6,-1,-1,-1 ],
                                      [ 7,-1,-1,-1 ],

                                      [ 7,-1,-1,-1 ],
                                      [ 8,-1,-1,-1 ],
                                      [ 1,-1,-1,-1 ] ])

  ##
  # Tests integration intervals for sub-cells.
  ##
  def test_intSc(self):
    # second order
    l_intIn, l_intSend, l_intSurf = Quad.intSc( 1, ['xi0', 'xi1'] )

    self.assertEqual( l_intIn, [ [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                   ( 'xi1', Fra(1,3), Fra(2,3) ) ] ] )

    self.assertEqual( l_intSend, [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ) ],

                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ) ] ] )

    self.assertEqual( l_intSurf[0], [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                        ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                      [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                        ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                      [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                        ( 'xi1', Fra(0,3), Fra(1,3) ) ] ] )

    self.assertEqual( l_intSurf[1], [ [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                        ( 'xi1', Fra(0,3), Fra(1,3) ) ],

                                      [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                        ( 'xi1', Fra(1,3), Fra(2,3) ) ],

                                      [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                        ( 'xi1', Fra(2,3), Fra(3,3) ) ] ] )

    self.assertEqual( l_intSurf[2], [ [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                        ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                      [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                        ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                      [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                        ( 'xi1', Fra(2,3), Fra(3,3) ) ] ] )

    self.assertEqual( l_intSurf[3], [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                        ( 'xi1', Fra(2,3), Fra(3,3) ) ],

                                      [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                        ( 'xi1', Fra(1,3), Fra(2,3) ) ],

                                      [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                        ( 'xi1', Fra(0,3), Fra(1,3) ) ] ] )

  ##
  # Tests integration intervals for DG sub-faces.
  ##
  def test_intSfDg(self):
    l_subs, l_intSfDg = Quad.intSfDg( 1, ['chi0'], ['xi0', 'xi1'] )

    # check substitutions
    self.assertEqual( l_subs, (
                                ( ( 'xi0', 'chi0' ),
                                  ( 'xi1', 0      ) ),
                                ( ( 'xi0', 1      ),
                                  ( 'xi1', 'chi0' ) ),
                                ( ( 'xi0', 'chi0' ),
                                  ( 'xi1', 1      ) ),
                                ( ( 'xi0', 0      ),
                                  ( 'xi1', 'chi0' ) )
                              ) )

    # check intervals
    self.assertEqual( l_intSfDg, [ [ # DG face0
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ) ],
                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ) ],
                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ) ]
                                   ],
                                   [ # DG face1
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ) ],
                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ) ],
                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ) ]
                                   ],
                                   [ # DG face2
                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ) ],
                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ) ],
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ) ]
                                   ],
                                   [ # DG face3
                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ) ],
                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ) ],
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ) ]
                                   ] ] )

  ##
  # Tests the derivation of sub-face types.
  #
  # Sub-faces at the DG-surface have types 0-3.
  # Sub-faces not at the DG-surface have types 8-11.
  ##
  def test_scTySf(self):
    # get sub-face types for order 2
    l_scTySfIn, l_scTySfSend = Quad.scTySf( 1 )

    # check inner types
    self.assertEqual( l_scTySfIn, [ [8, 9, 10, 11] ] )

    # check send types
    self.assertEqual( l_scTySfSend, [ [ 0,  9, 10,  3], # bottom-left corner
                                      [ 0,  9, 10, 11],
                                      [ 0,  1, 10, 11], # bottom-right corner
                                      [ 8,  1, 10, 11],
                                      [ 8,  1,  2, 11], # top-right corner
                                      [ 8,  9,  2, 11],
                                      [ 8,  9,  2,  3], # top-left corner
                                      [ 8,  9, 10,  3] ] )

    # get sub-face types for order 3
    l_scTySfIn, l_scTySfSend = Quad.scTySf( 2 )

    # check inner types
    self.assertEqual( l_scTySfIn, [ [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11],
                                    [8, 9, 10, 11] ] )

    # check send types
    self.assertEqual( l_scTySfSend, [ [0, 9, 10,  3], # bottom-left corner
                                      [0, 9, 10, 11],
                                      [0, 9, 10, 11],
                                      [0, 9, 10, 11],
                                      [0, 1, 10, 11], # bottom-right corner
                                      [8, 1, 10, 11],
                                      [8, 1, 10, 11],
                                      [8, 1, 10, 11],
                                      [8, 1,  2, 11], # top-right corner
                                      [8, 9,  2, 11],
                                      [8, 9,  2, 11],
                                      [8, 9,  2, 11],
                                      [8, 9,  2,  3], # top-left corner
                                      [8, 9, 10,  3],
                                      [8, 9, 10,  3],
                                      [8, 9, 10,  3] ] )

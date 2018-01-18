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
# Unit tests for the hex sub-grid.
##
import unittest
from fractions import Fraction as Fra
from . import Hex

# Reference grids:
#
# Elements are sorted by:
#   1) inner, send, receive
#   2) face if send or receive (bottom, front, right, back, left, top)
#   3) dimension x, y, z
#
# 2nd order
#                        _____________________
#                       /      /      /      /|
#                      /______/______/______/ |
#          .........  /      /      /      /| |
#            layer . /______/______/______/ | |
#            3     ./      /      /      /| |/|
#  slice 3 ......../______/______/______/ | / |
#            layer |      |      |      | |/| |
#            2     |      |      |      | / |/|
#  slice 2 ........|______|______|______|/| / |
#            layer |      |      |      | |/| |
#            1     |      |      |      | / |/
#  slice 1 ........|______|______|______|/| /
#            layer |      |      |      | |/ 
#            0     |      |      |      | /
#  slice 0 ........|______|______|______|/
#            layer .
#            -1    .
#          .........
#
#                     first | size
#   inner elements:   0     | 1
#   send elements:    1     | 26
#   receive elements: 27    | 54
#
#   slice 0 (first/second element id is below/above):
#
#           --/54     --/55     --/56
#        12-------13--------14--------15
#        |         |         |         |
#  --/65 |  33/07  |  34/08  |  35/09  | --/47
#        |         |         |         |
#        8---------9--------10--------11
#        |         |         |         |
#  --/64 |  30/04  |  31/05  |  32/06  | --/46
#        |         |         |         |
#        4---------5---------6---------7
#        |         |         |         |
#  --/63 |  27/01  |  28/02  |  29/03  | --/45
#        |         |         |         |
#        0---------1---------2---------3
#           --/36     --/37     --/38
#
#
#   slice 1 (first/second element id is below/above):
#
#           54/57     55/58     56/59
#        28-------29--------30--------31
#        |         |         |         |
#  65/68 |  07/20  |  08/21  |  09/17  | 47/50
#        |         |         |         |
#        24-------25--------26--------27
#        |         |         |         |
#  64/67 |  04/24  |  05/00  |  06/16  | 46/49
#        |         |         |         |
#        20-------21--------22--------23
#        |         |         |         |
#  63/66 |  01/10  |  02/11  |  03/12  | 45/48
#        |         |         |         |
#        16-------17--------18--------19
#           36/39     37/40     38/41
#
#
#   slice 2 (first/second element id is below/above):
#
#           57/60     58/61     59/62
#        44-------45--------46--------47
#        |         |         |         |
#  68/71 |  20/22  |  21/23  |  17/19  | 50/53
#        |         |         |         |
#        40-------41--------42--------43
#        |         |         |         |
#  67/70 |  24/25  |  00/27  |  16/18  | 49/52
#        |         |         |         |
#        36-------37--------38--------39
#        |         |         |         |
#  66/69 |  10/13  |  11/14  |  12/15  | 48/51
#        |         |         |         |
#        32-------33--------34--------35
#           39/42     40/43     41/44
#
#
#   slice 3 (first/second element id is below/above):
#
#           60/--     61/--      62/--
#        60-------61--------62--------63
#        |         |         |         |
#  71/-- |  22/78  |  23/79  |   19/80 | 53/--
#        |         |         |         |
#        56-------57--------58--------59
#        |         |         |         |
#  70/-- |  25/75  |  26/76  |  18/77  | 52/--
#        |         |         |         |
#        52-------53--------54--------55
#        |         |         |         |
#  69/-- |  13/72  |  14/73  |  15/74  | 51/--
#        |         |         |         |
#        48-------49--------50--------51
#           42/--     43/--      44/--
#
class TestGridHex( unittest.TestCase ):
  ##
  # Tests generation of vertices.
  ##
  def test_svs( self ):
    # 1st order
    l_svs = Hex.svs( 0 )
    self.assertEqual( l_svs, [ [0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1], [1,0,1], [0,1,1], [1,1,1] ] )

    # 2nd order
    l_svs = Hex.svs( 1 )
    self.assertEqual( l_svs, [ [0,        0,        0], [Fra(1,3),        0,        0], [Fra(2,3),        0,        0], [Fra(3,3),        0,        0],
                               [0, Fra(1,3),        0], [Fra(1,3), Fra(1,3),        0], [Fra(2,3), Fra(1,3),        0], [Fra(3,3), Fra(1,3),        0],
                               [0, Fra(2,3),        0], [Fra(1,3), Fra(2,3),        0], [Fra(2,3), Fra(2,3),        0], [Fra(3,3), Fra(2,3),        0],
                               [0, Fra(3,3),        0], [Fra(1,3), Fra(3,3),        0], [Fra(2,3), Fra(3,3),        0], [Fra(3,3), Fra(3,3),        0],

                               [0,        0, Fra(1,3)], [Fra(1,3),        0, Fra(1,3)], [Fra(2,3),        0, Fra(1,3)], [Fra(3,3),        0, Fra(1,3)],
                               [0, Fra(1,3), Fra(1,3)], [Fra(1,3), Fra(1,3), Fra(1,3)], [Fra(2,3), Fra(1,3), Fra(1,3)], [Fra(3,3), Fra(1,3), Fra(1,3)],
                               [0, Fra(2,3), Fra(1,3)], [Fra(1,3), Fra(2,3), Fra(1,3)], [Fra(2,3), Fra(2,3), Fra(1,3)], [Fra(3,3), Fra(2,3), Fra(1,3)],
                               [0, Fra(3,3), Fra(1,3)], [Fra(1,3), Fra(3,3), Fra(1,3)], [Fra(2,3), Fra(3,3), Fra(1,3)], [Fra(3,3), Fra(3,3), Fra(1,3)],

                               [0,        0, Fra(2,3)], [Fra(1,3),        0, Fra(2,3)], [Fra(2,3),        0, Fra(2,3)], [Fra(3,3),        0, Fra(2,3)],
                               [0, Fra(1,3), Fra(2,3)], [Fra(1,3), Fra(1,3), Fra(2,3)], [Fra(2,3), Fra(1,3), Fra(2,3)], [Fra(3,3), Fra(1,3), Fra(2,3)],
                               [0, Fra(2,3), Fra(2,3)], [Fra(1,3), Fra(2,3), Fra(2,3)], [Fra(2,3), Fra(2,3), Fra(2,3)], [Fra(3,3), Fra(2,3), Fra(2,3)],
                               [0, Fra(3,3), Fra(2,3)], [Fra(1,3), Fra(3,3), Fra(2,3)], [Fra(2,3), Fra(3,3), Fra(2,3)], [Fra(3,3), Fra(3,3), Fra(2,3)],

                               [0,        0, Fra(3,3)], [Fra(1,3),        0, Fra(3,3)], [Fra(2,3),        0, Fra(3,3)], [Fra(3,3),        0, Fra(3,3)],
                               [0, Fra(1,3), Fra(3,3)], [Fra(1,3), Fra(1,3), Fra(3,3)], [Fra(2,3), Fra(1,3), Fra(3,3)], [Fra(3,3), Fra(1,3), Fra(3,3)],
                               [0, Fra(2,3), Fra(3,3)], [Fra(1,3), Fra(2,3), Fra(3,3)], [Fra(2,3), Fra(2,3), Fra(3,3)], [Fra(3,3), Fra(2,3), Fra(3,3)],
                               [0, Fra(3,3), Fra(3,3)], [Fra(1,3), Fra(3,3), Fra(3,3)], [Fra(2,3), Fra(3,3), Fra(3,3)], [Fra(3,3), Fra(3,3), Fra(3,3)],
                             ] )

  ##
  # Tests scSv adjacency.
  ##
  def test_scSv( self ):
    # first order
    l_scSvIn, l_scSvSend, l_scSvRecv = Hex.scSv( 0 )

    self.assertEqual( l_scSvIn, [] )

    self.assertEqual( l_scSvSend, [ [ 0, 1, 2, 3,
                                      4, 5, 6, 7 ]
                                  ] )

    self.assertEqual( l_scSvRecv, [ [ -1, -1, -1, -1,
                                       0,  1,  2,  3 ],
                                    [ -1, -1,  0,  1,
                                      -1, -1,  4,  5 ],
                                    [  1, -1,  3, -1,
                                       5, -1,  7, -1 ],
                                    [  2,  3, -1, -1,
                                       6,  7, -1, -1 ],
                                    [ -1,  0, -1,  2,
                                      -1,  4, -1,  6 ],
                                    [  4,  5,  6,  7,
                                      -1, -1, -1, -1 ]
                                  ] )

    # 2nd order
    l_scSvIn, l_scSvSend, l_scSvRecv = Hex.scSv( 1 )
  
    self.assertEqual( l_scSvIn, [ [21, 22, 25, 26,
                                   37, 38, 41, 42]
                                ] )

    self.assertEqual( l_scSvSend, [ # bottom
                                    [  0,  1,  4,  5,
                                      16, 17, 20, 21 ],
                                    [  1,  2,  5,  6,
                                      17, 18, 21, 22 ],
                                    [  2,  3,  6,  7,
                                      18, 19, 22, 23 ],
                                    [  4,  5,  8,  9,
                                      20, 21, 24, 25 ],
                                    [  5,  6,  9, 10,
                                      21, 22, 25, 26 ],
                                    [  6,  7, 10, 11,
                                      22, 23, 26, 27 ],
                                    [  8,  9, 12, 13,
                                      24, 25, 28, 29 ],
                                    [  9, 10, 13, 14,
                                      25, 26, 29, 30 ],
                                    [ 10, 11, 14, 15,
                                      26, 27, 30, 31 ],
                                    # front
                                    [ 16, 17, 20, 21,
                                      32, 33, 36, 37 ],
                                    [ 17, 18, 21, 22,
                                      33, 34, 37, 38 ],
                                    [ 18, 19, 22, 23,
                                      34, 35, 38, 39 ],
                                    [ 32, 33, 36, 37,
                                      48, 49, 52, 53 ],
                                    [ 33, 34, 37, 38,
                                      49, 50, 53, 54 ],
                                    [ 34, 35, 38, 39,
                                      50, 51, 54, 55 ],
                                    # right
                                    [ 22, 23, 26, 27,
                                      38, 39, 42, 43 ],
                                    [ 26, 27, 30, 31,
                                      42, 43, 46, 47 ],
                                    [ 38, 39, 42, 43,
                                      54, 55, 58, 59 ],
                                    [ 42, 43, 46, 47,
                                      58, 59, 62, 63 ],
                                    # back
                                    [ 24, 25, 28, 29,
                                      40, 41, 44, 45 ],
                                    [ 25, 26, 29, 30,
                                      41, 42, 45, 46 ],
                                    [ 40, 41, 44, 45,
                                      56, 57, 60, 61 ],
                                    [ 41, 42, 45, 46,
                                      57, 58, 61, 62 ],
                                    # left
                                    [ 20, 21, 24, 25,
                                      36, 37, 40, 41 ],
                                    [ 36, 37, 40, 41,
                                      52, 53, 56, 57 ],
                                   # top
                                    [ 37, 38, 41, 42,
                                      53, 54, 57, 58 ]
                                  ] )

    self.assertEqual( l_scSvRecv, [ # bottom
                                    [ -1, -1, -1, -1,
                                       0,  1,  4,  5 ],
                                    [ -1, -1, -1, -1,
                                       1,  2,  5,  6 ],
                                    [ -1, -1, -1, -1,
                                       2,  3,  6,  7 ],
                                    [ -1, -1, -1, -1,
                                       4,  5,  8,  9 ],
                                    [ -1, -1, -1, -1,
                                       5,  6,  9, 10 ],
                                    [ -1, -1, -1, -1,
                                       6,  7, 10, 11 ],
                                    [ -1, -1, -1, -1,
                                       8,  9, 12, 13 ],
                                    [ -1, -1, -1, -1,
                                       9, 10, 13, 14 ],
                                    [ -1, -1, -1, -1,
                                      10, 11, 14, 15 ],
                                    # front
                                    [ -1, -1,  0,  1,
                                      -1, -1, 16, 17 ],
                                    [ -1, -1,  1,  2,
                                      -1, -1, 17, 18 ],
                                    [ -1, -1,  2,  3,
                                      -1, -1, 18, 19 ],
                                    [ -1, -1, 16, 17,
                                      -1, -1 , 32, 33 ],
                                    [ -1, -1, 17, 18,
                                      -1, -1, 33, 34 ],
                                    [ -1, -1, 18, 19,
                                      -1, -1, 34, 35 ],
                                    [ -1, -1, 32, 33,
                                      -1, -1, 48, 49 ],
                                    [ -1, -1, 33, 34,
                                      -1, -1, 49, 50 ],
                                    [ -1, -1, 34, 35,
                                      -1, -1, 50, 51 ],
                                    # right
                                    [  3, -1,  7, -1,
                                      19, -1, 23, -1 ],
                                    [  7, -1, 11, -1,
                                      23, -1, 27, -1 ],
                                    [ 11, -1, 15, -1,
                                      27, -1, 31, -1 ],
                                    [ 19, -1, 23, -1,
                                      35, -1, 39, -1 ],
                                    [ 23, -1, 27, -1,
                                      39, -1, 43, -1 ],
                                    [ 27, -1, 31, -1,
                                      43, -1, 47, -1 ],
                                    [ 35, -1, 39, -1,
                                      51, -1, 55, -1 ],
                                    [ 39, -1, 43, -1,
                                      55, -1, 59, -1 ],
                                    [ 43, -1, 47, -1,
                                      59, -1, 63, -1 ],
                                    # back
                                    [ 12, 13, -1, -1,
                                      28, 29, -1, -1 ],
                                    [ 13, 14, -1, -1,
                                      29, 30, -1, -1 ],
                                    [ 14, 15, -1, -1,
                                      30, 31, -1, -1 ],
                                    [ 28, 29, -1, -1,
                                      44, 45, -1, -1 ],
                                    [ 29, 30, -1, -1,
                                      45, 46, -1, -1 ],
                                    [ 30, 31, -1, -1,
                                      46, 47, -1, -1 ],
                                    [ 44, 45, -1, -1,
                                      60, 61, -1, -1 ],
                                    [ 45, 46, -1, -1,
                                      61, 62, -1, -1 ],
                                    [ 46, 47, -1, -1,
                                      62, 63, -1, -1 ],
                                    # left
                                    [ -1,  0, -1,  4,
                                      -1, 16, -1, 20 ],
                                    [ -1,  4, -1,  8,
                                      -1, 20, -1, 24 ],
                                    [ -1,  8, -1, 12,
                                      -1, 24, -1, 28 ],
                                    [ -1, 16, -1, 20,
                                      -1, 32, -1, 36 ],
                                    [ -1, 20, -1, 24,
                                      -1, 36, -1, 40 ],
                                    [ -1, 24, -1, 28,
                                      -1, 40, -1, 44 ],
                                    [ -1, 32, -1, 36,
                                      -1, 48, -1, 52 ],
                                    [ -1, 36, -1, 40,
                                      -1, 52, -1, 56 ],
                                    [ -1, 40, -1, 44,
                                      -1, 56, -1, 60 ],
                                    # top
                                    [ 48, 49, 52, 53,
                                      -1, -1, -1, -1 ],
                                    [ 49, 50, 53, 54,
                                      -1, -1, -1, -1 ],
                                    [ 50, 51, 54, 55,
                                      -1, -1, -1, -1 ],
                                    [ 52, 53, 56, 57,
                                      -1, -1, -1, -1 ],
                                    [ 53, 54, 57, 58,
                                      -1, -1, -1, -1 ],
                                    [ 54, 55, 58, 59,
                                      -1, -1, -1, -1 ],
                                    [ 56, 57, 60, 61,
                                      -1, -1, -1, -1 ],
                                    [ 57, 58, 61, 62,
                                      -1, -1, -1, -1 ],
                                    [ 58, 59, 62, 63,
                                      -1, -1, -1, -1 ]
                                  ] )

  ##
  # Tests scSfSv adjacency.
  ##
  def test_scSfSv( self ):
    l_scSfSvIn, l_scSfSvSend, l_scSfSvRecv = Hex.scSfSv( 0 )
    self.assertEqual( l_scSfSvIn, [] )

    self.assertEqual( l_scSfSvSend, [ [ [ 0, 2, 3, 1],
                                        [ 0, 1, 5, 4],
                                        [ 1, 3, 7, 5],
                                        [ 2, 6, 7, 3],
                                        [ 0, 4, 6, 2],
                                        [ 4, 5, 7, 6] ] ] )

    self.assertEqual( l_scSfSvRecv, [ [ [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [ 0,  1,  3,  2] ],

                                      [ [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [ 0,  4,  5,  1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1] ],

                                      [ [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [ 1,  5,  7,  3],
                                        [-1, -1, -1, -1] ],

                                      [ [-1, -1, -1, -1],
                                        [ 2,  3,  7,  6],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1] ],

                                      [ [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [ 0,  2,  6,  4],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1] ],

                                      [ [ 4,  6,  7,  5],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1],
                                        [-1, -1, -1, -1] ],
                                         ] )

  ##
  # Tests scSfSc adjacency.
  ##
  def test_scSfSc( self ):
    # first order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Hex.scSfSc( 0 )

    self.assertEqual( l_scSfScIn, [] )

    self.assertEqual( l_scSfScSend, [[1, 2, 3, 4, 5, 6]] )

    self.assertEqual( l_scSfScRecv, [ [-1, -1, -1, -1, -1,  0],
                                      [-1, -1, -1,  0, -1, -1],
                                      [-1, -1, -1, -1,  0, -1],
                                      [-1,  0, -1, -1, -1, -1],
                                      [-1, -1,  0, -1, -1, -1],
                                      [ 0, -1, -1, -1, -1, -1] ] )

    # second order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Hex.scSfSc( 1 )

    self.assertEqual( l_scSfScIn, [ [ 5, 11, 16, 21, 24, 26] ] )

    self.assertEqual( l_scSfScSend[0], [27, 36,  2,  4, 63, 10] )
    self.assertEqual( l_scSfScSend[1], [28, 37,  3,  5,  1, 11] )
    self.assertEqual( l_scSfScSend[6], [33,  4,  8, 54, 65, 20] )

    self.assertEqual( l_scSfScRecv[33], [-1, 22, -1, -1, -1, -1] )

  ##
  # Tests integration intervals for sub-cells.
  ##
  def test_intSc( self ):
    # first order
    l_intIn, l_intSend, l_intSurf = Hex.intSc( 0, ['xi0', 'xi1', 'xi2'] )

    self.assertEqual( l_intIn, [] )

    self.assertEqual( l_intSend, [ [ ( 'xi0', 0, 1 ),
                                     ( 'xi1', 0, 1 ),
                                     ( 'xi2', 0, 1 ) ] ] )

    self.assertEqual( l_intSurf, [ [ [ ( 'xi0', 0, 1 ),
                                       ( 'xi1', 0, 1 ),
                                       ( 'xi2', 0, 1 )  ] ] for l_sf in range(6) ] )

    # second order
    l_intIn, l_intSend, l_intSurf = Hex.intSc( 1, ['xi0', 'xi1', 'xi2'] )

    self.assertEqual( l_intIn, [ [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                   ( 'xi1', Fra(1,3), Fra(2,3) ),
                                   ( 'xi2', Fra(1,3), Fra(2,3) ) ] ] )

    self.assertEqual( l_intSend, [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],

                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                     
                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(0,3), Fra(1,3) ) ],


                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],

                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(0,3), Fra(1,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],


                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],

                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                   [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],


                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],

                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(2,3), Fra(3,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],


                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                   [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ],


                                   [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                     ( 'xi1', Fra(1,3), Fra(2,3) ),
                                     ( 'xi2', Fra(2,3), Fra(3,3) ) ] ] )

    self.assertEqual( l_intSurf[0:2], [ [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                            ( 'xi1', Fra(0,3), Fra(1,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                            ( 'xi1', Fra(0,3), Fra(1,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                            ( 'xi1', Fra(0,3), Fra(1,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],

                                          [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                            ( 'xi1', Fra(1,3), Fra(2,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                            ( 'xi1', Fra(1,3), Fra(2,3) ),
                                           ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                            ( 'xi1', Fra(1,3), Fra(2,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                            
                                          [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                            ( 'xi1', Fra(2,3), Fra(3,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                            ( 'xi1', Fra(2,3), Fra(3,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                          [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                            ( 'xi1', Fra(2,3), Fra(3,3) ),
                                            ( 'xi2', Fra(0,3), Fra(1,3) ) ] ],


                                         [ [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                           [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(0,3), Fra(1,3) ) ],
                                           [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(0,3), Fra(1,3) ) ],

                                           [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                           [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(1,3), Fra(2,3) ) ],
                                           [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(1,3), Fra(2,3) ) ],

                                           [ ( 'xi0', Fra(0,3), Fra(1,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                           [ ( 'xi0', Fra(1,3), Fra(2,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(2,3), Fra(3,3) ) ],
                                           [ ( 'xi0', Fra(2,3), Fra(3,3) ),
                                             ( 'xi1', Fra(0,3), Fra(1,3) ),
                                             ( 'xi2', Fra(2,3), Fra(3,3) ) ] ] ] )

  ##
  # Tests integration intervals for DG sub-faces.
  ##
  def test_intSfDg(self):
    # 1st order
    l_subs, l_intSfDg = Hex.intSfDg( 0, ['chi0', 'chi1'], ['xi0', 'xi1', 'xi2'] )

    self.assertEqual( l_subs,  ( ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 'chi1' ),
                                   ( 'xi2', 0      ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 0      ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 1      ),
                                   ( 'xi1', 'chi0' ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 1 ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 0      ),
                                   ( 'xi1', 'chi0' ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 'chi1' ),
                                   ( 'xi2', 1      ) ) ) )

    self.assertEqual( l_intSfDg, [ [ # DG face 0
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ],
                                   [ # DG face 1
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ],
                                   [ # DG face 2
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ],
                                   [ # DG face 3
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ],
                                   [ # DG face 4
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ],
                                   [ # DG face 5
                                     [ ( 'chi0', 0, 1 ),
                                       ( 'chi1', 0, 1 ) ]
                                   ] ] )

    # 2nd order
    l_subs, l_intSfDg = Hex.intSfDg( 1, ['chi0', 'chi1'], ['xi0', 'xi1', 'xi2'] )

    self.assertEqual( l_subs,  ( ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 'chi1' ),
                                   ( 'xi2', 0      ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 0      ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 1      ),
                                   ( 'xi1', 'chi0' ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 1 ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 0      ),
                                   ( 'xi1', 'chi0' ),
                                   ( 'xi2', 'chi1' ) ),

                                 ( ( 'xi0', 'chi0' ),
                                   ( 'xi1', 'chi1' ),
                                   ( 'xi2', 1      ) ) ) )

    self.assertEqual( l_intSfDg, [ [ # DG face 0
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ],


                                   [ # DG face 1
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ],


                                   [ # DG face 2
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ],


                                   [ # DG face 3
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ],


                                   [ # DG face 4
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ],


                                   [ # DG face 5
                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(0,3), Fra(1,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(1,3), Fra(2,3) ) ],


                                     [ ( 'chi0', Fra(0,3), Fra(1,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(1,3), Fra(2,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ],

                                     [ ( 'chi0', Fra(2,3), Fra(3,3) ),
                                       ( 'chi1', Fra(2,3), Fra(3,3) ) ]
                                   ] ] )

  ##
  # Tests the derivation of sub-face types.
  #
  # Sub-faces at the DG-surface have types 0-5.
  # Sub-faces not at the DG-surface have types 12-7.
  ##
  def test_scTySf(self):
    l_scTySfIn, l_scTySfSend = Hex.scTySf( 0 )

    self.assertEqual( l_scTySfIn, [] )
    self.assertEqual( l_scTySfSend, [ [0, 1, 2, 3, 4, 5] ] )

    l_scTySfIn, l_scTySfSend = Hex.scTySf( 1 )

    self.assertEqual( l_scTySfIn,   [ [12, 13, 14, 15, 16, 17] ] )
    self.assertEqual( l_scTySfSend, [ [   0,    1, 12+2, 12+3,    4, 12+5],
                                      [   0,    1, 12+2, 12+3, 12+4, 12+5],
                                      [   0,    1,    2, 12+3, 12+4, 12+5],
                                      
                                      [   0, 1+12, 12+2, 12+3,    4, 12+5],
                                      [   0, 1+12, 12+2, 12+3, 12+4, 12+5],
                                      [   0, 1+12,    2, 12+3, 12+4, 12+5],
                                      
                                      [   0, 1+12, 12+2,    3,    4, 12+5],
                                      [   0, 1+12, 12+2,    3, 12+4, 12+5],
                                      [   0, 1+12,    2,    3, 12+4, 12+5],


                                      [0+12,    1, 12+2, 12+3,    4, 12+5],
                                      [0+12,    1, 12+2, 12+3, 12+4, 12+5],
                                      [0+12,    1,    2, 12+3, 12+4, 12+5],

                                      [0+12,    1, 12+2, 12+3,    4,    5],
                                      [0+12,    1, 12+2, 12+3, 12+4,    5],
                                      [0+12,    1,    2, 12+3, 12+4,    5],


                                      [0+12, 1+12,    2, 12+3, 12+4, 12+5],
                                      [0+12, 1+12,    2,    3, 12+4, 12+5],
                                      [0+12, 1+12,    2, 12+3, 12+4,    5],
                                      [0+12, 1+12,    2,    3, 12+4,    5],

                                      [0+12, 1+12, 12+2,    3,    4, 12+5],
                                      [0+12, 1+12, 12+2,    3, 12+4, 12+5],
                                      [0+12, 1+12, 12+2,    3,    4,    5],
                                      [0+12, 1+12, 12+2,    3, 12+4,    5],


                                      [0+12, 1+12, 12+2, 12+3,    4, 12+5],
                                      [0+12, 1+12, 12+2, 12+3,    4,    5],


                                      [0+12, 1+12, 12+2, 12+3, 12+4,    5] ] )
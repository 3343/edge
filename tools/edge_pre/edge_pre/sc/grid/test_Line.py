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
# Unit tests for the line sub-grid.
##
import unittest
import fractions
from . import Line

class TestGridLine( unittest.TestCase ):
  ##
  # Tests generation of vertices.
  ##
  def test_svs(self):
    # FV
    l_svs = Line.svs( 0 )
    self.assertEqual( len(l_svs), 2 )
    for l_ve in l_svs:
      self.assertEqual( len(l_ve), 1 )

    self.assertEqual( l_svs[0][0], 0 )
    self.assertEqual( l_svs[1][0], 1 )

    # 2nd order
    l_svs = Line.svs( 1 )
    self.assertEqual( len(l_svs), 4 )
    for l_ve in l_svs:
      self.assertEqual( len(l_ve), 1 )

    self.assertEqual( l_svs[0][0], 0 )
    self.assertEqual( l_svs[1][0], fractions.Fraction(1, 3) )
    self.assertEqual( l_svs[2][0], fractions.Fraction(2, 3) )
    self.assertEqual( l_svs[3][0], 1 )

    # 3rd order
    l_svs = Line.svs( 2 )
    self.assertEqual( len(l_svs), 6 )
    for l_ve in l_svs:
      self.assertEqual( len(l_ve), 1 )

    self.assertEqual( l_svs[0][0], 0 )
    self.assertEqual( l_svs[1][0], fractions.Fraction(1, 5) )
    self.assertEqual( l_svs[2][0], fractions.Fraction(2, 5) )
    self.assertEqual( l_svs[3][0], fractions.Fraction(3, 5) )
    self.assertEqual( l_svs[4][0], fractions.Fraction(4, 5) )
    self.assertEqual( l_svs[5][0], 1 )

  ##
  # Tests scSv adjacency.
  ##
  def test_scSv(self):
    # 2nd order
    l_scSvIn, l_scSvSend, l_scSvRecv = Line.scSv( 1 )
    self.assertEqual( len(l_scSvIn), 1 )
    self.assertEqual( len(l_scSvIn[0]), 2 )

    self.assertEqual( len(l_scSvSend), 2 )
    self.assertEqual( l_scSvSend[0], [0,1] )
    self.assertEqual( l_scSvSend[1], [2,3] )

    self.assertEqual( len(l_scSvRecv), 2 )
    self.assertEqual( l_scSvRecv[0], [0,-1] )
    self.assertEqual( l_scSvRecv[1], [3,-1] )

    # 3rd order
    l_scSvIn, l_scSvSend, l_scSvRecv = Line.scSv( 2 )
    self.assertEqual( len(l_scSvIn), 3 )
    for l_sc in l_scSvIn:
      self.assertEqual( len(l_sc), 2 )

    self.assertEqual( l_scSvIn[0], [1,2] )
    self.assertEqual( l_scSvIn[1], [2,3] )
    self.assertEqual( l_scSvIn[2], [3,4] )

    self.assertEqual( len(l_scSvSend), 2 )
    self.assertEqual( l_scSvSend[0], [0,1] )
    self.assertEqual( l_scSvSend[1], [4,5] )

    self.assertEqual( len(l_scSvRecv), 2 )
    self.assertEqual( l_scSvRecv[0], [0,-1] )
    self.assertEqual( l_scSvRecv[1], [5,-1] )

  ##
  # Tests scSfSc adjacency.
  ##
  def test_scSfSc(self):
    # second order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Line.scSfSc( 1 )

    self.assertEqual( len(l_scSfScIn), 1 )
    self.assertEqual( l_scSfScIn[0], [1,2] )

    self.assertEqual( len(l_scSfScSend), 2 )
    self.assertEqual( l_scSfScSend[0], [3,0] )
    self.assertEqual( l_scSfScSend[1], [0,4] )

    self.assertEqual( len(l_scSfScRecv), 2 )
    self.assertEqual( l_scSfScRecv[0], [1,-1] )
    self.assertEqual( l_scSfScRecv[1], [2,-1] )

    # third order
    l_scSfScIn, l_scSfScSend, l_scSfScRecv = Line.scSfSc( 2 )

    self.assertEqual( len(l_scSfScIn), 3 )
    self.assertEqual( l_scSfScIn[0], [3,1] )
    self.assertEqual( l_scSfScIn[1], [0,2] )
    self.assertEqual( l_scSfScIn[2], [1,4] )

    self.assertEqual( len(l_scSfScSend), 2 )
    self.assertEqual( l_scSfScSend[0], [5,0] )
    self.assertEqual( l_scSfScSend[1], [2,6] )

    self.assertEqual( len(l_scSfScRecv), 2 )
    self.assertEqual( l_scSfScRecv[0], [3,-1] )
    self.assertEqual( l_scSfScRecv[1], [4,-1] )

  ##
  # Tests the derivation of sub-face types.
  #
  # Sub-faces at the DG-surface have types 0-1.
  # Sub-faces not at the DG-surface have type 2.
  ##
  def test_scTySf(self):
    # get sub-face types for order 2
    l_scTySfIn, l_scTySfSend = Line.scTySf( 1 )

    self.assertEqual( l_scTySfIn, [ [4,5] ] )
    self.assertEqual( l_scTySfSend, [ [0,5], [4,1] ] )

    # get sub-face types for order 3
    l_scTySfIn, l_scTySfSend = Line.scTySf( 2 )

    self.assertEqual( l_scTySfIn, [ [4,5], [4,5], [4,5] ] )
    self.assertEqual( l_scTySfSend, [ [0,5], [4,1] ] )

    # get sub-face types for order 4
    l_scTySfIn, l_scTySfSend = Line.scTySf( 3 )

    self.assertEqual( l_scTySfIn, [ [4,5], [4,5], [4,5], [4,5], [4,5] ] )
    self.assertEqual( l_scTySfSend, [ [0,5], [4,1] ] )

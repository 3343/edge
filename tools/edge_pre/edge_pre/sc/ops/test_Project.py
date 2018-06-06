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
# Unit tests for the projection operators.
##
import unittest
from . import Project
import edge_pre.sc.grid.Line
import edge_pre.dg.basis.Line
import edge_pre.sc.grid.Tria
import edge_pre.dg.basis.Tria
import edge_pre.sc.grid.Tet
import edge_pre.dg.basis.Tet
from fractions import Fraction as Fra
import sympy

class TestProject( unittest.TestCase ):
  ##
  # Tests derivation of the scatter operator for line elements.
  ##
  def test_scatterLine(self):
    #
    # 2nd order
    #
    l_syms, l_basis = edge_pre.dg.basis.Line.gen( 1 )
    l_int = [ (l_syms[0], 0, 1) ]
    l_mapsSc, l_detsSc = edge_pre.sc.grid.Line.intSc( 1, l_syms )

    # inner sub-cells
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[0],
                                 l_detsSc[0] )
    l_scatterUt = sympy.Matrix([[1], [0]])
    self.assertEqual( l_scatter, l_scatterUt )

    # send sub-cells
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[1],
                                 l_detsSc[1] )
    l_scatterUt = sympy.Matrix( [ [1,          1        ],
                                  [-Fra(2, 3), Fra(2, 3)]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    # surf sub-cells, face 1
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[2][0],
                                 l_detsSc[2][0] )
    l_scatterUt = sympy.Matrix( [ [1,         ],
                                  [-Fra(2, 3) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    # surf sub-cells, face2
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[2][1],
                                 l_detsSc[2][1] )
    l_scatterUt = sympy.Matrix( [ [1,                        ],
                                  [Fra(2, 3) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    #
    # 3rd order
    #
    l_syms, l_basis = edge_pre.dg.basis.Line.gen( 2 )
    l_int = [ (l_syms[0], 0, 1) ]
    l_mapsSc, l_detsSc = edge_pre.sc.grid.Line.intSc( 2, l_syms )

    # inner sub-cells
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[0],
                                 l_detsSc[0] )
    l_scatterUt = sympy.Matrix( [ [  1,           1,           1         ],
                                  [ -Fra(2, 5),   0,           Fra(2, 5) ],
                                  [ -Fra(6, 25), -Fra(12,25), -Fra(6,25) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    # send sub-cells
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[1],
                                 l_detsSc[1] )
    l_scatterUt = sympy.Matrix( [ [  1,           1          ],
                                  [ -Fra(4, 5),   Fra(4, 5)  ],
                                  [  Fra(12, 25), Fra(12,25) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    # surf sub-cells, face 1
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[2][0],
                                 l_detsSc[2][0] )
    l_scatterUt = sympy.Matrix( [ [  1           ],
                                  [ -Fra(4, 5)   ],
                                  [  Fra(12, 25) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

    # surf sub-cells, face 2
    l_scatter = Project.scatter( l_syms,
                                 l_basis,
                                 l_int,
                                 l_mapsSc[2][1],
                                 l_detsSc[2][1] )
    l_scatterUt = sympy.Matrix( [ [  1           ],
                                  [  Fra(4, 5)   ],
                                  [  Fra(12, 25) ]
                                ])
    self.assertEqual( l_scatter, l_scatterUt )

  ##
  # Tests derivation of the gather operator for line elements.
  ##
  def test_gatherLine(self):
    #
    # 2nd order
    #
    l_syms, l_basis = edge_pre.dg.basis.Line.gen( 1 )
    l_int = [ (l_syms[0], 0, 1) ]
    l_mapsSc, l_detsSc = edge_pre.sc.grid.Line.intSc( 1, l_syms )

    l_gather = Project.gather( l_syms,
                               l_basis,
                               l_int,
                               l_mapsSc[0]+l_mapsSc[1],
                               l_detsSc[0]+l_detsSc[1] )

    l_gatherUt = sympy.Matrix( [ [  Fra(1, 3),          0 ],
                                 [  Fra(1, 3), -Fra(3, 4) ],
                                 [  Fra(1, 3),  Fra(3, 4) ]
                               ])
    self.assertEqual( l_gather, l_gatherUt )

    #
    # 3rd order
    #
    l_syms, l_basis = edge_pre.dg.basis.Line.gen( 2 )
    l_int = [ (l_syms[0], 0, 1) ]
    l_mapsSc, l_detsSc = edge_pre.sc.grid.Line.intSc( 2, l_syms )

    l_gather = Project.gather( l_syms,
                               l_basis,
                               l_int,
                               l_mapsSc[0]+l_mapsSc[1],
                               l_detsSc[0]+l_detsSc[1] )

    l_gatherUt = sympy.Matrix( [ [  Fra(1, 5), -Fra(1, 4), -Fra(25, 84) ],
                                 [  Fra(1, 5),  0,         -Fra(25, 42) ],
                                 [  Fra(1, 5),  Fra(1, 4), -Fra(25, 84) ],
                                 [  Fra(1, 5), -Fra(1, 2),  Fra(25, 42) ],
                                 [  Fra(1, 5),  Fra(1, 2),  Fra(25, 42) ]
                               ])
    self.assertEqual( l_gather, l_gatherUt )

  ##
  # Tests that the gather operator inverts the scatter operator for line elements.
  ##
  def test_gatherScatterInvLine(self):
    for l_de in range(1,5):
      l_syms, l_basis = edge_pre.dg.basis.Line.gen( l_de )
      l_int = [ (l_syms[0], 0, 1) ]
      l_mapsSc, l_detsSc = edge_pre.sc.grid.Line.intSc( l_de, l_syms )

      l_scatter = Project.scatter( l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )

      l_gather  = Project.gather(  l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )

      self.assertEqual( l_scatter*l_gather, sympy.eye(len(l_basis)) )

  ##
  # Tests that the gather operator inverts the scatter operator for triangles.
  ##
  def test_gatherScatterInvTria(self):
    for l_de in range(1,3):
      l_syms, l_basis = edge_pre.dg.basis.Tria.gen( l_de )
      l_int = [ (l_syms[0], 0, 1-l_syms[1]), (l_syms[1], 0, 1) ]
      l_mapsSc, l_detsSc = edge_pre.sc.grid.Tria.intSc( l_de, l_syms )

      l_scatter = Project.scatter( l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )
      l_gather  = Project.gather(  l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )

      self.assertEqual( l_scatter*l_gather, sympy.eye(len(l_basis)) )


  ##
  # Tests that the gather operator inverts the scatter operator for tets.
  ##
  def test_gatherScatterInvTet(self):
    for l_de in range(0,1):
      l_syms, l_basis = edge_pre.dg.basis.Tet.gen( l_de )
      l_int = [ (l_syms[0], 0, 1-l_syms[1]-l_syms[2]), (l_syms[1], 0, 1-l_syms[2]), (l_syms[2], 0, 1) ]
      l_mapsSc, l_detsSc = edge_pre.sc.grid.Tet.intSc( l_de, l_syms )

      l_scatter = Project.scatter( l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )

      l_gather  = Project.gather(  l_syms,
                                   l_basis,
                                   l_int,
                                   l_mapsSc[0]+l_mapsSc[1],
                                   l_detsSc[0]+l_detsSc[1] )

      self.assertEqual( l_scatter*l_gather, sympy.eye(len(l_basis)) )
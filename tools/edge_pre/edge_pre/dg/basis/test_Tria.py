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
# Tests derivation of basis functions for triangles.
##
import unittest
from . import Tria
import sympy
from fractions import Fraction as Fra

class TestTria( unittest.TestCase ):
  ##
  # Tests generation of basis functions.
  ##
  def test_gen(self):
    l_xi1 = sympy.symbols('xi_1')
    l_xi2 = sympy.symbols('xi_2')


    # P4, hierarchical storage
    l_basisUt = [ sympy.sympify(1),
                  3*l_xi2 - 1,
                  2*l_xi1 + l_xi2 - 1,
                  10*l_xi2**2 - 8*l_xi2 + 1,
                  (5*l_xi2 - 1)*(2*l_xi1 + l_xi2 - 1),
                  -(l_xi2 - 1)**2/2 + 3*(2*l_xi1 + l_xi2 - 1)**2/2,
                  35*l_xi2**3 - 45*l_xi2**2 + 15*l_xi2 - 1,
                  (2*l_xi1 + l_xi2 - 1)*(36*l_xi2 + 21*(2*l_xi2 - 1)**2 - 17)/4,
                  (7*l_xi2 - 1)*(-(l_xi2 - 1)**2 + 3*(2*l_xi1 + l_xi2 - 1)**2)/2,
                  -(6*l_xi1*(l_xi2 - 1)**3 + 3*(l_xi2 - 1)**4 - 5*(l_xi2 - 1)*(2*l_xi1 + l_xi2 - 1)**3)/(2*l_xi2 - 2),
                  126*l_xi2**4 - 224*l_xi2**3 + 126*l_xi2**2 - 24*l_xi2 + 1,
                  (2*l_xi1 + l_xi2 - 1)*(21*(2*l_xi2 - 1)**3 + 21*(2*l_xi2 - 1)**2 - 2)/2,
                  (-(l_xi2 - 1)**2 + 3*(2*l_xi1 + l_xi2 - 1)**2)*(20*l_xi2 + 9*(2*l_xi2 - 1)**2 - 8)/2,
                  -(9*l_xi2 - 1)*(6*l_xi1*(l_xi2 - 1)**3 + 3*(l_xi2 - 1)**4 - 5*(l_xi2 - 1)*(2*l_xi1 + l_xi2 - 1)**3)/(2*l_xi2 - 2),
                  3*(l_xi2 - 1)**4/8 - 15*(l_xi2 - 1)**2*(2*l_xi1 + l_xi2 - 1)**2/4 + 35*(2*l_xi1 + l_xi2 - 1)**4/8 ]

    l_syms, l_basis = Tria.gen( 0 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( l_basis, [ 1 ] )

    l_syms, l_basis = Tria.gen( 1 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( len(l_basis), 3 )
    for l_ba in range(3):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

    l_syms, l_basis = Tria.gen( 2 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( len(l_basis), 6 )
    for l_ba in range(6):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

    l_syms, l_basis = Tria.gen( 3 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( len(l_basis), 10 )
    for l_ba in range(10):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

    l_syms, l_basis = Tria.gen( 4 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( len(l_basis), 15 )
    for l_ba in range(15):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

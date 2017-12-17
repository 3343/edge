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
# Tests derivation of basis functions for tetrahedrons.
##
import unittest
from . import Tet
import sympy
from fractions import Fraction as Fra

class TestTria( unittest.TestCase ):
  ##
  # Tests generation of basis functions.
  ##
  def test_gen(self):
    l_xi1 = sympy.symbols('xi_1')
    l_xi2 = sympy.symbols('xi_2')
    l_xi3 = sympy.symbols('xi_3')

    # P3, hierarchical storage
    l_basisUt = [ sympy.sympify(1),
                  4*l_xi3 - 1,
                  3*l_xi2 + l_xi3 - 1,
                  2*l_xi1 + l_xi2 + l_xi3 - 1,
                  15*l_xi3**2 - 10*l_xi3 + 1,
                  (6*l_xi3 - 1)*(3*l_xi2 + l_xi3 - 1),
                  10*l_xi2**2 + 8*l_xi2*l_xi3 - 8*l_xi2 + l_xi3**2 - 2*l_xi3 + 1,
                  (6*l_xi3 - 1)*(2*l_xi1 + l_xi2 + l_xi3 - 1),
                  (5*l_xi2 + l_xi3 - 1)*(2*l_xi1 + l_xi2 + l_xi3 - 1),
                  -Fra(1,2)*(l_xi2 + l_xi3 - 1)**2 + Fra(3,2)*(2*l_xi1 + l_xi2 + l_xi3 - 1)**2,
                  56*l_xi3**3 - 63*l_xi3**2 + 18*l_xi3 - 1,
                  (3*l_xi2 + l_xi3 - 1)*(14*l_xi3 + 7*(2*l_xi3 - 1)**2 - 6),
                  -(8*l_xi3 - 1)*(4*l_xi2*(l_xi3 - 1)**2 + 3*(l_xi3 - 1)**3 - 5*(l_xi3 - 1)*(2*l_xi2 + l_xi3 - 1)**2)/(2*l_xi3 - 2),
                  35*l_xi2**3 + 45*l_xi2**2*l_xi3 - 45*l_xi2**2 + 15*l_xi2*l_xi3**2 - 30*l_xi2*l_xi3 + 15*l_xi2 + l_xi3**3 - 3*l_xi3**2 + 3*l_xi3 - 1,
                  (14*l_xi3 + 7*(2*l_xi3 - 1)**2 - 6)*(2*l_xi1 + l_xi2 + l_xi3 - 1),
                  (8*l_xi3 - 1)*(5*l_xi2 + l_xi3 - 1)*(2*l_xi1 + l_xi2 + l_xi3 - 1),
                  -(36*l_xi2*(l_xi3 - 1)**2 + 17*(l_xi3 - 1)**3 - 21*(l_xi3 - 1)*(2*l_xi2 + l_xi3 - 1)**2)*(2*l_xi1 + l_xi2 + l_xi3 - 1)/(4*l_xi3 - 4),
                  Fra(1,2)*(8*l_xi3 - 1)*(-(l_xi2 + l_xi3 - 1)**2 + 3*(2*l_xi1 + l_xi2 + l_xi3 - 1)**2),
                  Fra(1,2)*(-(l_xi2 + l_xi3 - 1)**2 + 3*(2*l_xi1 + l_xi2 + l_xi3 - 1)**2)*(7*l_xi2 + l_xi3 - 1),
                  -(6*l_xi1*(l_xi2 + l_xi3 - 1)**3 + 3*(l_xi2 + l_xi3 - 1)**4 - 5*(l_xi2 + l_xi3 - 1)*(2*l_xi1 + l_xi2 + l_xi3 - 1)**3)/(2*l_xi2 + 2*l_xi3 - 2) ]

    # check FV basis
    l_syms, l_basis = Tet.gen( 0 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2, l_xi3 ] )
    self.assertEqual( l_basis, [ sympy.sympify(1) ] )

    # check degree 1
    l_syms, l_basis = Tet.gen( 1 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2, l_xi3 ] )
    self.assertEqual( len(l_basis), 4 )
    for l_ba in range(4):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

    # check degree 2
    l_syms, l_basis = Tet.gen( 2 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2, l_xi3 ] )
    self.assertEqual( len(l_basis), 10 )
    for l_ba in range(10):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

    # check degree 3
    l_syms, l_basis = Tet.gen( 3 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2, l_xi3 ] )
    self.assertEqual( len(l_basis), 20 )
    for l_ba in range(20):
      self.assertEqual( sympy.simplify( l_basis[l_ba]-l_basisUt[l_ba] ), 0 )

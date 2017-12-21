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
# Tests derivation of basis functions for quadrilaterals
##
import unittest
from . import Quad
import sympy
from fractions import Fraction as Fra

class TestQuad( unittest.TestCase ):
  ##
  # Tests generation of basis functions.
  ##
  def test_gen(self):
    l_xi1 = sympy.symbols('xi_1')
    l_xi2 = sympy.symbols('xi_2')

    l_syms, l_basis = Quad.gen( 0 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( l_basis, [ 1 ] )

    l_syms, l_basis = Quad.gen( 1 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( l_basis, [ 1, 2*l_xi1-1, 2*l_xi2-1, (2*l_xi1-1)*(2*l_xi2-1) ] )

    l_syms, l_basis = Quad.gen( 2 )
    self.assertEqual( l_syms,  [ l_xi1, l_xi2 ] )
    self.assertEqual( l_basis, [ sympy.sympify(1),
                                 2*l_xi1-1,
                                 Fra(3,2)*(2*l_xi1-1)**2 - Fra(1,2),
                                 2*l_xi2-1,
                                 (2*l_xi1-1)*(2*l_xi2-1),
                                 (2*l_xi2-1)*( Fra(3,2) * (2*l_xi1-1)**2 - Fra(1,2) ),
                                 Fra(3,2)*(2*l_xi2-1)**2 - Fra(1,2),
                                 (2*l_xi1-1)*( Fra(3,2) * (2*l_xi2-1)**2 - Fra(1,2) ),
                                   ( ( Fra(3,2) * (2*l_xi2-1)**2 - Fra(1,2) ) )
                                 * ( ( Fra(3,2) * (2*l_xi1-1)**2 - Fra(1,2) ) ) ] )

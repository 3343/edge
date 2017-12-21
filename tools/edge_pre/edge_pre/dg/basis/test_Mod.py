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
# Tests derivation of basis functions for line elements.
##
import unittest
from . import Mod
import sympy

class TestMod( unittest.TestCase ):
  ##
  # Tests generation of basis functions.
  ##
  def test_unify(self):
    l_xi1 = sympy.symbols('xi_1')
    l_xi2 = sympy.symbols('xi_2')
    l_chi1 = sympy.symbols('xi_1')
    l_chi2 = sympy.symbols('chi_2')
    l_chi3 = sympy.symbols('xi_2')

    l_basis1 = [ sympy.sympify(1),
                 2*l_xi1-1,
                 (2*l_xi2 -1)**2,
                 (2*l_xi1 -1)**3  ]

    l_basis2 = [ sympy.sympify(1),
                 l_chi1,
                 l_chi1 + 1,
                 (2*l_chi3 -1),
                 l_chi2  ]

    # unify
    l_su, l_bu = Mod.unify( [l_xi1, l_xi2], l_basis1,
                            [l_chi1, l_chi2, l_chi3],       l_basis2 )
    
    # check the results
    self.assertEqual( l_su[0], sympy.symbols('xi_1_unique') )
    self.assertEqual( l_su[1], sympy.symbols('chi_2') )
    self.assertEqual( l_bu, [ sympy.sympify(1),
                              sympy.symbols('xi_1_unique'),
                              sympy.symbols('xi_1_unique') + 1,
                              2* sympy.symbols('xi_2_unique') - 1,
                              sympy.symbols('chi_2') ] )
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018, Regents of the University of California
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
# Tests the scalar integration.
##
import unittest
from . import Scalar
import sympy

class TestScalar( unittest.TestCase ):
  ##
  # Tests the inexact integration.
  ##
  def test_quad(self):
    l_x = sympy.symbols('x')
    l_y = sympy.symbols('y')
    l_z = sympy.symbols('z')

    l_fun = 1-l_x*l_x-l_y*l_y
    l_int = [ (l_x, 0, sympy.sqrt(1-l_y**2)), (l_y, 0, 1) ]

    self.assertAlmostEqual( (sympy.pi/8).evalf(), Scalar.intQuad( l_fun, l_int )[0] )

    # Wolfram Alpha:
    # integrate 1-x2+y4+x-y3 dx dy for x from 0 to sin(1-y2) for y from 1 to 3
    l_fun = 1+l_x-l_x*l_x-l_y*l_y*l_y+l_y*l_y*l_y*l_y
    l_int = [ (l_x, 0, sympy.sin(1-l_y**2)), (l_y, 1, 3) ]
  
    self.assertAlmostEqual( Scalar.intQuad( l_fun, l_int )[0], -2.55741, 5 )


    # Wolfram Alpha:
    # integrate 1-x2+y4+x-y3+z2-z dx dy dz for x from 0 to sin(1-y2+z) for y from 1 to 3-z for z from 0 to 1
    l_fun = 1+l_x-l_x*l_x-l_y*l_y*l_y+l_y*l_y*l_y*l_y+l_z*l_z-l_z
    l_int = [ (l_x, 0, sympy.sin(1-l_y**2+l_z)), (l_y, 1, 3-l_z), (l_z, 0, 1) ]   

    self.assertAlmostEqual( Scalar.intQuad( l_fun, l_int )[0], 1.26587, 5 )

  ##
  # Tests the mixed integration (using character cut-off to switch between exact and inexact).
  ##
  def test_mix(self):
    l_x = sympy.symbols('x')
    l_y = sympy.symbols('y')

    l_fun = 1-l_x*l_x-l_y*l_y
    l_int = [ (l_x, 0, sympy.sqrt(1-l_y**2)), (l_y, 0, 1) ]

    l_res = Scalar.int( l_fun, l_int, 20 )
    self.assertEqual( l_res[0], sympy.pi/8 )
    self.assertEqual( l_res[1], True )

    l_res = Scalar.int( l_fun, l_int, 10 )
    self.assertAlmostEqual( l_res[0], (sympy.pi/8).evalf() )
    self.assertEqual( l_res[1], False )
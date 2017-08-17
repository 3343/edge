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
# Unit tests for the conversion of arrays to string.
##
import unittest
from . import ArrStr

class TestArrStr( unittest.TestCase ):
  ##
  # Tests conversion of 2d float.
  ##
  def test_float2d(self):
    l_dat = [ [0.0, 1.0], [2.0, 3.0], [0.0, 1.0], [4.0, 5.0] ]
    l_str = ArrStr.float2d( l_dat )
    l_strUt =\
"""0.0,1.0,
2.0,3.0,
0.0,1.0,
4.0,5.0,
"""

    l_dat = [ [0.0, 1.0, 5.0], [2.0, 3.0, 7.0], [0.0, 1.0, 8.0], [4.0, 5.0, 9.0] ]
    l_str = ArrStr.float2d( l_dat )
    l_strUt =\
"""0.0,1.0,5.0,
2.0,3.0,7.0,
0.0,1.0,8.0,
4.0,5.0,9.0,
"""
    self.assertEqual( l_str, l_strUt )

  ##
  # Tests conversion of 2d int.
  ##
  def test_int2d(self):
    l_dat = [ [1,2,3,4], [4,3,2,1], [1,3,4,2], [4,5,1,5], [5,9,2,4] ]
    l_str = ArrStr.int2d( l_dat )
    l_strUt =\
"""1,2,3,4,
4,3,2,1,
1,3,4,2,
4,5,1,5,
5,9,2,4,
"""
    self.assertEqual( l_str, l_strUt )

  ##
  # Tests conversion of 3d int.
  ##
  def test_int3d(self):
    l_dat = [
              [ [1,2,3], [3,4,5], [8,1,5], [9,2,2] ],
              [ [9,6,2], [9,5,2], [9,5,1], [5,1,9] ],
              [ [9,5,2], [8,5,2], [5,1,4], [2,4,2] ]
            ];
    l_str = ArrStr.int3d( l_dat )
    l_strUt =\
"""1,2,3,
3,4,5,
8,1,5,
9,2,2,

9,6,2,
9,5,2,
9,5,1,
5,1,9,

9,5,2,
8,5,2,
5,1,4,
2,4,2,

"""
    self.assertEqual( l_str, l_strUt )

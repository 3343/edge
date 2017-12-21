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
# Basis for hexes.
##
import sympy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages
from . import Line

##
# Generates basis functions for hexes.
#
# @param i_deg polynomial degree.
# @return two lists, 1) symbols 2) basis functions
##
def gen( i_deg ):
  # get line basis
  l_lineS, l_lineB = Line.gen( i_deg )

  l_syms = [ sympy.symbols('xi_1'), sympy.symbols('xi_2'), sympy.symbols('xi_3') ]

  # build hex basis functions as tensor product
  l_basis = []
  for l_d3 in range(i_deg+1):
    for l_d2 in range(i_deg+1):
      for l_d1 in range(i_deg+1):
        l_basis = l_basis + [ l_lineB[l_d1].subs( l_lineS[0], l_syms[0] ) *
                              l_lineB[l_d2].subs( l_lineS[0], l_syms[1] ) *
                              l_lineB[l_d3].subs( l_lineS[0], l_syms[2] ) ]

  return l_syms, l_basis
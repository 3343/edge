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
# Basis for triangles.
##
import sympy
import fractions
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages

##
# Generates basis functions for triangles.
#
# @param i_deg polynomial degree.
# @return two lists, 1) symbols 2) basis functions
##
def gen( i_deg ):
  l_xi1 = sympy.symbols('xi_1') # first triangular coords
  l_xi2 = sympy.symbols('xi_2') # second triangular coord

  l_eta1 = sympy.symbols('l_eta_1') # first collapsed coord
  l_eta2 = sympy.symbols('l_eta_2') # scond collapsed coord

  # basis
  l_basis = []

  # derive basis functions
  for l_de in range(0, i_deg+1):
    for l_p in range(0, l_de+1):
      # first principal function
      l_psiA = sympy.jacobi( l_p, 0, 0, l_eta1 )
      
      for l_q in range(0, l_de+1):
        if( l_p + l_q == l_de ): # build hierarchical basis following pascals triangle
          # second principal function
          l_psiB = ( ( (1-l_eta2)/2 )**l_p * sympy.jacobi( l_q, 2*l_p+1, 0, l_eta2 ) )

          l_basis = l_basis + [l_psiA * l_psiB]
          
          # insert triangular coords
          l_basis[-1] = l_basis[-1].subs( l_eta1, 2*(1+l_xi1)/(1-l_xi2) - 1 ).subs(l_eta2, l_xi2)
          
          # use triangle with xi1 in [0,1], xi2 in [0,1-xi1]
          l_basis[-1] = l_basis[-1].subs(l_xi1, 2*l_xi1-1).subs(l_xi2, 2*l_xi2 - 1)
          
          # simplify the basis
          l_basis[-1] = sympy.simplify( l_basis[-1] )

          # cancel nominators
          l_basis[-1] = sympy.cancel( l_basis[-1] )

  return [l_xi1, l_xi2], l_basis

##
# Plots basis functions for quadrilaterals
#
# @param i_out output path.
# @param i_syms symbols used in basis.
# @param i_funs basis functions.
##
def plot( i_out, i_syms, i_funs ):
  assert( len(i_syms) == 2 )

  # scaling for the parametric plot
  l_sca = (1-i_syms[1])

  # open pdf
  with PdfPages(i_out) as l_pdf:
    for l_fu in i_funs:
      sympy.plotting.plot3d_parametric_surface(  i_syms[0]*l_sca,
                                                 i_syms[1],
                                                 l_fu.subs( i_syms[0], i_syms[0]*l_sca ),
                                                (i_syms[0], 0, 1),
                                                (i_syms[1], 0, 1) )
      l_pdf.savefig()
      matplotlib.pyplot.close()

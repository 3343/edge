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
# Basis for tetrahedrons.
##
import sympy

##
# Generates basis functions for tetrahedrons.
#
# @param i_deg polynomial degree.
# @return two lists, 1) symbols 2) basis functions
##
def gen( i_deg ):
  l_xi1 = sympy.symbols('xi_1') # first tetrahedral coord
  l_xi2 = sympy.symbols('xi_2') # second tetrahedral coord
  l_xi3 = sympy.symbols('xi_3') # third tetrahedral coord

  l_eta1 = sympy.symbols('eta_1') # first collapsed coord
  l_eta2 = sympy.symbols('eta_2') # second collapsed coord
  l_eta3 = sympy.symbols('eta_3') # third collapsed coord

  # basis
  l_basis = []

  # derive basis functions
  for l_de in range(0, i_deg+1):
    for l_p1 in range(0, l_de+1):
      # first principal function
      l_psiA = sympy.jacobi( l_p1, 0, 0, l_eta1 )

      for l_p2 in range(0, l_de+1):
        # second principal function
        l_psiB = ( ( (1-l_eta2)/2 )**l_p1 * sympy.jacobi( l_p2, 2*l_p1+1, 0, l_eta2 ) )

        for l_p3 in range(0, l_de+1):
          if( l_p1 + l_p2 + l_p3 == l_de ): # build hierarchical basis following pascals triangle
            # third principal function
            l_psiC = ( (1-l_eta3)/2 )**(l_p1 + l_p2) * sympy.jacobi( l_p3, 2*l_p1+2*l_p2+2, 0, l_eta3 )

            l_basis = l_basis + [l_psiA * l_psiB * l_psiC]
            
            # insert tetrahedral coordinates
            l_basis[-1] = l_basis[-1].subs( l_eta1, 2*(1+l_xi1) / (-l_xi2-l_xi3) - 1 )
            l_basis[-1] = l_basis[-1].subs( l_eta2, 2*(1+l_xi2) / (1-l_xi3)      - 1 )
            l_basis[-1] = l_basis[-1].subs( l_eta3, l_xi3                            )

            # use tets with xi1, xi2, xi3 in [0,1]
            l_basis[-1] = l_basis[-1].subs(l_xi1, 2*l_xi1-1)
            l_basis[-1] = l_basis[-1].subs(l_xi2, 2*l_xi2-1)
            l_basis[-1] = l_basis[-1].subs(l_xi3, 2*l_xi3-1)
        
            # simplify basis
            l_basis[-1] = sympy.simplify( l_basis[-1] )

            # cancel denominators
            l_basis[-1] = sympy.cancel( l_basis[-1] )

  return [l_xi1, l_xi2, l_xi3], l_basis

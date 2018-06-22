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
# Scalar integration.
##
import sympy
import scipy.integrate

##
# Performs a numerical integration of the given sympy-function and intervals.
#
# @param i_fun function which is integrated.
# @param i_int multi-dimensional integration intervals. The least dependent intervals comes last.
##
def intQuad( i_fun,
             i_int ):
  # extract variables and intervals in descending order
  l_varNames = []
  l_varInts  = []

  for l_in in i_int:
    l_varNames = l_varNames + [ l_in[0]  ]
    l_varInts  = l_varInts  + [ l_in[1:] ]

  # lambdify the integrand
  l_funLam = sympy.lambdify( l_varNames, i_fun )

  # lambdify the integrals
  l_intLam = []
  for l_li in range( len(l_varNames) ):
    # assemble dependent variables
    l_varDeps = []
    for l_de in range( l_li+1, len(l_varNames) ):
      l_varDeps = l_varDeps + [ l_varNames[l_de] ]

    # set both bounds for this variable
    l_intLam = l_intLam + [ sympy.lambdify(l_varDeps, l_varInts[l_li]) ]

  # perform the quadrature
  return scipy.integrate.nquad( l_funLam,
                                l_intLam )

##
# Integrates the given function exactly or using quadrature.
#
# @param i_fun function, which is integrated.
# @param i_int integration intervals.
# @param i_quad if the number of characters in the string-representation of the function is large, inexact integration is used.
#
# @return 0: integration result, 1: true if exact, false otherwise
##
def int( i_fun,
         i_int,
         i_quad=20 ):
  l_exact = ( len(str(i_fun)) < i_quad )

  if( l_exact ):
    l_res = sympy.integrate( i_fun, *i_int )
  else:
    l_res = intQuad( i_fun, i_int )[0]

  return l_res, l_exact
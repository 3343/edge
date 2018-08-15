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
# Annotation of entities through expressions.
##
import numpy

##
# Evaluates the expression for the given points
#
# @param i_nVars number of variables.
# @param i_expr expression which gets evaluated.
# @param i_ptCrds coordinates of the points.
##
def evalPt( i_nVars,
            i_expr,
            i_ptCrds ):
  # determine result values
  l_vals = numpy.zeros( (i_nVars, len(i_ptCrds) ), dtype = 'float64')

  l_comp = compile( i_expr, '<string>', 'exec' )

  for l_pt in range( len(i_ptCrds) ):
    l_vars = { 'x': i_ptCrds[l_pt, 0],
               'y': i_ptCrds[l_pt, 1],
               'z': i_ptCrds[l_pt, 2] }

    exec( l_comp, l_vars )
    l_vals[:, l_pt] = l_vars['q']

  return l_vals

##
# Evaluates the expression for the given entities.
# The average values of the vertices will be used.
#
# @param i_nVars number of variables.
# @param i_expr expression which gets evaluated.
# @param i_enVe vertices adjacent to the entities.
# @param i_veCrds coordinates of the vertices.
##
def evalVe( i_nVars,
            i_expr,
            i_enVe,
            i_veCrds ):
  # determine result values
  l_vals = numpy.zeros( (i_nVars, len(i_enVe) ), dtype = 'float64')

  l_comp = compile( i_expr, '<string>', 'exec' )

  # iterate over the elements
  for l_el in range(len(i_enVe)):
    l_tmpVars = numpy.zeros( i_nVars )

    # determine values of vertices
    for l_ve in i_enVe[l_el]:
      l_vars = { 'x': i_veCrds[l_ve, 0],
                 'y': i_veCrds[l_ve, 1],
                 'z': i_veCrds[l_ve, 2] }

      exec( l_comp, l_vars )
      l_tmpVars += l_vars['q']

    # average vertex values for the elements
    l_vals[:, l_el] = l_tmpVars / len(i_enVe[l_el])
  
  return l_vals
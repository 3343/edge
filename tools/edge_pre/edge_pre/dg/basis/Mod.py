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
# Modifies a given basis.
##
import sympy

##
# Unifies the symbols in the given two bases.
#
# @param i_syms1 symbols used in first basis.
# @param i_basis1 first basis.
# @param i_syms2 symbols used in second basis.
# @param i_basis2 second basis.
# @return modified second symbols and basis.
##
def unify( i_syms1, i_basis1,
           i_syms2, i_basis2 ):
  # init new symbols
  l_syms  = []

  # check that string unique is not part of second syms
  for l_sy in i_syms2:
    assert( 'unique' not in str(l_sy) )

  # iterate over second basis symbols
  for l_s2 in i_syms2:
    # default: symbol is not unique
    l_un = False
    # continue until unique
    while( l_un == False ):
      # check for uniqueness
      l_un = True
      for l_s1 in i_syms1:
        if( l_s1 == l_s2 ):
          l_un = False
      
      # modify if still not unique
      if( l_un == False ): l_s2 = sympy.symbols( str(l_s2)+"_unique" )
    # attach unique symbol
    l_syms = l_syms + [l_s2]

  # generate modified basis
  l_basis = []
  for l_ba in i_basis2:
    # replace symbols
    for l_sy in range( len(i_syms2) ):
      l_ba = l_ba.subs( i_syms2[l_sy], l_syms[l_sy] )
    l_basis = l_basis + [l_ba]

  return l_syms, l_basis
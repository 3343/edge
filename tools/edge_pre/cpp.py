#!/usr/bin/env python
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
# Standalone tool which converts an EDGEpre CSV-array to C++, defining the data in the namespace edge::pre.
##
import argparse

##
# Replaces the given sub-strings in the input string.
#
# @param i_str input string.
# @param i_reps replacements.
##
def replace( i_str,
             i_reps ):
  l_str = i_str

  for l_re in i_reps:
    assert( len(l_re) == 2 )
    l_str = l_str.replace( l_re[0], l_re[1] )

  return l_str

##
# Converts the given CSV to C++.
#
# @param i_csvFile CSV file which gets converted.
# @param i_objSpace name space (below edge::pre) which will be used.
# @param i_objName name of the C++ matrix (will be appended by 'g_'). Additionally a raw pointer to the first entry and a size entry are created ('g_' + [...] + 'Raw', 'g_' + [...] + 'Size').
# @param i_objType type of the C++ matrix.
##
def convert( i_csvFile,
             i_objSpace,
             i_objName,
             i_objType ):
  # read the file contents
  with open( i_csvFile, 'r' ) as l_reader:
    # determine the dimensions
    l_dim = [0,0,0]
    for l_ro in l_reader:
      if l_dim[1] == 0:
        l_dim[0] = len( l_ro.split(',') )-1

      if l_ro != '\n' and l_dim[2] == 0:
        l_dim[1] = l_dim[1]+1

      if l_ro == '\n':
        l_dim[2] = l_dim[2]+1

    # check for dims 0 and 1
    assert( l_dim[0] > 0 )
    assert( l_dim[1] > 0 )

    # reset reader
    l_reader.seek(0)

    # generate code for initilization of the matrix
    l_init = '{\n'
    for l_d2 in range(max(l_dim[2],1)):
      if( l_dim[2] > 0 ): l_init = l_init + '  {\n'

      for l_d1 in range(l_dim[1]):
        l_init = l_init + '    {  '
        l_init = l_init + l_reader.readline()[:-1]
        l_init = l_init + '  },\n'

      if( l_dim[2] > 0 ):
        l_reader.readline()
        l_init = l_init + '  },\n'
    l_init = l_init + '}'

    # assemble C++ variable
    l_var = 'g_' + i_objName
    if( l_dim[2] > 0 ):
      l_var = l_var + '[' + str(l_dim[2]) + ']'
    l_var = l_var + '[' + str(l_dim[1]) + '][' + str(l_dim[0]) + ']'

    # assemble namespace
    l_ns = 'edge::pre::' + i_objSpace

    # assemble C++ header
    l_head = '#include <cstddef>\n'
    l_head = l_head + '#include <limits>\n'
    l_head = l_head + 'namespace edge {\n'
    l_head = l_head + '  namespace pre {\n'
    l_head = l_head + '    namespace ' + i_objSpace + ' {\n'
    l_head = l_head + '      extern ' + i_objType + ' const   '    + l_var + ';\n'
    l_head = l_head + '      extern ' + i_objType + ' const * g_' + i_objName + 'Raw;\n'
    l_head = l_head + '      extern std::size_t const g_' + i_objName + 'Size;\n'
    l_head = l_head + '    }\n'
    l_head = l_head + '  }\n'
    l_head = l_head + '}\n'

    # assemble C++ wrapper
    l_wrap = l_head + i_objType + ' const';
    l_wrap = l_wrap + ' ' + l_ns + '::'
    l_wrap = l_wrap + l_var + ' = '

    l_wrap = l_wrap + l_init + ';\n'
    l_wrap = l_wrap + i_objType + ' const * ' + l_ns + '::g_' + i_objName + 'Raw = &' + l_ns + '::g_' + i_objName + '[0]'
    if( l_dim[1] > 0 ): l_wrap = l_wrap + '[0]'
    if( l_dim[2] > 0 ): l_wrap = l_wrap + '[0]'
    l_wrap = l_wrap + ';\n'
    l_nEnt = l_dim[0] * l_dim[1] * max(l_dim[2], 1)
    l_wrap = l_wrap + 'std::size_t const ' + l_ns + '::g_' + i_objName + 'Size = ' + str(l_nEnt) + ';\n'

    return l_wrap

##
# Main
##
# parse command line options
l_parser = argparse.ArgumentParser( description='Converts EDGEpre\'s CSV-files to C++.' )

l_parser.add_argument( '-i', '--in_csv',
                       dest     = 'in_csv',
                       required = True,
                       type     = str,
                       help     = 'Input: CSV-file of EDGEpre')

l_parser.add_argument( '-s', '--obj_name_space',
                       dest     = 'obj_name_space',
                       required = True,
                       type     = str,
                       help     = 'C++ namespace of the object')

l_parser.add_argument( '-n', '--obj_name',
                       dest     = 'obj_name',
                       required = True,
                       type     = str,
                       help     = 'C++ name of the object')

l_parser.add_argument( '-t', '--obj_type',
                       dest     = 'obj_type',
                       required = True,
                       type     = str,
                       help     = 'C++ type of the object')

l_parser.add_argument( '-r', '--rep',
                       dest     = 'rep',
                       required = False,
                       type     = str,
                       nargs    = '*',
                       help     = 'Optional string replacements applied as last step to the C++-content. Expected is a list of tuples. For example, \'--reg "-1" abc "-2" def\' would replace all occurences of "-1" with "abc" and all of "-2" with "def".' )

l_parser.add_argument( '-o', '--out_cpp',
                       dest     = 'out_cpp',
                       required = True,
                       type     = str,
                       help     = 'Output: C++-file')

l_args = vars(l_parser.parse_args())

# derive replacements
l_reps = []
if l_args['rep']:
  for l_re in range(0, len(l_args['rep']), 2 ):
    l_reps = l_reps + [[ l_args['rep'][l_re], l_args['rep'][l_re+1] ]]

with open( l_args['out_cpp'], 'w' ) as l_fi:
  l_fi.write( replace(
                convert( l_args['in_csv'],
                         l_args['obj_name_space'],
                         l_args['obj_name'],
                         l_args['obj_type'] ),
                l_reps )
            )

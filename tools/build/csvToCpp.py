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
# Converts a matrix given as CSV to a C++ source file, defining the matrix in the namespace edge::pre.
##

##
# Converts the given CSV to C++.
#
# @param i_csvFile CSV file which gets converted.
# @oaram i_objSpace name space (below edge::pre) which will be used.
# @param i_objName name of the C++ matrix (will be appended by 'g_'). Additionally a raw pointer to the first entry and a size entry are created ('g_' + [...] + 'Raw', 'g_' + [...] + 'Size').
# @param i_objType type of the C++ matrix.
#
# @return C++ code as a string.
##
def convert( i_csvFile,
             i_objSpace,
             i_objName,
             i_objType ):
  # read the file contents
  with open( i_csvFile, 'rb' ) as l_reader:
    # generate code for initilization of the matrix
    l_init = '{\n'
    for l_ro in l_reader:
      l_init = l_init + '  {\n  '
      l_init = l_init + l_ro
      l_init = l_init + '  },\n'

    l_init = l_init + '}'

    # reset reader
    l_reader.seek(0)

    # get the number of rows and columns
    l_nRo =  0;
    l_nCo = -1;

    for l_ro in l_reader:
      l_nRo = l_nRo + 1
      # count number of columns and make sure all rows have the same number of columns
      l_nCoTmp = len( l_ro.split(',') )
      assert( l_nCo == -1 or l_nCo == l_nCoTmp )
      l_nCo = l_nCoTmp

    assert( l_nRo != 0 and l_nCo != 0 )

    # assemble C++ variable
    l_var = 'g_' + i_objName +  '[' + str(l_nRo) + '][' + str(l_nCo) + ']'

    # assemble namespace
    l_ns = 'edge::pre::' + i_objSpace

    # assemble C++ header
    l_head = '#include <cstddef>\n'
    l_head = l_head + 'namespace edge {\n'
    l_head = l_head + '  namespace pre {\n'
    l_head = l_head + '    namespace ' + i_objSpace + ' {\n'
    l_head = l_head + '      extern ' + i_objType + '        '    + l_var + ';\n'
    l_head = l_head + '      extern ' + i_objType + '      * g_' + i_objName + 'Raw;\n'
    l_head = l_head + '      extern std::size_t const   g_' + i_objName + 'Size;\n'
    l_head = l_head + '    }\n'
    l_head = l_head + '  }\n'
    l_head = l_head + '}\n'

    # assemble C++ wrapper
    l_wrap = l_head + i_objType;
    l_wrap = l_wrap + ' ' + l_ns + '::'
    l_wrap = l_wrap + l_var + ' = '

    l_wrap = l_wrap + l_init + ';\n'
    l_wrap = l_wrap + i_objType + ' * ' + l_ns + '::g_' + i_objName + 'Raw = &' + l_ns + '::g_' + i_objName + '[0][0];\n'
    l_wrap = l_wrap + 'std::size_t const ' + l_ns + '::g_' + i_objName + 'Size = ' + str(l_nRo*l_nCo) + ';\n'

    return l_wrap

##
# Generates the cpp-file from a csv-file in SCons.
# The generated cpp-file uses double precision.
#
# @param target SCons target files (cpp-files).
# @param source SSons sourve fiels (csv-files).
# @param env SCons environment, which defines a dict for key 'edge_pre', which contains for the matrix name for the absolute path of every csv-file as key.
##
def csvToCpp( target, source, env ):
  assert( len(target) == len(source) )

  for l_en in range( len(source) ):
    l_contents = convert( source[l_en].get_abspath(),
                          'subcell',
                          env['edge_pre'][source[l_en].get_abspath()],
                          'double const' )

    with open( target[l_en].get_abspath(), 'w+' ) as l_writer:
      l_writer.write( l_contents )

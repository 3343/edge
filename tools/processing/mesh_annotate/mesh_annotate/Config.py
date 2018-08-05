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
# Configuration of the mesh annotations.
##
import xml.etree.ElementTree
import logging

class Config:
  ##
  # Constructor, which parses the XML-config.
  #
  # @param i_xml xml config.
  ##
  def __init__( self, i_xml ):
    # parse the xml
    l_xml = xml.etree.ElementTree.parse( i_xml )

    # read element config
    l_conf = l_xml.getroot()

    self.m_vars = []
    for l_va in l_conf.find('vars').findall('var'):
      self.m_vars.append( l_va.text )

    # expression for the variables
    self.m_expr = l_conf.find('expr').text

    # output mesh
    self.m_meshes = {}
    for l_ty in ['in', 'out']:
      self.m_meshes[l_ty] = l_conf.find('meshes').find(l_ty).text

    self.log()

  ##
  # Prints the configuration.
  #
  def log( self ):
    logging.info( '  XML config:' )
    for l_va in range( len(self.m_vars) ):
      logging.info( '    var #' + str(l_va) + ': ' + self.m_vars[l_va] )

    logging.info('  expr:')
    for l_li in self.m_expr.split('\n'):
      if l_li != '':
        logging.info( '    ' + l_li )
    
    logging.info(  'meshes:')
    for l_ty in ['in', 'out']:
      logging.info( '    ' + l_ty + ': ' + self.m_meshes[l_ty] )
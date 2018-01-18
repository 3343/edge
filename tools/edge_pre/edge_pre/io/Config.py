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
# Configuration of the pre-processing.
##
import xml.etree.ElementTree
import logging

class Config:
  ##
  # Constructor, which parses the XML-config.
  #
  # @param i_pathToXml xml config of the study.
  ##
  def __init__(self, i_pathToXml):
    # parse the xml
    l_xml = xml.etree.ElementTree.parse( i_pathToXml )

    # read element config
    self.m_degs = []
    self.m_types = []
    self.m_out = {}
    l_conf = l_xml.getroot()

    # degrees
    for l_de in l_conf.find('degs').findall('deg'):
      self.m_degs.append( int(l_de.text) )

    # types
    for l_ty in l_conf.find('types').findall('type'):
      self.m_types.append( l_ty.text )

    # output directories
    l_out = l_conf.find('out_dirs')
    if l_out is not None:
      for l_ot in ['dg', 'subcell', 'plots']:
        if l_out.find(l_ot) is not None:
          self.m_out[l_ot] = l_out.find(l_ot).text

    self.print()

  ##
  # Prints the configuration.
  ##
  def print(self):
    logging.info( 'EDGEpre\'s config:' )
    logging.info( '  elements:'        )
    l_degs = ''
    for l_de in self.m_degs:
      l_degs = l_degs + ' ' + str(l_de)
    logging.info( '    degs: ' + l_degs )
    l_types = ''
    for l_ty in self.m_types:
      l_types = l_types + ' ' + l_ty
    logging.info( '    types:' + l_types )
    logging.info( '    out_dirs:' )
    for l_ot in ['dg', 'subcell', 'plots']:
      if l_ot in self.m_out:
        logging.info( '      ' + l_ot + ': ' + self.m_out[l_ot] )

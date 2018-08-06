#!/usr/bin/env python
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
# Annotates a mesh.
##
import argparse
import logging
import numpy
import mesh_annotate.Config
import mesh_annotate.Mesh
import mesh_annotate.Expression

# parse command line options
l_parser = argparse.ArgumentParser( description='Mesh annotations in EDGE.' )

l_parser.add_argument( '-x', '--xml',
                       dest     = 'xml',
                       required = True,
                       type     = str,
                       help     = 'XML configuration of the mesh annotation.')
l_args = vars( l_parser.parse_args() )

# setup logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# welcome our users
logging.info( 'performing the mesh annotation, please hold..' )

# parse xml config
logging.info('parsing XML')
l_config = mesh_annotate.Config.Config( l_args['xml'] )

# read mesh
logging.info( 'reading mesh' )
l_mesh = mesh_annotate.Mesh.Mesh( l_config.m_meshes['in'] )

if( l_mesh.m_nDis == 2 ):
  l_elTy = 'tria'
elif( l_mesh.m_nDis == 3 ):
  l_elTy = 'tet'
else: assert( False )

logging.info( 'mesh info:' )
logging.info( '  #vertices: '  + str( l_mesh.nEn( 'vertex' ) ) )
logging.info( '  #edges: '     + str( l_mesh.nEn( 'edge' )   ) )
logging.info( '  #triangles: ' + str( l_mesh.nEn( 'tria' )   ) )
if( l_elTy == 'tet' ):
  logging.info( '  #tetrahedrons: ' + str( l_mesh.nEn( 'tet' ) ) )

# get vertex coordinates
logging.info('querying mesh for vertex coordinates')
l_veCrds = l_mesh.veCrds()
logging.info('querying mesh for element to vertex adjacency')

# get vertices adjacent to the element
l_elVe =  l_mesh.enAd(2, 0)

# reshape the array according to the number of vertices per element
if( l_mesh.m_nDis == 2 ):
  assert( l_mesh.nEn('tria')*3 == len(l_elVe) )
  l_elVe = numpy.reshape( l_elVe, (-1, 3) )
elif( l_mesh.m_nDis == 3 ):
  assert( l_mesh.nEn('tet')*4 == len(l_elVe) )
  l_elVe = numpy.reshape( l_elVe, (-1, 4) )
else:
  assert( False )

# evaluate expression
logging.info( 'evaluating expression' )
if( l_config.m_type == 'vertices' ):
  l_vals = mesh_annotate.Expression.evalVe( len(l_config.m_vars),
                                            l_config.m_expr,
                                            l_elVe,
                                            l_veCrds )
elif( l_config.m_type == 'centroids' ):
  l_cens = l_mesh.cen( l_elTy )
  l_vals = mesh_annotate.Expression.evalPt( len(l_config.m_vars),
                                            l_config.m_expr,
                                            l_cens )
else: assert( False )

# wipe global ids tag
logging.info( 'wiping entity ids' )
l_mesh.deleteIds()

# store the data
logging.info( 'writing annotated mesh' )
for l_va in range(len(l_config.m_vars)):
  l_mesh.store( l_elTy,
                l_config.m_vars[l_va],
                l_vals[l_va, :].flatten() )

l_mesh.write( l_config.m_meshes['out'] )

logging.info( 'done' )
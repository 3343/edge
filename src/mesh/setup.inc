/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Includes of the mesh representations.
 **/

#if defined PP_T_MESH_REGULAR

///////////////////////////////////////////////////////
// TODO: Replace line, quad, hex8r analogue to tet4. //
///////////////////////////////////////////////////////
#if defined PP_T_ELEMENTS_LINE
edge::mesh::Regular l_mesh( edge::mesh::Regular::Line,
                            l_config.m_nElementsX,
                            l_config.m_sizeX );
#elif defined PP_T_ELEMENTS_QUAD4R
edge::mesh::Regular l_mesh( edge::mesh::Regular::Quadrilateral,
                            l_config.m_nElementsX, l_config.m_sizeX,
                            l_config.m_nElementsY, l_config.m_sizeY );
#elif defined PP_T_ELEMENTS_HEX8R
edge::mesh::Regular l_mesh( edge::mesh::Regular::Hexahedral,
                            l_config.m_nElementsX, l_config.m_sizeX,
                            l_config.m_nElementsY, l_config.m_sizeY,
                            l_config.m_nElementsZ, l_config.m_sizeZ );
#elif defined PP_T_ELEMENTS_TET4
unsigned int l_nXHex[3];
l_nXHex[0] = l_config.m_nElementsX;
l_nXHex[1] = l_config.m_nElementsY;
l_nXHex[2] = l_config.m_nElementsZ;
double l_corner[3] = {0,0,0};
double l_dX[3];
l_dX[0] = l_config.m_sizeX / l_config.m_nElementsX;
l_dX[1] = l_config.m_sizeY / l_config.m_nElementsY;
l_dX[2] = l_config.m_sizeZ / l_config.m_nElementsZ;

edge::mesh::regular::Tet l_mesh;
l_mesh.init( l_nXHex,
             edge::parallel::g_rank,
             edge::parallel::g_nRanks,
             l_corner,
             l_dX );
#else
assert( false );
#endif

#elif defined PP_T_MESH_UNSTRUCTURED
int *l_bndVals = NULL;
if( l_config.m_bndConId.size() > 0 ) l_bndVals =  &l_config.m_bndConId[0];

edge::mesh::Unstructured l_mesh( l_config.m_bndConId.size(), l_bndVals, l_config.m_periodic );
l_mesh.read( l_config.m_meshFileIn, l_config.m_meshOptRead );

if( l_config.m_meshFileOut != "" ) {
  l_mesh.write( l_config.m_meshFileOut.c_str() );
}
#else
#error mesh type not supported
#endif

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
 * Finalize for the advection equation.
 **/

double l_norms[N_CRUNS][3][1];

edge::io::ErrorNorms l_errorWriter( l_config.m_errorNormsType,
                                    l_config.m_errorNormsFile );

for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
  real_mesh l_sizeX, l_sizeY, l_sizeZ;
#if defined PP_T_MESH_REGULAR
  l_sizeX = l_config.m_sizeX;
  l_sizeY = l_config.m_sizeY;
  l_sizeZ = l_config.m_sizeZ;
#else
  l_sizeX = l_sizeY = l_sizeZ = 9.424777960769379715387930149838508652591508198125317462924;
#endif

  if( l_errorWriter.outEnabled() )
    edge::advection::setups::Convergence::getSineErrorNorms( l_basis,
                                                             l_cfr,
                                                             l_mesh.getInMap(),
                                                             l_mesh.getElLayout(),
                                                             l_internal.m_connect,
                                                             l_internal.m_vertexChars,
                                                             l_internal.m_elementModePrivate1,
                                                             l_norms[l_cfr],
                                                             0,
                                                             l_sizeX,
                                                             0,
                                                             l_sizeY,
                                                             0,
                                                             l_sizeZ,
                                                             l_cfr * 2.0,
                                                             0,
                                                             0,
                                                             4 / N_DIM );
}


l_errorWriter.write( l_norms );

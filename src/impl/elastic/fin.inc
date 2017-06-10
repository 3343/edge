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
 * Finalize for elastics.
 **/

// get the error norms
double l_errorNorms[N_CRUNS][3][N_QUANTITIES];

// create output
edge::io::ErrorNorms l_errorWriter( l_config.m_errorNormsType,
                                    l_config.m_errorNormsFile );

bool l_convergence = false;

for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
  if( l_internal.m_globalPrivate1[0][l_run].initialDofs == setupChars::InitialDofs::planeWaves ) {
    // we assume global time stepping by passing nEntsOwn below
    EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

    // get plane wave error if output is requested
    if( l_errorWriter.outEnabled() )
      edge::elastic::setups::Convergence::getPlaneErrorNorms( l_run,
                                                              l_basis,
                                                              l_mesh.getInMap(),
                                                              l_mesh.getElLayout(),
                                                              l_internal.m_connect,
                                                              l_internal.m_vertexChars,
                                                              l_internal.m_elementChars,
                                                              l_internal.m_elementModePrivate1,
                                                              l_errorNorms[l_run],
                                                              -50+l_run*5,
                                                              -50+l_run*5,
                                                              -50+l_run*5 );
    l_convergence = true;
  }
}

if( l_convergence ) {
  l_errorWriter.write( l_errorNorms );
}

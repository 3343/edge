/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2019, Alexander Breuer
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
 * Setup for the shallow water equations.
 **/

// get bathymetry and displacement model from mesh
EDGE_LOG_INFO << "    querying mesh for bathymetry and displacement data";
bool l_meshData = true;
std::size_t l_nBath = l_edgeV.nVas( "/edge_v/bath" );
if( l_nBath != l_edgeV.nEls() ) l_meshData = false;

std::size_t l_nDisp[N_CRUNS] = {0};
for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
  l_nDisp[l_cr] = l_edgeV.nVas( "edge_v/disp_" + std::to_string(l_cr) );
  if( l_nDisp[l_cr] != l_edgeV.nEls() ) l_meshData = false;
}

if( l_meshData ) {
  float * l_bath = new float[l_nBath];
  l_edgeV.getData( "/edge_v/bath",
                   l_bath );

  float * l_disp[N_CRUNS];
  for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
    l_disp[l_cr] = new float[l_nDisp[l_cr]];
    l_edgeV.getData( "/edge_v/disp_" + std::to_string(l_cr),
                     l_disp[l_cr] );
  }

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < l_edgeV.nEls(); l_el++ ) {
    l_internal.m_elementModeShared1[l_el][0][0] = l_bath[l_el];
    for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
      // add displacement to water column since we share the bathymetry
      l_internal.m_elementModePrivate1[l_el][0][0][l_cr] = std::max( -l_bath[l_el] + l_disp[l_cr][l_el], real_base(0) );
      l_internal.m_elementModePrivate1[l_el][1][0][l_cr] = 0;
      l_internal.m_elementModePrivate1[l_el][2][0][l_cr] = 0;
    }
  }

  delete[] l_bath;
  for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
    delete[] l_disp[l_cr];
  }
}

if( l_meshData ) {
  EDGE_LOG_INFO << "    successfully obtained data from mesh";
}
else {
  EDGE_LOG_INFO << "    failed obtaining data from mesh, continuing w/o";
}

// setup shared memory parallelization
l_shared.regWrkRgn( 0, 0, 0,
                    0,
                    l_edgeV.nFas(),
                    0 );

l_shared.regWrkRgn( 0, 1, 1,
                    0,
                    l_edgeV.nEls(),
                    0 );

// initialize Bathymetry and DOFs
if( !l_meshData ) {
  for( std::size_t l_el = 0; l_el < l_edgeV.nEls(); l_el++ )
    l_internal.m_elementModeShared1[l_el][0][0] = 0;

  edge::setups::InitialDofs<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS >::dg( std::size_t(0),
                   l_edgeV.nEls(),
                   l_config.m_initValsExprStrs,
                   l_basis,
                   l_internal.m_connect.elVe,
                   l_internal.m_vertexChars,
                   l_internal.m_elementModePrivate1 );
}

// get time steps
// Warning: Those are only valid if the allow time step doesn't lower significantly
edge::swe::solvers::FiniteVolume<
  T_SDISC.ELEMENT,
  N_CRUNS
>::getTimeStepStatistics( l_edgeV.nEls(),
                          l_internal.m_elementChars,
                          l_internal.m_elementModePrivate1,
                          l_dT[0],
                          l_dT[1],
                          l_dT[2] );

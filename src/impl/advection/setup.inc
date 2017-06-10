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
 * Setup for the advection equation.
 **/
#ifdef PP_USE_MPI
  // init mpi layout
  EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

  l_mpi.initLayout( l_enLayouts[2],
                    l_internal.m_elementModePrivate2[0][0][0],
                    N_ELEMENT_MODES*N_CRUNS*sizeof(t_elementModePrivate2),
                    0,
                    1 );
#endif

// setup shared memory parallelization
for( int_tg l_tg = 0; l_tg < l_enLayouts[2].timeGroups.size(); l_tg++ ) {
  // local inner-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg );

  // local send-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg );

  // neigh inner-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg );

  // neigh send-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg );
}

// initialize DOFs
for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
  CHECK( l_config.m_setups[l_cfr] == "sine");

  real_mesh l_sizeX, l_sizeY, l_sizeZ;
#if defined PP_T_MESH_REGULAR
  l_sizeX = l_config.m_sizeX;
  l_sizeY = l_config.m_sizeY;
  l_sizeZ = l_config.m_sizeZ;
#else
  l_sizeX = l_sizeY = l_sizeZ = 9.424777960769379715387930149838508652591508198125317462924;
#endif

  edge::advection::setups::Convergence::setSineWaves(  l_cfr,
                                                       l_basis,
                                                       l_internal.m_nElements,
                                                       l_internal.m_connect.elVe,
                                                       l_internal.m_vertexChars,
                                                       l_internal.m_elementChars,
                                                       l_internal.m_elementModePrivate1,
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

// initialize wave speeds
edge::advection::setups::Convergence::setConstantSpeed( l_internal.m_nElements,
                                                        l_internal.m_elementShared1,
                                                        1,
                                                        0.5,  // ignored for #dim==1
                                                        -1 ); // ignored for #dim!=3

#if PP_ORDER > 1
// set up star "matrices"
edge::advection::solvers::AderDg::setupStarM( l_internal.m_nElements,
                                              l_internal.m_vertexChars,
                                              l_internal.m_connect.elVe,
                                              l_internal.m_elementShared1,
                                              l_internal.m_elementShared2 );
#endif

// set up flux solvers
edge::advection::solvers::common::setupSolvers( l_internal.m_nElements,
                                                l_internal.m_nFaces,
                                                l_mesh.getInMap()->faMeDa,
                                                l_mesh.getInMap()->faDaMe,
                                                l_mesh.getInMap()->elMeDa,
                                                l_mesh.getInMap()->elDaMe,
                                                l_internal.m_connect.elVe,
                                                l_internal.m_connect.faEl,
                                                l_internal.m_connect.elFa,
                                                l_internal.m_vertexChars,
                                                l_internal.m_faceChars,
                                                l_internal.m_elementShared1,
                                                l_internal.m_elementShared3 );

// get time step stats
edge::advection::solvers::common::getTimeStepStatistics( l_internal.m_nElements,
                                                         l_internal.m_elementChars,
                                                         l_internal.m_elementShared1,
                                                         l_dT[0],
                                                         l_dT[1],
                                                         l_dT[2] );

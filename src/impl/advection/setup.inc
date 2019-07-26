/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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

// allocate flex data for time buffers and DOFs
{
  int_spType  l_spTypes[1] = { C_LTS_EL[t_ltsEl::EL_DOFS] };
  std::size_t l_spSizes[1] = { 1 };

  real_base (**l_raw)[N_ELEMENT_MODES][N_CRUNS];

  l_raw = l_dynMem.flex<
    real_base [N_ELEMENT_MODES][N_CRUNS] >( l_internal.m_nElements,
                                            1,
                                            1,
                                            l_spTypes,
                                            l_spSizes,
                                            l_internal.m_elementChars );

  l_internal.m_globalShared3[0] = l_raw;
  l_internal.m_globalShared3[1] = l_raw + l_internal.m_nElements;
}

#ifdef PP_USE_MPI
  // init mpi layout
  EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

  l_mpi.addDefault( l_enLayouts[2],
                    l_internal.m_globalShared3[0][0][0][0],
                    N_ELEMENT_MODES*N_CRUNS*sizeof(real_base),
                    0,
                    1 );
#endif

// setup shared memory parallelization
for( int_tg l_tg = 0; l_tg < l_enLayouts[2].timeGroups.size(); l_tg++ ) {
  int_spType l_spTypes[2] = { LIMIT, EXTREMA };

  // local, inner-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      1, l_spTypes, l_internal.m_elementChars );

  // local, send-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      1, l_spTypes, l_internal.m_elementChars );

  // neigh, inner-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      2, l_spTypes, l_internal.m_elementChars );

  // neigh, send-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      2, l_spTypes, l_internal.m_elementChars );
}

// initialize DOFs
{
  // DG
  edge::setups::InitialDofs<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS >::dg( 0,
                   l_internal.m_nElements,
                   l_config.m_initValsExprStrs,
                   l_basis,
                   l_internal.m_connect.elVe,
                   l_internal.m_vertexChars,
                   l_internal.m_elementModePrivate1 );
}

// initialize wave speeds
edge::advection::setups::Convergence::setConstantSpeed( l_internal.m_nElements,
                                                        l_internal.m_elementShared1,
                                                        1,
                                                        0.5,  // ignored for #dim==1
                                                        -1 ); // ignored for #dim!=3

#if PP_ORDER > 1
// initialize ADER-DG solver
edge::advection::solvers::AderDg< real_base,
                                  T_SDISC.ELEMENT,
                                  ORDER,
                                  ORDER,
                                  N_CRUNS > l_aderDg;
l_internal.m_globalShared4[0] = &l_aderDg;

// set up star "matrices"
edge::advection::solvers::AderDg<
    real_base,
    T_SDISC.ELEMENT,
    ORDER,
    ORDER,
    N_CRUNS >::setupStarM( l_internal.m_nElements,
                           l_internal.m_vertexChars,
                           l_internal.m_connect.elVe,
                           l_internal.m_elementShared1,
                           l_internal.m_elementShared2 );
#endif

// get time step stats
edge::advection::solvers::common::getTimeStepStatistics( l_internal.m_nElements,
                                                         l_internal.m_elementChars,
                                                         l_internal.m_elementShared1,
                                                         l_dT[0],
                                                         l_dT[1],
                                                         l_dT[2] );

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

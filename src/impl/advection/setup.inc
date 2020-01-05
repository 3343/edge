/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
  int_spType  l_spTypes[2] = { C_LTS_EL[t_ltsEl::EL_INT_LT], C_LTS_EL[t_ltsEl::EL_INT_GT] };
  std::size_t l_spSizes[2] = { 1,                            1                            };

  real_base (**l_raw)[N_ELEMENT_MODES][N_CRUNS];

  l_raw = l_dynMem.flex<
    real_base [N_ELEMENT_MODES][N_CRUNS] >( l_edgeV.nEls(),
                                            2,
                                            1,
                                            l_spTypes,
                                            l_spSizes,
                                            l_internal.m_elementChars );

  l_internal.m_globalShared3[0] = l_raw;
  l_internal.m_globalShared3[1] = l_raw +   l_edgeV.nEls();
  l_internal.m_globalShared3[2] = l_raw + 2*l_edgeV.nEls();
}

#ifdef PP_USE_MPI
  l_distributed.init( l_edgeV.nTgs(),
                      C_ENT[T_SDISC.ELEMENT].N_FACES,
                      l_edgeV.nEls(),
                      N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS*sizeof(real_base),
                      l_edgeV.getCommStruct(),
                      l_edgeV.getSendFa(),
                      l_edgeV.getSendEl(),
                      l_edgeV.getRecvFa(),
                      l_edgeV.getRecvEl(),
                      l_dynMem );
#endif

// setup shared memory parallelization
std::size_t l_firstDe = 0;

for( std::size_t l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
  // local, inner
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_firstDe,
                      l_edgeV.nTgElsIn()[l_tg],
                      l_edgeV.nTgs() - l_tg );

  // neigh, inner-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 2,
                      l_firstDe,
                      l_edgeV.nTgElsIn()[l_tg],
                      l_edgeV.nTgs() - l_tg );

  l_firstDe  += l_edgeV.nTgElsIn()[l_tg];
}

for( std::size_t l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
  // local, send
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_firstDe,
                      l_edgeV.nTgElsSe()[l_tg],
                      l_edgeV.nTgs()*2 - l_tg );

  // neigh, send-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_firstDe,
                      l_edgeV.nTgElsSe()[l_tg],
                      l_edgeV.nTgs()*2 - l_tg );

  l_firstDe  += l_edgeV.nTgElsSe()[l_tg];
}

// initialize DOFs
{
  // DG
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

// initialize wave speeds
edge::advection::setups::Convergence::setConstantSpeed( l_edgeV.nEls(),
                                                        l_internal.m_elementShared1,
                                                        1, 0.5, -1 );

// initialize ADER-DG solver
edge::advection::solvers::AderDg< real_base,
                                  T_SDISC.ELEMENT,
                                  ORDER,
                                  ORDER,
                                  N_CRUNS > l_aderDg( l_dynMem );
l_internal.m_globalShared4[0] = &l_aderDg;

// set up star "matrices"
edge::advection::solvers::AderDg<
    real_base,
    T_SDISC.ELEMENT,
    ORDER,
    ORDER,
    N_CRUNS >::setupStarM( l_edgeV.nEls(),
                           l_internal.m_vertexChars,
                           l_internal.m_connect.elVe,
                           l_internal.m_elementShared1,
                           l_internal.m_elementShared2 );

// get time step stats
edge::advection::solvers::common::getTimeStepStatistics( l_edgeV.nEls(),
                                                         l_internal.m_elementChars,
                                                         l_internal.m_elementShared1,
                                                         l_dT[0],
                                                         l_dT[1],
                                                         l_dT[2] );

// sync background parameters and set up flux solvers
{
  std::size_t l_nComm = l_edgeV.nCommElFa();

  t_bgPars *l_bgParsIn  = new t_bgPars[l_nComm];
  t_bgPars *l_bgParsOut = new t_bgPars[l_nComm];

  for( std::size_t l_co = 0; l_co < l_nComm; l_co++ ) {
    std::size_t l_el = l_edgeV.getSendEl()[l_co];
    l_bgParsOut[l_co] = l_internal.m_elementShared1[l_el][0];
  }

  l_distributed.syncData( 0,
                          sizeof(t_bgPars),
                          (unsigned char *) l_bgParsOut,
                          (unsigned char *) l_bgParsIn );
  delete[] l_bgParsOut;

  edge::advection::solvers::common::setupSolvers( l_edgeV.nElsIn(),
                                                  l_edgeV.nElsSe(),
                                                  l_edgeV.nFas(),
                                                  l_edgeV.nCommElFa(),
                                                  l_edgeV.getRecvFa(),
                                                  l_edgeV.getRecvEl(),
                                                  l_internal.m_connect.elVe,
                                                  l_internal.m_connect.faEl,
                                                  l_internal.m_connect.elFa,
                                                  l_internal.m_connect.elFaEl,
                                                  l_internal.m_vertexChars,
                                                  l_internal.m_faceChars,
                                                  l_internal.m_elementShared1,
                                                  l_bgParsIn,
                                                  l_internal.m_elementShared3 );
  delete[] l_bgParsIn;
}

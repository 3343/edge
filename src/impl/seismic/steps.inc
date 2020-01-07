/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Invidual steps of a time group for seismic settings.
 **/
if( i_step == 0 ) {
  // ADER-DG local
  real_base (** l_sendPtrs)[N_FACE_MODES][N_CRUNS] = (real_base (**) [N_FACE_MODES][N_CRUNS]) m_sendPtrs[getUpdatesSync()%m_nCommBuffers];
  m_internal.m_globalShared4[0][0].local( i_first,
                                          i_size,
                                          getUpdatesSync()%2 == 0,
                                          m_covSimTime,
                                          m_dt,
                                          i_enSp[0],
                                          m_internal.m_elementChars,
                                          m_internal.m_connect.vIdElFaEl,
                                          m_internal.m_elementModePrivate1,
#if (PP_N_RELAXATION_MECHANISMS > 0)
                                          m_internal.m_elementModePrivate2[0],
#else
                                          nullptr,
#endif
                                          m_internal.m_globalShared6,
                                          l_sendPtrs,
                                          io_recvs );
}
else if( i_step == 1 ) {
#ifdef PP_HAS_HDF5
  // apply point sources
  m_internal.m_globalShared3[0].apply( i_first,
                                       i_size,
                          (real_base) (m_covSimTime),
                          (real_base) (m_covSimTime+m_dt),
                                       m_internal.m_elementModePrivate1 );
#endif
}
else if( i_step == 2 ) {
  // ADER-DG: neigh contrib
  real_base (** l_recvPtrs)[N_FACE_MODES][N_CRUNS] = (real_base (**) [N_FACE_MODES][N_CRUNS]) m_recvPtrs[getUpdatesSync()%m_nCommBuffers];
  m_internal.m_globalShared4[0][0].neigh( i_first,
                                          i_size,
                                          getUpdatesSync()%2 == 0,
                                          m_internal.m_faceChars,
                                          m_internal.m_elementChars,
                                          m_internal.m_connect.elFa,
                                          m_internal.m_connect.elFaEl,
                                          m_internal.m_connect.fIdElFaEl,
                                          m_internal.m_connect.vIdElFaEl,
                                          m_internal.m_globalShared6,
                                          m_internal.m_elementModePrivate1,
#if (PP_N_RELAXATION_MECHANISMS > 0)
                                          m_internal.m_elementModePrivate2[0],
#else
                                          nullptr,
#endif
                                          l_recvPtrs
                                        );
}
else EDGE_LOG_FATAL << "step not supported in seismic implementation: " << i_step;
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
 * Invidual steps of a time group update for the elastic wave equations.
 **/

// get admissibility ids
unsigned short l_admId[4];
edge::sc::Steering::getAdmIds( m_updatesSync,
                               l_admId );

// get sc tDofs ids
unsigned short l_tDofsIds[2];
edge::sc::Steering::getDofsIds( m_updatesSync,
                                l_tDofsIds );

// get extrema ids
unsigned short l_extIds[2];
edge::sc::Steering::getExtIds( m_updatesSync,
                               l_extIds );

if( i_step == 0 ) {
  // ADER-DG local
  edge::elastic::solvers::AderDg<
    T_SDISC.ELEMENT,
    N_QUANTITIES,
    ORDER,
    ORDER,
    N_CRUNS
  >::local( i_first,
            i_size,
            m_covSimTime,
            m_dT,
            i_enSp[0],
            m_internal.m_elementChars,
            m_internal.m_globalShared1[0],
            m_internal.m_elementShared4,
            m_internal.m_elementShared2,
            m_internal.m_elementModePrivate1,
            m_internal.m_globalShared6,
            io_recvs,
            m_internal.m_mm );
}
else if ( i_step == 1 ) {
  edge::elastic::solvers::AderDg<
    T_SDISC.ELEMENT,
    N_QUANTITIES,
    ORDER,
    ORDER,
    N_CRUNS
  >::rupture( i_first,
              i_size,
              i_enSp[0],
  (real_base) m_covSimTime,
  (real_base) m_dT,
              m_internal.m_globalShared2[0].connect.scDgAd,
              m_internal.m_globalShared2[0].connect.liDoLiDu,
              m_internal.m_globalShared2[0].connect.liLp,
              m_internal.m_elementShared4,
              m_internal.m_globalShared7[0],
              m_internal.m_globalShared5[0].gl,
              m_internal.m_globalShared5[0].fa,
              m_internal.m_globalShared5[0].sf,
              m_internal.m_globalShared2[0].tDofs[ l_tDofsIds[0] ],
              m_internal.m_globalShared2[0].adm[ l_admId[1] ],
              m_internal.m_globalShared2[0].lock,
              io_recvsSf,
              m_internal.m_mm
            );
}
else if( i_step == 2 ) {
  // ADER-DG: neigh contrib
  edge::elastic::solvers::AderDg<
    T_SDISC.ELEMENT,
    N_QUANTITIES,
    ORDER,
    ORDER,
    N_CRUNS
  >::neigh( i_first,
            i_size,
            i_enSp[0],
            i_enSp[1],
            i_enSp[2],
            m_internal.m_globalShared1[0],
            m_internal.m_globalShared2[0].ops.scatter,
            m_internal.m_faceChars,
            m_internal.m_elementChars,
            m_internal.m_elementShared3,
            m_internal.m_globalShared2[0].connect.faSfSc,
            m_internal.m_globalShared2[0].connect.scDgAd,
            m_internal.m_globalShared2[0].connect.lpFaLp,
            m_internal.m_connect.elFa,
            m_internal.m_connect.elFaEl,
            m_internal.m_connect.fIdElFaEl,
            m_internal.m_connect.vIdElFaEl,
            m_internal.m_globalShared6,
            m_internal.m_elementModePrivate1,
            m_internal.m_globalShared2[0].adm[ l_admId[1] ],
            m_internal.m_globalShared2[0].ext[ l_extIds[0] ],
            m_internal.m_globalShared2[0].ext[ l_extIds[1] ],
            m_internal.m_globalShared2[0].tDofs[ l_tDofsIds[1] ],
            m_internal.m_mm
          );
}
else if( i_step == 3 ) {
#ifdef PP_HAS_HDF5
  // apply point sources
  m_internal.m_globalShared3[0].apply( i_first,
                                       i_size,
                          (real_base) (m_covSimTime),
                          (real_base) (m_covSimTime+m_dT),
                                       m_internal.m_elementModePrivate1 );
#endif
}
else if( i_step == 4 ) {
  // reset counters in first step since sync
  if( m_updatesSync != 0 ){}
  else {
    sc::Steering::resetLimSync( i_first,
                                i_size,
                                m_internal.m_globalShared2[0].connect.lpLi,
                                m_internal.m_globalShared2[0].limSync );
  }

  edge::sc::Limiter<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS,
    edge::elastic::sc::Boundary<
      T_SDISC.ELEMENT,
      ORDER,
      N_CRUNS > >::aPost( i_first,
                          i_size,
              (real_base) m_dT,
                          m_internal.m_globalShared1[0].mat,
                          m_internal.m_globalShared2[0].ops,
                          m_internal.m_globalShared2[0].connect,
                          m_internal.m_globalShared2[0].iBndSu,
                          m_internal.m_faceChars,
    (real_base (*)[C_ENT[T_SDISC.ELEMENT].N_FACES][N_QUANTITIES][N_QUANTITIES]) m_internal.m_elementShared2, // TODO: fix type
    (real_base (*)[C_ENT[T_SDISC.ELEMENT].N_FACES][N_QUANTITIES][N_QUANTITIES]) m_internal.m_elementShared3, // TODO: fix type
                          m_internal.m_connect.elFa,
                          m_internal.m_connect.elFaEl,
                          m_internal.m_connect.vIdElFaEl,
                          m_internal.m_connect.fIdElFaEl,
                          m_internal.m_globalShared2[0].adm[ l_admId[0] ],
                          m_internal.m_globalShared2[0].adm[ l_admId[1] ],
                          m_internal.m_globalShared2[0].adm[ l_admId[2] ],
                          m_internal.m_globalShared2[0].adm[ l_admId[3] ],
                          m_internal.m_globalShared2[0].lock,
                          m_internal.m_globalShared2[0].limSync,
                          m_internal.m_globalShared6,
                          m_internal.m_elementModePrivate1,
                          m_internal.m_globalShared2[0].tDofs[ l_tDofsIds[0] ],
                          m_internal.m_globalShared2[0].tDofs[ l_tDofsIds[1] ],
                          m_internal.m_globalShared2[0].dofs,
                          m_internal.m_globalShared2[0].ext[ l_extIds[0] ],
                          m_internal.m_globalShared2[0].ext[ l_extIds[1] ],
                          m_internal.m_mm,
                          m_internal.m_globalShared8[0] );
}
else EDGE_LOG_FATAL << "step not supported in elastic implementation: " << i_step;

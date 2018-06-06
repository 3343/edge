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
 * Steps for the advection equation.
 **/
// assemble admissibility ids through permutations
#if PP_ORDER > 1
unsigned short l_admId[3];
if( (m_updatesPer%3) == 0 ) {
  l_admId[0] = 0;
  l_admId[1] = 1;
  l_admId[2] = 2;
}
else if( (m_updatesPer%3) == 1 ) {
  l_admId[0] = 2;
  l_admId[1] = 0;
  l_admId[2] = 1;
}
else {
  l_admId[0] = 1;
  l_admId[1] = 2;
  l_admId[2] = 0;
}
#endif

if( i_step == 0 ) {
#if PP_ORDER == 1
  edge::advection::solvers::FiniteVolume::timeInt( i_first,
                                                   i_size,
                                                   m_dT,
                                                   m_internal.m_elementModePrivate1,
                                                   m_internal.m_elementModePrivate2 );
#else
  edge::advection::solvers::AderDg<
    T_SDISC.ELEMENT,
    ORDER,
    ORDER,
    N_CRUNS >::local(  i_first,
                       i_size,
                       i_enSp[0].first,
                       m_dT,
                       m_internal.m_globalShared3[0].connect.faSfSc,
                       m_internal.m_globalShared1[0],
                       m_internal.m_elementChars,
                       m_internal.m_elementShared2,
                       m_internal.m_elementShared3,
                       m_internal.m_elementModePrivate1,
                       m_internal.m_globalShared3[0].dofs,
                       m_internal.m_globalShared2,
                       m_internal.m_globalShared3[0].surf );
#endif
}
else if( i_step == 1 ) {
#if PP_ORDER == 1
  edge::advection::solvers::FiniteVolume::update( i_first,
                                                  i_size,
                                                  m_internal.m_connect.elFaEl,
                                                  m_internal.m_elementShared3,
                                                  m_internal.m_elementModePrivate2,
                                                  m_internal.m_elementModePrivate1 );
#else
  edge::advection::solvers::AderDg<
    T_SDISC.ELEMENT,
    ORDER,
    ORDER,
    N_CRUNS >::neigh( i_first,
                      i_size,
                      i_enSp[0].first,
                      i_enSp[1].first,
                      m_internal.m_globalShared1[0],
                      m_internal.m_globalShared3[0].ops.scatter,
                      m_internal.m_faceChars,
                      m_internal.m_elementChars,
                      m_internal.m_elementShared3,
                      m_internal.m_globalShared3[0].connect.liVeEx,
                      m_internal.m_connect.elFa,
                      m_internal.m_connect.elFaEl,
                      m_internal.m_connect.fIdElFaEl,
                      m_internal.m_connect.vIdElFaEl,
                      m_internal.m_globalShared2,
                      m_internal.m_elementModePrivate1,
                      m_internal.m_globalShared3[0].adm[ l_admId[0] ],
                      m_internal.m_globalShared3[0].adm[ l_admId[1] ],
                      m_internal.m_globalShared3[0].ext[m_updatesPer%2],
                      m_internal.m_globalShared3[0].ext[(m_updatesPer+1)%2] );
#endif
}
else if( i_step == 2 ) {
#if PP_ORDER > 1
  edge::sc::Limiter<
    T_SDISC.ELEMENT,
    ORDER,
    N_CRUNS >::aPost( i_first,
                      i_size,
                      m_dT,
                      m_internal.m_globalShared1[0].mat.flux,
                      m_internal.m_globalShared3[0].ops,
                      m_internal.m_globalShared3[0].connect,
                      m_internal.m_elementShared3,
                      m_internal.m_connect.elFaEl,
                      m_internal.m_connect.vIdElFaEl,
                      m_internal.m_connect.fIdElFaEl,
                      m_internal.m_globalShared3[0].adm[ l_admId[0] ],
                      m_internal.m_globalShared3[0].adm[ l_admId[1] ],
                      m_internal.m_globalShared3[0].adm[ l_admId[2] ],
                      m_internal.m_globalShared2,
                      m_internal.m_elementModePrivate1,
                      m_internal.m_globalShared3[0].surf,
                      m_internal.m_globalShared3[0].dofs,
                      m_internal.m_globalShared3[0].ext[m_updatesPer%2],
                      m_internal.m_globalShared3[0].ext[(m_updatesPer+1)%2] );
#endif
}
else EDGE_LOG_FATAL << "step not supported in advection implementation: " << i_step;

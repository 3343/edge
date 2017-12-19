/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, Regents of the University of California
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
 * LTS time group.
 **/

#include "TimeGroupStatic.h"
#include "monitor/instrument.hpp"

#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/ts_dep.inc"
#elif defined PP_T_EQUATIONS_ELASTIC
#include "impl/elastic/ts_dep.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/ts_dep.inc"
#endif

edge::time::TimeGroupStatic::TimeGroupStatic(       int_ts          i_rate,
                                                    int_ts          i_funMult,
                                              const data::Internal &i_internal ):
// m_rate(     i_rate     ),
 m_funMult(  i_funMult  ),
 m_internal( i_internal )
{
  m_covSimTime = 0;
  m_updatesPer = 0;
  m_updatesReq = 0;
}

void edge::time::TimeGroupStatic::setUp( double i_dTfun,
                                         double i_time ) {
  // set general time step
  m_dTgen = i_dTfun * m_funMult;

  // reset number of updates since sync
  m_updatesSync = 0;

  // derive number of required updates
  m_updatesReq  = i_time / m_dTgen;
  m_updatesReq += 1;

  // derive final time step
  m_dTfin = std::max( 0.0, i_time - ( m_dTgen * (m_updatesReq-1) ) );

  // set time step of first update
  setDt();
}

void edge::time::TimeGroupStatic::setDt() {
  if( m_updatesReq > 1 ) m_dT = m_dTgen;
  else                   m_dT = m_dTfin;
}

void edge::time::TimeGroupStatic::updateTsInfo() {
  // update counters
  m_updatesSync++;
  m_updatesPer++;
  m_updatesReq--;

  m_covSimTime += m_dT;

  setDt();
}

void edge::time::TimeGroupStatic::computeStep( unsigned short                              i_step,
                                               int_el                                      i_first,
                                               int_el                                      i_size,
                                               t_timeRegion                        const * i_enSp,
                                               io::Receivers                             & io_recvs,
                                               io::ReceiversQuad< real_base,
                                                                  T_SDISC.ELEMENT,
                                                                  ORDER,
                                                                  N_CRUNS >              & io_recvsQuad  ) {
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/inc/time/tgs_steps.inc"
#elif defined PP_T_EQUATIONS_ELASTIC
#include "impl/elastic/steps.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/inc/time/tgs_steps.inc"
#else
#error "steps not defined"
#endif
}

void edge::time::TimeGroupStatic::limSync() {
  sc::Steering::resetAdm( m_updatesSync,
                          m_internal.m_globalShared2[0].adm );

  sc::Steering::resetExt( m_updatesSync,
                          m_internal.m_globalShared2[0].ext );
}
/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2015-2018, Regents of the University of California
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
#elif defined PP_T_EQUATIONS_SEISMIC
#include "impl/seismic/ts_dep.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/ts_dep.inc"
#endif

#include "sc/Steering.hpp"
#include "sc/Limiter.hpp"

edge::time::TimeGroupStatic::TimeGroupStatic( unsigned short    i_nTgs,
                                              unsigned short    i_tgId,
                                              data::Internal & i_internal ): m_internal( i_internal ) {
  // derive multiple of fundamental time step (rate-2 LTS)
  m_funMul = 1;
  for( unsigned short l_tg = 0; l_tg < i_tgId; l_tg++ ) {
    m_funMul *= 2;
  }
  // derive divisor of max time step
  m_maxDiv = 1;
  for( unsigned short l_tg = i_tgId+1; l_tg < i_nTgs; l_tg++ ) {
    m_maxDiv *= 2;
  }

  m_covSimTime = 0;
  m_nTsSync = 0;
  m_nTsPer = 0;
  m_nTsReqFull = 0;
  m_dt     = std::numeric_limits< double >::max();
  m_dtFull = std::numeric_limits< double >::max();
  m_dtSync = std::numeric_limits< double >::max();
}

void edge::time::TimeGroupStatic::setUp( double i_dtFun,
                                         double i_time ) {
  // set full time step of the group
  m_dtFull = i_dtFun * m_funMul;

  // reset number of updates since sync
  m_nTsSync = 0;

  // derive number of full time steps
  double l_dtMax = i_dtFun * (m_funMul * m_maxDiv);
  std::size_t l_nUpdatesMaxFull = i_time / l_dtMax;
  m_nTsReqFull = l_nUpdatesMaxFull * m_maxDiv;

  // set sync time step
  m_dtSync = i_time - l_nUpdatesMaxFull * l_dtMax;
  m_dtSync /= m_maxDiv;
  m_dtSync = std::max( 0.0, m_dtSync );

  // set time step of first update
  setDt();
}

void edge::time::TimeGroupStatic::setDt() {
  if( m_nTsSync < m_nTsReqFull ) {
    m_dt = m_dtFull;
  }
  else {
    m_dt = m_dtSync;
  }
}

void edge::time::TimeGroupStatic::updateTsInfo() {
  // update counters
  m_nTsSync++;
  m_nTsPer++;

  m_covSimTime += m_dt;

  setDt();
}

void edge::time::TimeGroupStatic::computeStep( unsigned short                              i_step,
                                               int_el                                      i_first,
                                               int_el                                      i_size,
                                               int_el                              const * i_enSp,
                                               io::Receivers                             & io_recvs,
                                               io::ReceiversSf< real_base,
                                                                T_SDISC.ELEMENT,
                                                                ORDER,
                                                                N_CRUNS >                & io_recvsSf  ) {
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/inc/time/tgs_steps.inc"
#elif defined PP_T_EQUATIONS_SEISMIC
#include "impl/seismic/steps.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/inc/time/tgs_steps.inc"
#else
#error "steps not defined"
#endif
}

void edge::time::TimeGroupStatic::limSync() {
  sc::Steering::resetAdm( m_nTsSync,
                          m_internal.m_globalShared2[0].adm );

  sc::Steering::resetDofs( m_nTsSync,
                           m_internal.m_globalShared2[0].tDofs );

  sc::Steering::resetExt( m_nTsSync,
                          m_internal.m_globalShared2[0].ext );
}
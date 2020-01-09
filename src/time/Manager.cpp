/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Management of the time stepping.
 **/

#include "Manager.h"
#include "monitor/instrument.hpp"

void edge::time::Manager::schedule() {
#if defined PP_T_EQUATIONS_ADVECTION
#include "src/impl/advection/inc/time/man_sched.inc"
#elif defined PP_T_EQUATIONS_SEISMIC
#include "src/impl/seismic/inc/time/man_sched.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "src/impl/swe/inc/time/man_sched.inc"
#else
#error "scheduling not defined"
#endif
}

edge::time::Manager::Manager( double                            i_dt,
                              parallel::Shared                & i_shared,
                              parallel::Distributed           & i_distributed,
                              std::vector< TimeGroupStatic  > & i_timeGroups,
                              io::Receivers                   & i_recvs ): m_dTfun(       i_dt          ),
                                                                           m_shared(      i_shared      ),
                                                                           m_distributed( i_distributed ),
                                                                           m_recvs(       i_recvs       ) {
  for( std::size_t l_tg = 0; l_tg < i_timeGroups.size(); l_tg++ ) {
    m_timeGroups.push_back( &i_timeGroups[l_tg] );
  }
  m_cflow = ( unsigned short(*)[N_ENTRIES_CONTROL_FLOW] ) new unsigned short[i_timeGroups.size() * N_ENTRIES_CONTROL_FLOW];
}

edge::time::Manager::~Manager() {
  delete[] m_cflow;
}

bool edge::time::Manager::getTimePredAvailable( unsigned short i_tg ) {
  // time group with smaller time step
  bool l_left = (i_tg == 0) ?
                true :
                   m_timeGroups[i_tg-1]->nTimePredInner() == m_timeGroups[i_tg]->nTimePredInner()*2
                && m_timeGroups[i_tg-1]->nTimePredSend()  == m_timeGroups[i_tg]->nTimePredSend()*2;

  // time group with larger time step
  bool l_right = (i_tg == m_timeGroups.size() - 1 ) ?
                 true :
                      m_timeGroups[i_tg]->nTimePredInner() <= m_timeGroups[i_tg+1]->nTimePredInner()*2
                   && m_timeGroups[i_tg]->nTimePredSend()  <= m_timeGroups[i_tg+1]->nTimePredSend()*2;

  return l_left && l_right;
}

bool edge::time::Manager::getTimePredConsumed( unsigned short i_tg ) {
  // time group with smaller time step
  bool l_left = (i_tg == 0) ?
                true :
                      m_timeGroups[i_tg-1]->nDofUpInner() == m_timeGroups[i_tg]->nTimePredInner()*2
                   && m_timeGroups[i_tg-1]->nDofUpSend()  == m_timeGroups[i_tg]->nTimePredSend()*2;

  // time group with larger time step: last time group, odd time step, even time step (accumulated data consumed)
  bool l_right = (i_tg == m_timeGroups.size() - 1 ) ?
                 true :    m_timeGroups[i_tg]->getUpdatesSync()%2 == 0
                        || (    m_timeGroups[i_tg]->nTimePredInner() == m_timeGroups[i_tg+1]->nDofUpInner()*2
                             && m_timeGroups[i_tg]->nTimePredSend()  == m_timeGroups[i_tg+1]->nDofUpSend()*2 );

  return l_left && l_right;
}

void edge::time::Manager::communicate() {
  m_distributed.comm();
}

void edge::time::Manager::compute() {
  PP_INSTR_FUN("compute")

  // scheduling and communicating workers have other duties, pure workers stay where they are
  bool l_schdCmm = m_shared.isSched() || m_shared.isComm();

  while( m_finished == false ) {
    bool l_wrk;
    unsigned short l_tg;
    unsigned short l_st;
    unsigned int   l_id;
    int_el         l_first;
    int_el         l_size;
    int_el         l_enSp[128];

    // check for work
    l_wrk = m_shared.getWrkTd( l_tg, l_st, l_id, l_first, l_size, l_enSp );

    if( l_wrk == true ) {
      // set status to "in progress"
      m_shared.setStatusTd( parallel::Shared::IPR, l_id );

      // compute
      PP_INSTR_REG_DEF(step)
      PP_INSTR_REG_BEG(step,"step")
#pragma warning push
#pragma warning(disable:68)
      PP_INSTR_PAR_UINT64("step_id",  (uint64_t) l_st )
      PP_INSTR_PAR_UINT64("cflow_id", (uint64_t) l_id )
#pragma warning pop

      m_timeGroups[l_tg]->computeStep( l_st, l_first, l_size, l_enSp, m_recvs );

      PP_INSTR_REG_END(step)

      // set status to "finished"
      m_shared.setStatusTd( parallel::Shared::FIN, l_id );
    }

    // non-pure workers are allowed to exit
    if( l_schdCmm == true ) break;
  }
}

void edge::time::Manager::simulate( double i_time ) {
  PP_INSTR_FUN("simulate")

  // propagate sync time to all time groups
  for( unsigned short l_tg = 0; l_tg < m_timeGroups.size(); l_tg++ ) {
    m_timeGroups[l_tg]->setUp( m_dTfun,
                               i_time );
  }

  // reset all statuses to wait
  m_shared.resetStatus( parallel::Shared::WAI );

  // reset the distributed memory interface
  m_distributed.reset();

  // reset control flow
  for( unsigned short l_tg = 0; l_tg < m_timeGroups.size(); l_tg++ ) {
    for( unsigned short l_cf = 0; l_cf < N_ENTRIES_CONTROL_FLOW; l_cf++ ) {
      m_cflow[l_tg][l_cf] = std::numeric_limits< unsigned short >::max();
    }
  }

  // we are not finished until the scheduling threads decides so
  m_finished = false;

  // jump into respective tasks
#ifdef PP_USE_OMP
#pragma omp parallel
{
#endif
  while( m_finished == false ) {
    if( m_shared.isSched() ) schedule();
    if( m_shared.isComm()  ) communicate();
    if( m_shared.isWrk()   ) compute();
  }
#ifdef PP_USE_OMP
}
#endif

  // (re-)balance work regions
  m_shared.balance();
}

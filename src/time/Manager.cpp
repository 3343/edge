/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
#ifdef PP_USE_INSTR
      PP_INSTR_REG_DEF( wrk_rgn_0  ) PP_INSTR_REG_DEF( wrk_rgn_1  )
      PP_INSTR_REG_DEF( wrk_rgn_2  ) PP_INSTR_REG_DEF( wrk_rgn_3  )
      PP_INSTR_REG_DEF( wrk_rgn_4  ) PP_INSTR_REG_DEF( wrk_rgn_5  )
      PP_INSTR_REG_DEF( wrk_rgn_6  ) PP_INSTR_REG_DEF( wrk_rgn_7  )
      PP_INSTR_REG_DEF( wrk_rgn_8  ) PP_INSTR_REG_DEF( wrk_rgn_9  )
      PP_INSTR_REG_DEF( wrk_rgn_10 ) PP_INSTR_REG_DEF( wrk_rgn_11 )
      PP_INSTR_REG_DEF( wrk_rgn_12 ) PP_INSTR_REG_DEF( wrk_rgn_13 )
      PP_INSTR_REG_DEF( wrk_rgn_14 ) PP_INSTR_REG_DEF( wrk_rgn_15 )
      PP_INSTR_REG_DEF( wrk_rgn_16 ) PP_INSTR_REG_DEF( wrk_rgn_17 )
      PP_INSTR_REG_DEF( wrk_rgn_18 ) PP_INSTR_REG_DEF( wrk_rgn_19 )
      PP_INSTR_REG_DEF( wrk_rgn_20 ) PP_INSTR_REG_DEF( wrk_rgn_21 )
      PP_INSTR_REG_DEF( wrk_rgn_22 ) PP_INSTR_REG_DEF( wrk_rgn_23 )
      PP_INSTR_REG_DEF( wrk_rgn_24 ) PP_INSTR_REG_DEF( wrk_rgn_25 )
      PP_INSTR_REG_DEF( wrk_rgn_26 ) PP_INSTR_REG_DEF( wrk_rgn_27 )
      PP_INSTR_REG_DEF( wrk_rgn_28 ) PP_INSTR_REG_DEF( wrk_rgn_29 )
      PP_INSTR_REG_DEF( wrk_rgn_30 ) PP_INSTR_REG_DEF( wrk_rgn_31 )
      PP_INSTR_REG_DEF( wrk_rgn_32 ) PP_INSTR_REG_DEF( wrk_rgn_33 )
      PP_INSTR_REG_DEF( wrk_rgn_34 ) PP_INSTR_REG_DEF( wrk_rgn_35 )
      PP_INSTR_REG_DEF( wrk_rgn_36 ) PP_INSTR_REG_DEF( wrk_rgn_37 )
      PP_INSTR_REG_DEF( wrk_rgn_38 ) PP_INSTR_REG_DEF( wrk_rgn_39 )
      PP_INSTR_REG_DEF( wrk_rgn_undef )
      if(      l_id ==  0  ) { PP_INSTR_REG_BEG( wrk_rgn_0,  "wrk_rgn_0"  ) } else if( l_id ==  1  ) { PP_INSTR_REG_BEG( wrk_rgn_1,  "wrk_rgn_1"  ) }
      else if( l_id ==  2  ) { PP_INSTR_REG_BEG( wrk_rgn_2,  "wrk_rgn_2"  ) } else if( l_id ==  3  ) { PP_INSTR_REG_BEG( wrk_rgn_3,  "wrk_rgn_3"  ) }
      else if( l_id ==  4  ) { PP_INSTR_REG_BEG( wrk_rgn_4,  "wrk_rgn_4"  ) } else if( l_id ==  5  ) { PP_INSTR_REG_BEG( wrk_rgn_5,  "wrk_rgn_5"  ) }
      else if( l_id ==  6  ) { PP_INSTR_REG_BEG( wrk_rgn_6,  "wrk_rgn_6"  ) } else if( l_id ==  7  ) { PP_INSTR_REG_BEG( wrk_rgn_7,  "wrk_rgn_7"  ) }
      else if( l_id ==  8  ) { PP_INSTR_REG_BEG( wrk_rgn_8,  "wrk_rgn_8"  ) } else if( l_id ==  9  ) { PP_INSTR_REG_BEG( wrk_rgn_9,  "wrk_rgn_9"  ) }
      else if( l_id == 10  ) { PP_INSTR_REG_BEG( wrk_rgn_10, "wrk_rgn_10" ) } else if( l_id == 11  ) { PP_INSTR_REG_BEG( wrk_rgn_11, "wrk_rgn_11" ) }
      else if( l_id == 12  ) { PP_INSTR_REG_BEG( wrk_rgn_12, "wrk_rgn_12" ) } else if( l_id == 13  ) { PP_INSTR_REG_BEG( wrk_rgn_13, "wrk_rgn_13" ) }
      else if( l_id == 14  ) { PP_INSTR_REG_BEG( wrk_rgn_14, "wrk_rgn_14" ) } else if( l_id == 15  ) { PP_INSTR_REG_BEG( wrk_rgn_15, "wrk_rgn_15" ) }
      else if( l_id == 16  ) { PP_INSTR_REG_BEG( wrk_rgn_16, "wrk_rgn_16" ) } else if( l_id == 17  ) { PP_INSTR_REG_BEG( wrk_rgn_17, "wrk_rgn_17" ) }
      else if( l_id == 18  ) { PP_INSTR_REG_BEG( wrk_rgn_18, "wrk_rgn_18" ) } else if( l_id == 19  ) { PP_INSTR_REG_BEG( wrk_rgn_19, "wrk_rgn_19" ) }
      else if( l_id == 20  ) { PP_INSTR_REG_BEG( wrk_rgn_20, "wrk_rgn_20" ) } else if( l_id == 21  ) { PP_INSTR_REG_BEG( wrk_rgn_21, "wrk_rgn_21" ) }
      else if( l_id == 22  ) { PP_INSTR_REG_BEG( wrk_rgn_22, "wrk_rgn_22" ) } else if( l_id == 23  ) { PP_INSTR_REG_BEG( wrk_rgn_23, "wrk_rgn_23" ) }
      else if( l_id == 24  ) { PP_INSTR_REG_BEG( wrk_rgn_24, "wrk_rgn_24" ) } else if( l_id == 25  ) { PP_INSTR_REG_BEG( wrk_rgn_25, "wrk_rgn_25" ) }
      else if( l_id == 26  ) { PP_INSTR_REG_BEG( wrk_rgn_26, "wrk_rgn_26" ) } else if( l_id == 27  ) { PP_INSTR_REG_BEG( wrk_rgn_27, "wrk_rgn_27" ) }
      else if( l_id == 28  ) { PP_INSTR_REG_BEG( wrk_rgn_28, "wrk_rgn_28" ) } else if( l_id == 29  ) { PP_INSTR_REG_BEG( wrk_rgn_29, "wrk_rgn_29" ) }
      else if( l_id == 30  ) { PP_INSTR_REG_BEG( wrk_rgn_30, "wrk_rgn_30" ) } else if( l_id == 31  ) { PP_INSTR_REG_BEG( wrk_rgn_31, "wrk_rgn_31" ) }
      else if( l_id == 32  ) { PP_INSTR_REG_BEG( wrk_rgn_32, "wrk_rgn_32" ) } else if( l_id == 33  ) { PP_INSTR_REG_BEG( wrk_rgn_33, "wrk_rgn_33" ) }
      else if( l_id == 34  ) { PP_INSTR_REG_BEG( wrk_rgn_34, "wrk_rgn_34" ) } else if( l_id == 35  ) { PP_INSTR_REG_BEG( wrk_rgn_35, "wrk_rgn_35" ) }
      else if( l_id == 36  ) { PP_INSTR_REG_BEG( wrk_rgn_36, "wrk_rgn_36" ) } else if( l_id == 37  ) { PP_INSTR_REG_BEG( wrk_rgn_37, "wrk_rgn_37" ) }
      else if( l_id == 38  ) { PP_INSTR_REG_BEG( wrk_rgn_38, "wrk_rgn_38" ) } else if( l_id == 39  ) { PP_INSTR_REG_BEG( wrk_rgn_39, "wrk_rgn_39" ) }
      else { PP_INSTR_REG_BEG( wrk_rgn_undef, "wrk_rgn_undef" ) }
#endif

      m_timeGroups[l_tg]->computeStep( l_st, l_first, l_size, l_enSp, m_recvs );

#ifdef PP_USE_INSTR
      if(      l_id ==  0  ) { PP_INSTR_REG_END( wrk_rgn_0  ) } else if( l_id ==   1  ) { PP_INSTR_REG_END( wrk_rgn_1  ) }
      else if( l_id ==  2  ) { PP_INSTR_REG_END( wrk_rgn_2  ) } else if( l_id ==   3  ) { PP_INSTR_REG_END( wrk_rgn_3  ) }
      else if( l_id ==  4  ) { PP_INSTR_REG_END( wrk_rgn_4  ) } else if( l_id ==   5  ) { PP_INSTR_REG_END( wrk_rgn_5  ) }
      else if( l_id ==  6  ) { PP_INSTR_REG_END( wrk_rgn_6  ) } else if( l_id ==   7  ) { PP_INSTR_REG_END( wrk_rgn_7  ) }
      else if( l_id ==  8  ) { PP_INSTR_REG_END( wrk_rgn_8  ) } else if( l_id ==   9  ) { PP_INSTR_REG_END( wrk_rgn_9  ) }
      else if( l_id == 10  ) { PP_INSTR_REG_END( wrk_rgn_10 ) } else if( l_id ==  11  ) { PP_INSTR_REG_END( wrk_rgn_11 ) }
      else if( l_id == 12  ) { PP_INSTR_REG_END( wrk_rgn_12 ) } else if( l_id ==  13  ) { PP_INSTR_REG_END( wrk_rgn_13 ) }
      else if( l_id == 14  ) { PP_INSTR_REG_END( wrk_rgn_14 ) } else if( l_id ==  15  ) { PP_INSTR_REG_END( wrk_rgn_15 ) }
      else if( l_id == 16  ) { PP_INSTR_REG_END( wrk_rgn_16 ) } else if( l_id ==  17  ) { PP_INSTR_REG_END( wrk_rgn_17 ) }
      else if( l_id == 18  ) { PP_INSTR_REG_END( wrk_rgn_18 ) } else if( l_id ==  19  ) { PP_INSTR_REG_END( wrk_rgn_19 ) }
      else if( l_id == 20  ) { PP_INSTR_REG_END( wrk_rgn_20 ) } else if( l_id ==  21  ) { PP_INSTR_REG_END( wrk_rgn_21 ) }
      else if( l_id == 22  ) { PP_INSTR_REG_END( wrk_rgn_22 ) } else if( l_id ==  23  ) { PP_INSTR_REG_END( wrk_rgn_23 ) }
      else if( l_id == 24  ) { PP_INSTR_REG_END( wrk_rgn_24 ) } else if( l_id ==  25  ) { PP_INSTR_REG_END( wrk_rgn_25 ) }
      else if( l_id == 26  ) { PP_INSTR_REG_END( wrk_rgn_26 ) } else if( l_id ==  27  ) { PP_INSTR_REG_END( wrk_rgn_27 ) }
      else if( l_id == 28  ) { PP_INSTR_REG_END( wrk_rgn_28 ) } else if( l_id ==  29  ) { PP_INSTR_REG_END( wrk_rgn_29 ) }
      else if( l_id == 30  ) { PP_INSTR_REG_END( wrk_rgn_30 ) } else if( l_id ==  31  ) { PP_INSTR_REG_END( wrk_rgn_31 ) }
      else if( l_id == 32  ) { PP_INSTR_REG_END( wrk_rgn_32 ) } else if( l_id ==  33  ) { PP_INSTR_REG_END( wrk_rgn_33 ) }
      else if( l_id == 34  ) { PP_INSTR_REG_END( wrk_rgn_34 ) } else if( l_id ==  35  ) { PP_INSTR_REG_END( wrk_rgn_35 ) }
      else if( l_id == 36  ) { PP_INSTR_REG_END( wrk_rgn_36 ) } else if( l_id ==  37  ) { PP_INSTR_REG_END( wrk_rgn_37 ) }
      else if( l_id == 38  ) { PP_INSTR_REG_END( wrk_rgn_38 ) } else if( l_id ==  39  ) { PP_INSTR_REG_END( wrk_rgn_39 ) }
      else { PP_INSTR_REG_END( wrk_rgn_undef ) }
#endif

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

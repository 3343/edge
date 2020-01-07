/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Scheduling for seismic simulations with point sources.
 **/

/*
 * Local control flow ids:
 *
 * max: first touch since synchronization
 * 0: ready to be progressed
 * 1: in progress
 * 2+: done in various steps
 */

// make sure we have our entries
static_assert( N_ENTRIES_CONTROL_FLOW == 8, "entries of control flow not matching" );

// update control flow if any of the work regions (local, src, neigh) finished
// outside of the regular loop for a consistent view in the next loop
for( unsigned short l_tg = 0; l_tg < m_timeGroups.size(); l_tg++ ) {
  unsigned short l_tgOff = N_ENTRIES_CONTROL_FLOW * l_tg;

  // initialize control flow if necessary
  if( m_cflow[l_tg][0] == std::numeric_limits< unsigned short >::max() ) {
    // reset all entries
    for( unsigned short l_en = 0; l_en < N_ENTRIES_CONTROL_FLOW; l_en++ ) {
      m_cflow[l_tg][l_en] = 0;
    }

    // mark initial sends as finished
    m_cflow[l_tg][6] = 2;

    // issue receives
    m_distributed.beginRecvs( m_timeGroups[l_tg]->nRecvSync()%2 == 0,
                              l_tg,
                              m_timeGroups[l_tg]->nRecvSync()%m_distributed.nCommBuffers() );
    m_cflow[l_tg][7] = 1;
  }

  // local, inner
  if( m_cflow[l_tg][0] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+0 ) ) {
    m_cflow[l_tg][0] = 2;
    m_timeGroups[l_tg]->updateTimePredInner();
  }
  // local, send
  if( m_cflow[l_tg][1] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+1 ) ) {
    m_cflow[l_tg][1] = 2;
    m_timeGroups[l_tg]->updateTimePredSend();
  }
  // src, inner
  if( m_cflow[l_tg][2] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+2 ) ) {
    m_cflow[l_tg][2] = 2;
  }
  // src, send
  if( m_cflow[l_tg][3] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+3 ) ) {
    m_cflow[l_tg][3] = 2;
  }
  // neigh, inner
  if( m_cflow[l_tg][4] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+4 ) ) {
    m_cflow[l_tg][4] = 2;
    m_timeGroups[l_tg]->updateDofUpInner();
  }
  // neigh, send
  if( m_cflow[l_tg][5] == 1 && m_shared.getStatusAll( parallel::Shared::FIN, l_tgOff+5 ) ) {
    m_cflow[l_tg][5] = 2;
    m_timeGroups[l_tg]->updateDofUpSend();
  }
  // MPI-send
  if( m_cflow[l_tg][6] == 1 && m_distributed.finSends( m_timeGroups[l_tg]->nSendSync()%2 == 1,
                                                       l_tg,
                                                       m_timeGroups[l_tg]->nSendSync()%m_distributed.nCommBuffers() ) ) {
    m_cflow[l_tg][6] = 2;
    m_timeGroups[l_tg]->updateSend();
  }
  // MPI-recv
  if( m_cflow[l_tg][7] == 1 && m_distributed.finRecvs( m_timeGroups[l_tg]->nRecvSync()%2 == 0,
                                                       l_tg,
                                                       m_timeGroups[l_tg]->nRecvSync()%m_distributed.nCommBuffers() ) ) {
    m_cflow[l_tg][7] = 2;
    m_timeGroups[l_tg]->updateRecv();
  }
}

// number of finished time groups
unsigned short l_nTgsFin = 0;

// resolve dependencies
for( unsigned short l_tg = 0; l_tg < m_timeGroups.size(); l_tg++ ) {
  unsigned short l_tgOff = N_ENTRIES_CONTROL_FLOW * l_tg;

  // local, send
  if( m_cflow[l_tg][1] == 0 && m_cflow[l_tg][6] == 2 ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+1 );

    m_cflow[l_tg][1] = 1;
    m_cflow[l_tg][6] = 0;
  }

  // local, inner
  if( m_cflow[l_tg][0] == 0 ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+0 );

    m_cflow[l_tg][0] = 1;
  }

  // schedule sends
  if( m_cflow[l_tg][1] == 2 && m_cflow[l_tg][6] == 0 ) {
    m_distributed.beginSends( m_timeGroups[l_tg]->nSendSync()%2 == 1,
                              l_tg,
                              m_timeGroups[l_tg]->nSendSync()%m_distributed.nCommBuffers() );

    m_cflow[l_tg][1] = 3;
    m_cflow[l_tg][6] = 1;
  }

  // schedule src, send
  if( m_cflow[l_tg][1] == 3 && m_cflow[l_tg][3] == 0 ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+3 );

    m_cflow[l_tg][1] = 4;
    m_cflow[l_tg][3] = 1;
  }

  // schedule src, inner
  if( m_cflow[l_tg][0] == 2 && m_cflow[l_tg][2] == 0 ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+2 );

    m_cflow[l_tg][0] = 4;
    m_cflow[l_tg][2] = 1;
  }

  // schedule neigh, send
  if( m_cflow[l_tg][7] == 2 &&
      m_cflow[l_tg][0] >= 4 && m_cflow[l_tg][1] >= 4 && m_cflow[l_tg][3] == 2 &&
      m_cflow[l_tg][5] == 0 && getTimePredAvailable(l_tg) ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+5 );

    m_cflow[l_tg][1] = 5;
    m_cflow[l_tg][3] = 3;
    m_cflow[l_tg][5] = 1;
    m_cflow[l_tg][7] = 0;
  }

  // schedule neigh, inner
  if( m_cflow[l_tg][0] >= 4 && m_cflow[l_tg][1] >= 4 && m_cflow[l_tg][2] == 2 &&
      m_cflow[l_tg][4] == 0 && getTimePredAvailable(l_tg) ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+4 );

    m_cflow[l_tg][0] = 5;
    m_cflow[l_tg][2] = 3;
    m_cflow[l_tg][4] = 1;

    // flush receivers if buffer size gets low
    m_recvs.flushIf();
  }

  // issue recvs
  if( m_cflow[l_tg][5] == 2 && m_cflow[l_tg][7] == 0 ) {
    // post next receive if not last time step
    if( !m_timeGroups[l_tg]->lastTimeStep() ) {
      m_distributed.beginRecvs( m_timeGroups[l_tg]->nRecvSync()%2 == 0,
                                l_tg,
                                m_timeGroups[l_tg]->nRecvSync()%m_distributed.nCommBuffers() );
      m_cflow[l_tg][7] = 1;
    }
    else {
      m_cflow[l_tg][7] = 3;
    }

    m_cflow[l_tg][5] = 3;
  }

  // time step for the time group complete
  if( m_cflow[l_tg][4] == 2 && m_cflow[l_tg][5] == 3 && getTimePredConsumed(l_tg) ) {
    m_timeGroups[l_tg]->updateTsInfo();

    if( !m_timeGroups[l_tg]->finished() ) {
      m_cflow[l_tg][0] = m_cflow[l_tg][1] = 0;
      m_cflow[l_tg][2] = m_cflow[l_tg][3] = 0;
      m_cflow[l_tg][4] = m_cflow[l_tg][5] = 0;
    }
    else {
      m_cflow[l_tg][4] = m_cflow[l_tg][5] = 4;
    }
  }

  // increase finished time group counter if appropiate
  if( m_timeGroups[l_tg]->finished() && m_cflow[l_tg][6] >= 2 && m_cflow[l_tg][7] >= 2 ) {
    l_nTgsFin++;
  }
}

// check if all time groups are finished
if( l_nTgsFin == m_timeGroups.size() ) {
  m_finished = true;
  return;
}
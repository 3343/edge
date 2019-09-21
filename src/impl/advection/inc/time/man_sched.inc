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
 * Scheduling for the advection equation.
 **/

/*
 * Local control flow ids:
 *
 * max: first touch since synchronization
 * 0: ready to be progressed
 * 1: in progress
 * 2: done
 * 3: done, processed by our own control flow
 */


// make sure we have our six entries
static_assert( N_ENTRIES_CONTROL_FLOW == 8, "entries of control flow not matching" );

#ifdef PP_USE_MPI
EDGE_LOG_FATAL << "fix me!";
#endif

// number of finished time groups
unsigned short l_nTgsFin = 0;

// iterate over the time groups
for( unsigned short l_tg = 0; l_tg < m_timeGroups.size(); l_tg++ ) {
  unsigned short l_tgOff = N_ENTRIES_CONTROL_FLOW * l_tg;

  // initialize control flow if necessary
  if( m_cflow[l_tg][0] == std::numeric_limits< unsigned short >::max() ) {
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+0 );
    m_shared.setStatusAll( parallel::Shared::RDY, l_tgOff+1 );
    m_cflow[l_tg][0] = m_cflow[l_tg][1] = 1;
    m_cflow[l_tg][6] = m_cflow[l_tg][7] = 0;
  }

  /*
   * update control flow with finished tasks
   */
  // local, inner
  if( m_cflow[l_tg][0] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, l_tgOff+0) ) {
    m_cflow[l_tg][0] = 2;
    m_timeGroups[l_tg]->updateTimePredInner();
  }

  // local, send
  if( m_cflow[l_tg][1] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, l_tgOff+1) ) {
    m_cflow[l_tg][1] = 2;
    m_timeGroups[l_tg]->updateTimePredSend();
  }

  // neigh, inner
  if( m_cflow[l_tg][6] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, l_tgOff+6) ) {
    m_cflow[l_tg][6] = 2;
    m_timeGroups[l_tg]->updateDofUpInner();
  }

  // neigh, send
  if( m_cflow[l_tg][7] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, l_tgOff+7) ) {
    m_cflow[l_tg][7] = 2;
    m_timeGroups[l_tg]->updateDofUpSend();
  }

  // schedule neighboring updates
  if( m_cflow[l_tg][0] == 2 && m_cflow[l_tg][1] == 2 && getTimePredAvailable(l_tg) ) {
    m_shared.setStatusAll(parallel::Shared::RDY, l_tgOff+6);
    m_shared.setStatusAll(parallel::Shared::RDY, l_tgOff+7);

    m_cflow[l_tg][0] = m_cflow[l_tg][1] = 3;
    m_cflow[l_tg][6] = m_cflow[l_tg][7] = 1;
  }

  // time step complete
  if( m_cflow[l_tg][6] == 2 && m_cflow[l_tg][7] == 2 && getTimePredConsumed(l_tg) ) {
    m_timeGroups[l_tg]->updateTsInfo();

    if( !m_timeGroups[l_tg]->finished() ) {
      m_shared.setStatusAll(parallel::Shared::RDY, l_tgOff+0);
      m_shared.setStatusAll(parallel::Shared::RDY, l_tgOff+1);
      m_cflow[l_tg][0] = m_cflow[l_tg][1] = 1;
    }
    else {
      m_cflow[l_tg][0] = m_cflow[l_tg][1] = 0;
    }

    m_cflow[l_tg][6] = m_cflow[l_tg][7] = 0;
  }

  // increase finished time group counter if appropiate
  if( m_timeGroups[l_tg]->finished() ) {
    l_nTgsFin++;
  }
}

// check if all time groups are finished
if( l_nTgsFin == m_timeGroups.size() ) {
  m_finished = true;
  return;
}
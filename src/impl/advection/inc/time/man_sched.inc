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
 * Scheduling for the advection equation.
 **/

/*
 * Local control flow ids:
 *
 * max: first touch since synchronization
 * 0: ready to be progressed
 * 1: in progress
 * 2: done
 */


// make sure we have our eight entries
static_assert( N_ENTRIES_CONTROL_FLOW == 7, "entries of control flow not matching" );

// initialize control flow if neccessary
if( m_cflow[0] == std::numeric_limits< unsigned short >::max() ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 0);
  m_shared.setStatusAll(parallel::Shared::RDY, 1);
  m_cflow[0] = 1;  m_cflow[1] = 1;

  m_mpi.beginRecvs(0, 0);
  m_cflow[5] = 1;

  m_cflow[2] = 0;
  m_cflow[3] = 0;
  m_cflow[4] = 0;
}

////////////////////////////////////////
// TODO: GLOBAL TIME STEPPING ONLY //
////////////////////////////////////////
/*
 * update control flow with finished tasks
 */

// inner-elements finished local computations
if( m_cflow[0] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 0) ) m_cflow[0] = 2;

// send-elements finished local computations
if( m_cflow[1] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 1) ) m_cflow[1] = 2;

// sends completed
if( m_cflow[2] == 1 && m_mpi.finSends(0, 0) ) m_cflow[2] = 2;

// inner-elements finished neighboring computations
if( m_cflow[3] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 3) ) m_cflow[3] = 2;

// send-elements finished neighboring computations
if( m_cflow[4] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 4) ) m_cflow[4] = 2;

// receive completed
if( m_cflow[5] == 1 && m_mpi.finRecvs(0, 0) ) m_cflow[5] = 2;

/*
 * Simple implementation: "Blocking"
 */
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[2] == 0 ) {
    m_mpi.beginSends(0, 0);
    m_cflow[2] = 1;
}

// check if we are finished
if( m_timeGroups[0]->finished() ) {
  // wait for communication to finish before signaling the other threads
  if( m_cflow[2] == 2 && m_cflow[5] == 2 ) {
    m_finished = true;
  }
  return;
}

// local work and receive done
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[5] == 2 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 3);
  m_shared.setStatusAll(parallel::Shared::RDY, 4);

  m_cflow[0] = m_cflow[1] = 0;
  m_cflow[3] = m_cflow[4] = 1;
  m_cflow[5] = 0;
}

// limiting complete
if( m_cflow[3] == 2 && m_cflow[4] == 2 ) {
  m_timeGroups[0]->updateTsInfo();

  if( !m_timeGroups[0]->finished() ) {
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_shared.setStatusAll(parallel::Shared::RDY, 1);
    m_cflow[0] = m_cflow[1] = 1;

    m_mpi.beginRecvs(0, 0);
    m_cflow[5] = 1;

    m_cflow[2] = 0;
  }
  else m_cflow[5] = 2;

  m_cflow[3] = m_cflow[4] = 0;
}

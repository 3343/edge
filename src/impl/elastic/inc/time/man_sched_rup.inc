/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Scheduling for elastics with rupture physics.
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
static_assert( N_ENTRIES_CONTROL_FLOW == 10, "entries of control flow not matching" );

// initialize control flow if neccessary
if( m_cflow[0] == std::numeric_limits< unsigned short >::max() ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 1);
  m_shared.setStatusAll(parallel::Shared::RDY, 0);
  m_cflow[0] = 1; m_cflow[1] = 1;

  m_mpi.beginRecvs(0);
  m_cflow[5] = 1;

  m_cflow[2] = 0;
  m_cflow[3] = 0;
  m_cflow[4] = 0;
  m_cflow[6] = 0;
  m_cflow[7] = 0;
  m_cflow[8] = 0;
  m_cflow[9] = 0;
}

/*
 * update control flow with finished tasks
 */

// inner-elements finished local computations
if( m_cflow[0] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 0) ) m_cflow[0] = 2;

// send-elements finished local computations
if( m_cflow[1] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 1) ) m_cflow[1] = 2;

// sends completed
if( m_cflow[2] == 1 && m_mpi.finSends(0) ) m_cflow[2] = 2;

// inner-faces finished rupture computations
if( m_cflow[6] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 6) ) m_cflow[6] = 2;

// send-/recv-faces finished rupture computations
if( m_cflow[7] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 7) ) m_cflow[7] = 2;

// receive completed
if( m_cflow[5] == 1 && m_mpi.finRecvs(0) ) m_cflow[5] = 2;

// inner-elements finished neighboring update
if( m_cflow[3] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 3) ) m_cflow[3] = 2;

// send-elements finished neighboring update
if( m_cflow[4] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 4) ) m_cflow[4] = 2;

// inner-elements finished limiting
if( m_cflow[8] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 8) ) m_cflow[8] = 2;

// send-elements finished limiting
if( m_cflow[9] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 9) ) m_cflow[9] = 2;

////////////////////////////////////////
// TODO: GLOBAL TIME STEPPING ONLY //
////////////////////////////////////////

// check if we are finished
if( m_timeGroups[0]->finished() ) {
  // wait for communication to finish before signaling the other threads
  EDGE_CHECK( m_cflow[2] == 2 ); // check that the sends are done
  EDGE_CHECK( m_cflow[5] == 2 ); // check that the recvs are done

  m_finished = true;
  return;
}

// finished local updates for send-elements, start sending the messages 
if( m_cflow[1] == 2 && m_cflow[2] == 0 ) {
    m_mpi.beginSends(0);
    m_cflow[2] = 1;
}

// local work and receive done, start rupture computations for send- and receive-faces
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[7] == 0 && m_cflow[5] == 2 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 7);
  m_cflow[7] = 1;

  // receive is done, new request might be posted once the current's info is processed
  m_cflow[5] = 0;
}

// local work done, start rupture computations for inner-faces
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[6] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 6);
  m_cflow[6] = 1;
}

// send-/recv-faces done with rupture computations, start neighboring updates for send-elements
if( m_cflow[7] == 2 && m_cflow[4] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 4);
  m_cflow[4] = 1;

  // both, 6 and 7 haven been scheduled at this point, prepare new execution of 0 and 1
  m_cflow[0] = m_cflow[1] = 0;
}

// inner-faces done with rupture computations, start neighboring updates for inner-elements
if( m_cflow[6] == 2 && m_cflow[3] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 3);
  m_cflow[3] = 1;

  // flush receivers if buffer size gets low
  m_recvs.flushIf();
  m_recvsSf.flushIf();
}

//  neighboring updates for send-elements done, get the MPI-recv out
if( m_cflow[4] == 2 && m_cflow[5] == 0 ) {
  if( !m_timeGroups[0]->lastTimeStep() ) {
    m_mpi.beginRecvs(0);
    m_cflow[5] = 1;
  }
  else {
    // default to a finished receive if no additional time steps are required after this one
    m_cflow[5] = 2;
  }
}

// neighboring updates for inner elements done, statrt limiting inner elements
if( m_cflow[3] == 2 && m_cflow[8] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 8);
  m_cflow[8] = 1;
}

// neighboring updates for send elements done, start limiting send elements
if( m_cflow[4] == 2 && m_cflow[9] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 9);
  m_cflow[9] = 1;
}

// neighboring updates for inner- and send-elements done, plus MPI-send is done, initiated local update of send elements
if( m_cflow[8] == 2 && m_cflow[9] == 2 && m_cflow[2] == 2 ) {
  EDGE_CHECK( m_cflow[1] == 0 );
  EDGE_CHECK( m_cflow[5] == 1 || m_cflow[5] == 2 );
  EDGE_CHECK( m_cflow[6] == 2 );
  EDGE_CHECK( m_cflow[7] == 2 );

  // this is the final and most restrictive condition of the time step, update the ts-info
  m_timeGroups[0]->updateTsInfo();

  if( !m_timeGroups[0]->finished() ) {
    // local update for send elements
    m_shared.setStatusAll(parallel::Shared::RDY, 1);
    m_cflow[1] = 1;

    // new send request can be posted once dependencies are resolved
    m_cflow[2] = 0;

    // reset the four cflow-entities
    m_cflow[6] = m_cflow[7] = 0;
    m_cflow[3] = m_cflow[4] = 0;
    m_cflow[8] = m_cflow[9] = 0;

    // also initiated the local updates of inner elements
    if( m_cflow[0] == 0 ) {
      m_shared.setStatusAll(parallel::Shared::RDY, 0);
      m_cflow[0] = 1;
    }
  }
}

// limiting for inner- and send-elements done, start the local part for the inner element update
// this condition is included in the statement above but less restrictive since local updates
// do not depend on sends
if( m_cflow[8] == 2 && m_cflow[9] == 2 && m_cflow[0] == 0 ) {
  // only schedule a new time step if don't we reach a synchronization point at the end of the time step
  if( !m_timeGroups[0]->finished() && !m_timeGroups[0]->lastTimeStep() ) {
    // not that condition above is also allowed to schedule 0
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_cflow[0] = 1;
  }
}

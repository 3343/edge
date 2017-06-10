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
 * Scheduling for elastics.
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
static_assert( N_ENTRIES_CONTROL_FLOW == 8, "entires of control flow not matching" );

// initialize control flow if neccessary
if( m_cflow[0] == std::numeric_limits< unsigned short >::max() ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 1);
  m_shared.setStatusAll(parallel::Shared::RDY, 0);
  m_cflow[0] = m_cflow[1] = 1;

  m_mpi.beginRecvs(0);
  m_cflow[5] = 1;

  m_cflow[2] = 0;
  m_cflow[3] = 0;
  m_cflow[4] = 0;
  m_cflow[6] = 0;
  m_cflow[7] = 0;
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
if( m_cflow[2] == 1 && m_mpi.finSends(0) ) m_cflow[2] = 2;

// inner-elements finished neighboring computations
if( m_cflow[3] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 3) ) m_cflow[3] = 2;

// send-elements finished neighboring computations
if( m_cflow[4] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 4) ) m_cflow[4] = 2;

// receive completed
if( m_cflow[5] == 1 && m_mpi.finRecvs(0) ) m_cflow[5] = 2;

// inner-elements finished source computations
if( m_cflow[6] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 6) ) m_cflow[6] = 2;

// send-elements finished source computations
if( m_cflow[7] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 7) ) m_cflow[7] = 2;

/*
 * Simple implementation: "Blocking"
 */
#ifdef PP_SCHED_BLOCKING
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[2] == 0 ) {
    m_mpi.beginSends(0);

    m_cflow[2] = 1;
}

// check if we are finished
if( m_timeGroups[0]->finished() ) {
  EDGE_CHECK( m_cflow[2] == 2 ); // check that the sends are done
  EDGE_CHECK( m_cflow[5] == 2 ); // check that the recvs are done

  m_finished = true;
  return;
}

// inner work and receive done
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[5] == 2 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 3);
  m_shared.setStatusAll(parallel::Shared::RDY, 4);

  // flush receivers if buffer size gets low
  m_recvs.flushIf();

  m_cflow[0] = m_cflow[1] = 0;
  m_cflow[3] = m_cflow[4] = 1;
  m_cflow[5] = 0;
}

// neighboring updates done
if( m_cflow[3] == 2 && m_cflow[4] == 2 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 6);
  m_shared.setStatusAll(parallel::Shared::RDY, 7);

  m_cflow[3] = m_cflow[4] = 0;
  m_cflow[6] = m_cflow[7] = 1;
}

// src updates and send done
if( m_cflow[6] == 2 && m_cflow[7] == 2 && m_cflow[2] == 2 ) {
  m_timeGroups[0]->updateTsInfo();

  if( !m_timeGroups[0]->finished() ) {
    m_mpi.beginRecvs(0);
    m_cflow[5] = 1;

    m_cflow[0] = m_cflow[1] = 1;
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_shared.setStatusAll(parallel::Shared::RDY, 1);
    m_cflow[2] = 0;
  }
  else m_cflow[5] = 2;

  m_cflow[6] = m_cflow[7] = 0;
}
#else
/*
 * Overlapping communication and computation
 */

// finished local updates for send-elements, start sending the messages 
if( m_cflow[1] == 2 && m_cflow[2] == 0 ) {
    m_mpi.beginSends(0);
    m_cflow[2] = 1;
}

// check if we are finished
if( m_timeGroups[0]->finished() ) {
  EDGE_CHECK( m_cflow[2] == 2 ); // check that the sends are done
  EDGE_CHECK( m_cflow[5] == 2 ); // check that the recvs are done

  m_finished = true;
  return;
}

// local work and receive done, start neighboring flux updates for send-elements
// the condition is identical to (***), however we post the work packages here,
// before the neighboring updates for the inner-elements, since the send-elements have priority
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[4] == 0 && m_cflow[5] == 2 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 4);
}

// local work done, start neighboring flux updates for inner-elements
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[3] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 3);
  m_cflow[3] = 1;
}

// local work and receive done: flush receivers if necessary and set dependent statuses (***)
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[4] == 0 && m_cflow[5] == 2 ) {
  // neigboring update for send-elements is in progress
  m_cflow[4] = 1;

  // flush receivers if buffer size gets low
  m_recvs.flushIf();

  // this call is more restrictive than the one before, therefore resets 0 and 1
  EDGE_CHECK( m_cflow[3] != 0 );
  m_cflow[0] = m_cflow[1] = 0;

  // receive is done, new request might be posted once the current's info is processed
  m_cflow[5] = 0;
}

// inner-elements done with neighboring flux updates, start src updates for inner-elements
if( m_cflow[3] == 2 && m_cflow[6] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 6);
  m_cflow[6] = 1;
}

// send-elements done with neighboring flux updates, start src updates for send-elements
if( m_cflow[4] == 2 && m_cflow[7] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 7);
  m_cflow[7] = 1;
}

// src updates for local- and inner-elements done plus MPI-send is done
// Remark: 0 depends on 4, considered as 7 here
if( m_cflow[6] == 2 && m_cflow[7] == 2 && m_cflow[2] == 2 ) {
  EDGE_CHECK( m_cflow[3] == 2 );
  EDGE_CHECK( m_cflow[4] == 2 );
  EDGE_CHECK( m_cflow[5] == 0 );

  // time step complete, update the ts-info
  m_timeGroups[0]->updateTsInfo();

  // only schedule a new time step if we didn't reach a synchronization point
  if( !m_timeGroups[0]->finished() ) {
    m_mpi.beginRecvs(0); // this could already be called after finishing 4
    m_cflow[5] = 1;

    m_shared.setStatusAll(parallel::Shared::RDY, 1);
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_cflow[0] = m_cflow[1] = 1;

    // new send request can be posted once dependencies are resolved
    m_cflow[2] = 0;
  }
  else m_cflow[5] = 2; // set receive to finished since we are not posting a new request

  // reset the four cflow-entities
  m_cflow[3] = m_cflow[4] = 0;
  m_cflow[6] = m_cflow[7] = 0;
}
#endif

/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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

// make sure we have our entries
static_assert( N_ENTRIES_CONTROL_FLOW == 16, "entries of control flow not matching" );

// get number of updates since sync
int_ts l_nUpsSync = m_timeGroups[0]->getUpdatesSync();

/*
 * derive MPI groups for current and upcoming (after time step)
 *   1) extrema
 *   2) candidate admissibility
 *   3) sub-cell tDOFs
 */
unsigned short l_tmpIds[4];
std::uintptr_t l_tmpIdData;

unsigned short l_ad[2], l_sc[2], l_ex[2];
edge::sc::Steering::getAdmIds( l_nUpsSync,
                               l_tmpIds );
l_tmpIdData = m_timeGroups[0]->getAdmDataId( l_tmpIds[1] );
l_ad[0] = m_mpi.getMg( l_tmpIdData );
l_tmpIdData = m_timeGroups[0]->getAdmDataId( l_tmpIds[3] );
l_ad[1] = m_mpi.getMg( l_tmpIdData );

edge::sc::Steering::getDofsIds( l_nUpsSync,
                                l_tmpIds );
l_tmpIdData = m_timeGroups[0]->getDofsScDataId( l_tmpIds[0] );
l_sc[0] = m_mpi.getMg( l_tmpIdData );
l_tmpIdData = m_timeGroups[0]->getDofsScDataId( l_tmpIds[1] );
l_sc[1] = m_mpi.getMg( l_tmpIdData );

edge::sc::Steering::getExtIds( l_nUpsSync,
                               l_tmpIds );
l_tmpIdData = m_timeGroups[0]->getExDataId( l_tmpIds[0] );
l_ex[0] = m_mpi.getMg( l_tmpIdData );
l_tmpIdData = m_timeGroups[0]->getExDataId( l_tmpIds[1] );
l_ex[1] = m_mpi.getMg( l_tmpIdData );

/*
 * initialize control flow if necessary
 */
if( m_cflow[0] == std::numeric_limits< unsigned short >::max() ) {
  // default for everyone: ready to be scheduled
  for( unsigned short l_cf = 0; l_cf < N_ENTRIES_CONTROL_FLOW; l_cf++ )
    m_cflow[l_cf] = 0;

  // compute local updates: send elements
  m_shared.setStatusAll(parallel::Shared::RDY, 1);
  m_cflow[1] = 1;

  // compute local updates: inner elements
  m_shared.setStatusAll(parallel::Shared::RDY, 0);
  m_cflow[0] = 1;

  // compute rupture net-updates: send/recv faces
  m_shared.setStatusAll(parallel::Shared::RDY, 4);
  m_cflow[4] = 1;

  // compute rupture net-updates: inner faces
  m_shared.setStatusAll(parallel::Shared::RDY, 3);
  m_cflow[3] = 1;

  // receive time integrated DG DOFs
  m_mpi.beginRecvs(0, 0);
  m_cflow[13] = 1;

  // receive candidate admissibility
  m_mpi.beginRecvs(0, l_ad[0] );
  m_cflow[15] = 1;

  // init communication
  m_cflow[ 2] = 0;
  m_cflow[ 7] = 2;
  m_cflow[ 8] = 2;
  m_cflow[11] = 2;
  m_cflow[12] = 2;
  m_cflow[14] = 2;
}

/*
 * update control flow with finished tasks
 */
// compute local, inner
if( m_cflow[ 0] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 0) ) m_cflow[0] = 2;

// compute local, send
if( m_cflow[ 1] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 1) ) m_cflow[1] = 2;

// MPI-send DG-DOFS
if( m_cflow[ 2] == 1 && m_mpi.finSends(0, 0) ) m_cflow[2] = 2;


// compute rupture, inner
if( m_cflow[ 3] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 3) ) m_cflow[3] = 2;

// compute rupture, send-/recv
if( m_cflow[ 4] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 4) ) m_cflow[4] = 2;


// compute neigh, inner-elements
if( m_cflow[ 5] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 5) ) m_cflow[5] = 2;

// compute neigh, send-elements
if( m_cflow[ 6] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 6) ) m_cflow[6] = 2;

// MPI-recv extrema
if( m_cflow[ 7] == 1 && m_mpi.finRecvs(0, l_ex[0]) ) m_cflow[7] = 2;

// MPI-send candidate admissibility
if( m_cflow[ 8] == 1 && m_mpi.finSends(0, l_ad[0]) ) m_cflow[8] = 2;


// compute limit+, inner
if( m_cflow[ 9] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 9) ) m_cflow[9] = 2;

// compute limit+, send
if( m_cflow[10] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 10) ) m_cflow[10] = 2;

// MPI-recv: SC-DOFs
if( m_cflow[11] == 1 && m_mpi.finRecvs(0, l_sc[0]) ) m_cflow[11] = 2;

// MPI-send: SC DOF
if( m_cflow[12] == 1 && m_mpi.finSends(0, l_sc[0]) ) m_cflow[12] = 2;

// MPI-recv: DG-DOF
if( m_cflow[13] == 1 &&m_mpi.finRecvs(0, 0) ) m_cflow[13] = 2;

// MPI-send: extrema
if( m_cflow[14] == 1 && m_mpi.finSends(0, l_ex[0]) ) m_cflow[14] = 2;

// MPI-receive: candidate admissibility
if( m_cflow[15] == 1 &&  m_mpi.finRecvs(0, l_ad[0]) ) m_cflow[15] = 2;

/*
 * resolve dependencies
 */
// MPI-send DG DOFs done: compute local send
if( m_cflow[2] == 2 && m_cflow[1] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 1);
  m_cflow[1] = 1;

  m_cflow[2] = 0;
}

// compute local, send done: MPI-send DG DOFs 
if( m_cflow[1] == 2 && m_cflow[2] == 0 ) {
    m_mpi.beginSends(0, 0);
    m_cflow[2] = 1;
}

// MPI-recv SC DOFs + MPI-send SC DOFs done: compute rupture for send/recv
if( m_cflow[11] == 2 && m_cflow[12] == 2 && m_cflow[4] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 4);
  m_cflow[4] = 1;
}

// compute inner local + compute send local + MPI-recv DG DOFs + MPI-recv extrema done: compute neigh, send
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[13] == 2 && m_cflow[7] == 2 && m_cflow[6] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 6);
  m_cflow[6] = 1;

  // allow for new MPI-recv extrema and MPI-send admissibility
  m_cflow[7] = 0;
  m_cflow[8] = 0;

  // flush DG receivers if buffer size gets low
  m_recvs.flushIf();
}

// compute local, inner + compute local, send done: compute neighboring for inner
if( m_cflow[0] == 2 && m_cflow[1] == 2 && m_cflow[5] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 5);
  m_cflow[5] = 1;
}

//  compute neigh, send done: MPI-recv extrema
if( m_cflow[6] == 2 && m_cflow[7] == 0 ) {
  m_mpi.beginRecvs(0, l_ex[1]);

  // this is fire and forget: we are only concerned with the extrema in the next time step (l_ex[0] then)
  m_cflow[7] = 2;
}

// compute rupture, send/recv + compute neigh, send done: send admissibility
if( m_cflow[4] == 2 && m_cflow[6] == 2 && m_cflow[8] == 0 ) {
  m_mpi.beginSends(0, l_ad[0]);
  m_cflow[8] = 1;

  // flush sub-cell receivers if buffer gets low
  m_recvsSf.flushIf();
}

//   compute rupture, inner + compute rupture, send/recv
// + compute neigh, inner + compute neigh, send
// + MPI-recv candidate admissibility done:
//     compute limit+, send
if(    m_cflow[3] == 2 && m_cflow[4] == 2
    && m_cflow[5] == 2 && m_cflow[6] == 2
    && m_cflow[15] == 2
    && m_cflow[10] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 10);
  m_cflow[10] = 1;

  // we obtained this time steps admissibility from adjacent ranks:
  // the limiter's sends have to be complete
  EDGE_CHECK_EQ( m_cflow[11], 2 );
  EDGE_CHECK_EQ( m_cflow[12], 2 );
  EDGE_CHECK_EQ( m_cflow[13], 2 );
  EDGE_CHECK_EQ( m_cflow[14], 2 );

  // reset communication
  m_cflow[11] = 0;
  m_cflow[12] = 0;
  m_cflow[13] = 0;
  m_cflow[14] = 0;
  m_cflow[15] = 0;
}

// compute rupture, inner + compute rupture, send/recv + compute neigh, inner + compute neigh, send done: compute limit+, inner
if( m_cflow[3] == 2 && m_cflow[4] == 2 && m_cflow[5] == 2 && m_cflow[6] == 2 && m_cflow[9] == 0 ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 9);
  m_cflow[9] = 1;
}

// compute limit+, send done:
//   MPI-recv SC-DOFs, MPI-recv SC-DOFs, MPI-recv DG-DOFs, MPI-send extrema, MPI-recv candidate admissibility
if( m_cflow[10] == 2 && m_cflow[11] == 0 ) {
  EDGE_CHECK_EQ( m_cflow[12], 0 );
  EDGE_CHECK_EQ( m_cflow[13], 0 );
  EDGE_CHECK_EQ( m_cflow[14], 0 );
  EDGE_CHECK_EQ( m_cflow[15], 0 );

  m_mpi.beginRecvs(0, l_sc[1]);
  m_mpi.beginSends(0, l_sc[1]);
  m_mpi.beginSends(0, l_ex[1]);

  // receive next time steps candidate admissibility and DG-DOFs only if there is another time step
  if( !m_timeGroups[0]->lastTimeStep() ) {
    m_mpi.beginRecvs(0, 0);
    m_mpi.beginRecvs(0, l_ad[1]);
  }

  // fire and forget comunication: relevant in next time step
  m_cflow[11] = 2;
  m_cflow[12] = 2;
  m_cflow[13] = 2;
  m_cflow[14] = 2;
  m_cflow[15] = 2;
}

// inner limit+ + send limit+ done
if( m_cflow[9] == 2 && m_cflow[10] == 2 ) {
  // perform safety checks
  EDGE_CHECK_EQ( m_cflow[ 0], 2 );
  EDGE_CHECK_EQ( m_cflow[ 1], 2 );

  EDGE_CHECK_EQ( m_cflow[ 3], 2 );
  EDGE_CHECK_EQ( m_cflow[ 4], 2 );
  EDGE_CHECK_EQ( m_cflow[ 5], 2 );
  EDGE_CHECK_EQ( m_cflow[ 6], 2 );
  EDGE_CHECK_EQ( m_cflow[ 8], 2 );
  // fire and forget
  EDGE_CHECK_EQ( m_cflow[ 7], 2 );
  EDGE_CHECK_EQ( m_cflow[11], 2 );
  EDGE_CHECK_EQ( m_cflow[12], 2 );
  EDGE_CHECK_EQ( m_cflow[13], 2 );
  EDGE_CHECK_EQ( m_cflow[14], 2 );
  EDGE_CHECK_EQ( m_cflow[15], 2 );

  // this is the final and most restrictive condition of the time step, update the ts-info
  m_timeGroups[0]->updateTsInfo();

  // schedule next time step if not finished
  if( !m_timeGroups[0]->finished() ) {
    // inner, local
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_cflow[0] = 1;

    // inner rupture
    m_shared.setStatusAll(parallel::Shared::RDY, 3);
    m_cflow[3] = 1;

    // reset compute tasks
    m_cflow[ 1] = 0;
    m_cflow[ 4] = 0;
    m_cflow[ 5] = 0;
    m_cflow[ 6] = 0;
    m_cflow[ 9] = 0;
    m_cflow[10] = 0;
  
    // set possible ongoing communication relevant for next time step
    m_cflow[ 7] = 1;
    m_cflow[11] = 1;
    m_cflow[12] = 1;
    m_cflow[13] = 1;
    m_cflow[14] = 1;
    m_cflow[15] = 1;

    // print a progress report
    if(   edge::parallel::g_rank==0
       && m_timeGroups[0]->getUpdatesPer()%25 == 0 ) {
      EDGE_LOG_INFO_ALL << "finished time step: #"
                        << m_timeGroups[0]->getUpdatesPer()
                        << ", time: "
                        << m_timeGroups[0]->getCovSimTime();
    }
  }
  else {
    // wait and progress possibly ongoing communication before leaving for synchronization
    while(    m_mpi.finSends( 0, 0       ) == false     // DG DOFs: send
           || m_mpi.finSends( 0, l_sc[1] ) == false     // SC DOFs: send
           || m_mpi.finRecvs( 0, l_sc[1] ) == false     // SC DOFS: recv
           || m_mpi.finSends( 0, l_ex[1] ) == false     // extrema: send
           || m_mpi.finRecvs( 0, l_ex[1] ) == false ) { // extrema: recv
      m_mpi.comm( m_shared.isSched(),
                  m_finished,
                  m_shared.isCommLead() );
    }

    m_finished = true;
    return;
  }
}
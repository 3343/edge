/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Scheduling for the shallow water equations.
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
static_assert( N_ENTRIES_CONTROL_FLOW == 2, "entries of control flow not matching" );

// initialize control flow if neccessary
if( m_cflow[0][0] == std::numeric_limits< unsigned short >::max() ) {
  m_shared.setStatusAll(parallel::Shared::RDY, 0);
  m_cflow[0][0] = 1;

  m_cflow[0][1] = 0;
}

// net-updates done
if( m_cflow[0][0] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 0) ) m_cflow[0][0] = 2;

// elements finished computations
if( m_cflow[0][1] == 1 && m_shared.getStatusAll(parallel::Shared::FIN, 1) ) m_cflow[0][1] = 2;

// check if we are finished
if( m_timeGroups[0]->finished() ) {
  m_finished = true;
  return;
}

if( m_cflow[0][0] == 2 && m_cflow[0][1] == 0 ) {
    m_shared.setStatusAll(parallel::Shared::RDY, 1);
    m_cflow[0][1] = 1;

    m_cflow[0][0] = 0;
}

// updates elements
if( m_cflow[0][1] == 2 && m_cflow[0][0] == 0 ) {
  m_timeGroups[0]->updateTsInfo();

  if( !m_timeGroups[0]->finished() ) {
    m_shared.setStatusAll(parallel::Shared::RDY, 0);
    m_cflow[0][0] = 1;
  }

  m_cflow[0][1] = 0;
}

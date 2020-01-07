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
 * MPI interface.
 **/
#include "MpiRemix.h"
#include "io/logging.h"
#include "global.h"

void edge::parallel::MpiRemix::init( unsigned short         i_nTgs,
                                     unsigned short         i_nElFas,
                                     std::size_t            i_nEls,
                                     std::size_t            i_nByFa,
                                     std::size_t    const * i_commStruct,
                                     unsigned short const * i_sendFa,
                                     std::size_t    const * i_sendEl,
                                     unsigned short const * i_recvFa,
                                     std::size_t    const * i_recvEl,
                                     data::Dynamic        & io_dynMem,
                                     std::size_t            i_nIterComm  ) {
  // call base init
  Distributed::init( i_nTgs,
                     i_nElFas,
                     i_nEls,
                     i_nByFa,
                     i_commStruct,
                     i_sendFa,
                     i_sendEl,
                     i_recvFa,
                     i_recvEl,
                     1,
                     io_dynMem );

  // init MPI-2 specifics
  m_nIterComm = i_nIterComm;

  std::size_t l_testSize = m_nChs * sizeof( int );
  m_sendTests = (int*) io_dynMem.allocate( l_testSize );
  m_recvTests = (int*) io_dynMem.allocate( l_testSize );

  std::size_t l_reqSize = m_nChs * sizeof( MPI_Request );
  m_sendReqs = (MPI_Request*) io_dynMem.allocate( l_reqSize );
  m_recvReqs = (MPI_Request*) io_dynMem.allocate( l_reqSize );

  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    m_sendTests[l_ch] = 0;
    m_sendReqs[l_ch]  = MPI_REQUEST_NULL;

    m_recvTests[l_ch] = 0;
    m_recvReqs[l_ch]  = MPI_REQUEST_NULL;
  }
}

void edge::parallel::MpiRemix::beginSends( bool           i_lt,
                                           unsigned short i_tg,
                                           unsigned short ) {
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    // only send if the message's time group matches
    bool l_match = checkSendTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // get the message on the way
    if( l_match ) {
      int l_err = MPI_Isend( m_sendBuffers + m_sendMsgs[l_ch].offL,
                             m_sendMsgs[l_ch].size,
                             MPI_BYTE,
                             m_sendMsgs[l_ch].rank,
                             m_sendMsgs[l_ch].tag,
                             MPI_COMM_WORLD,
                             &m_sendReqs[l_ch] );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      m_sendTests[l_ch] = 0;
    }
  }
}

void edge::parallel::MpiRemix::beginRecvs( bool           i_lt,
                                           unsigned short i_tg,
                                           unsigned short ) {
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // get the message on the way
    if( l_match ) {
      int l_err = MPI_Irecv( m_recvBuffers + m_recvMsgs[l_ch].offL,
                             m_recvMsgs[l_ch].size,
                             MPI_BYTE,
                             m_recvMsgs[l_ch].rank,
                             m_recvMsgs[l_ch].tag,
                             MPI_COMM_WORLD,
                             &m_recvReqs[l_ch] );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      m_recvTests[l_ch] = 0;
    }
  }
}

void edge::parallel::MpiRemix::comm() {
  for( std::size_t l_it = 0; l_it < m_nIterComm; l_it++ ) {
    std::size_t l_nFinSend = 0;
    std::size_t l_nFinRecv = 0;

    for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
      // test on send
      int l_test = -1;
      int l_err = MPI_Test( &m_sendReqs[l_ch],
                            &l_test,
                            MPI_STATUS_IGNORE );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      EDGE_CHECK_NE( l_test, -1 );
      if( l_test == 1 ) {
        l_nFinSend++;
        m_sendReqs[l_ch] = MPI_REQUEST_NULL;
        m_sendTests[l_ch] = l_test;
      }

      // test on receive
      l_test = -1;
      l_err = MPI_Test( &m_recvReqs[l_ch],
                        &l_test,
                        MPI_STATUS_IGNORE );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      EDGE_CHECK_NE( l_test, -1 );
      if( l_test == 1 ) {
        l_nFinRecv++;
        m_recvReqs[l_ch] = MPI_REQUEST_NULL;
        m_recvTests[l_ch] = l_test;
      }
    }

    // abort early if everything is finished already
    if( l_nFinSend == m_nChs && l_nFinRecv == m_nChs ) break;
  }
}

bool edge::parallel::MpiRemix::finSends( bool           i_lt,
                                         unsigned short i_tg,
                                         unsigned short ) const {
  // iterate over send messages of the time group
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkSendTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    if( l_match && m_sendTests[l_ch] != 1 ) return false;
  }

  return true;
}

bool edge::parallel::MpiRemix::finRecvs( bool           i_lt,
                                         unsigned short i_tg,
                                         unsigned short ) const {
  // iterate over recv messages of the time group
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    if( l_match && m_recvTests[l_ch] != 1 ) return false;
  }

  return true;
}
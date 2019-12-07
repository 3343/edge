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

void edge::parallel::MpiRemix::min( std::size_t      i_nVals,
                                    double         * i_vals,
                                    unsigned short * o_min ){
  // initialize to true
  for( std::size_t l_va = 0; l_va < i_nVals; l_va++ )
    o_min[l_va] = 1;

#ifdef PP_USE_MPI
  // global min variables
  std::vector< double > l_gVals;
  l_gVals.resize( i_nVals );

  // determine minimum values
  MPI_Allreduce( i_vals, &l_gVals[0], i_nVals, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

  std::vector< int > l_bidsIn( i_nVals );
  std::vector< int > l_bidsOut( i_nVals );

  // perform bidding
  for( std::size_t l_va = 0; l_va < i_nVals; l_va++ ) {
    // bid on the variable, if we match the minimum
    if( i_vals[l_va] == l_gVals[l_va] )
      l_bidsIn[l_va] = parallel::g_rank;
    else
      l_bidsIn[l_va] = std::numeric_limits< int >::max();
  }

  // determine minimum bidding rank
  MPI_Allreduce( &l_bidsIn[0], &l_bidsOut[0], i_nVals, MPI_INT, MPI_MIN, MPI_COMM_WORLD );

  // update results, if we won the bid
  for( std::size_t l_va = 0; l_va < i_nVals; l_va++ )
    if( l_bidsOut[l_va] != parallel::g_rank ) o_min[l_va] = 0;
#endif
}

bool edge::parallel::MpiRemix::checkSendTgLt( std::size_t    i_ch,
                                              bool           i_lt,
                                              unsigned short i_tg ) const {
#ifdef PP_USE_MPI
  // message's time group has to match
  bool l_match = (m_sendMsgs[i_ch].tg == i_tg);
  // either message is greater-equal or less-than was requested
  if( m_sendMsgs[i_ch].lt ) {
    l_match = l_match && i_lt;
  }
  return l_match;
#else
  return true;
#endif
}

bool edge::parallel::MpiRemix::checkRecvTgLt( std::size_t    i_ch,
                                              bool           i_lt,
                                              unsigned short i_tg ) const {
#ifdef PP_USE_MPI
  // message's time group has to match
  bool l_match = (m_recvMsgs[i_ch].tg == i_tg);
  // either message is greater-equal or less-than was requested
  if( m_recvMsgs[i_ch].lt ) {
    l_match = l_match && i_lt;
  }
  return l_match;
#else
  return true;
#endif
}

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
#ifdef PP_USE_MPI
  m_nIterComm = i_nIterComm;

  // derive the number of communication channels and communicating faces
  m_nChs = i_commStruct[0];
  m_nSeRe = (std::size_t *) io_dynMem.allocate( m_nChs * sizeof(std::size_t) );

  std::size_t l_sizeSend = 0;
  std::size_t l_sizeRecv = 0;
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    std::size_t l_tg    = i_commStruct[1 + l_ch*4 + 0];
    std::size_t l_tgAd  = i_commStruct[1 + l_ch*4 + 2];
    std::size_t l_nSeRe = i_commStruct[1 + l_ch*4 + 3];
    m_nSeRe[l_ch] = l_nSeRe;

    l_sizeSend += (l_tg > l_tgAd) ? 2 * l_nSeRe * i_nByFa : l_nSeRe * i_nByFa;
    l_sizeRecv += (l_tg < l_tgAd) ? 2 * l_nSeRe * i_nByFa : l_nSeRe * i_nByFa;
  }

  // allocate send and receive buffer
  l_sizeSend *= sizeof(unsigned char);
  l_sizeRecv *= sizeof(unsigned char);
  m_sendBuffer = (unsigned char*) io_dynMem.allocate( l_sizeSend );
  m_recvBuffer = (unsigned char*) io_dynMem.allocate( l_sizeRecv );

  // allocate pointer data structure
  m_sendPtrs = (unsigned char**) io_dynMem.allocate( i_nEls * i_nElFas * sizeof(unsigned char*) );
  m_recvPtrs = (unsigned char**) io_dynMem.allocate( i_nEls * i_nElFas * sizeof(unsigned char*) );

  // init with null pointers (no communication)
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < i_nElFas; l_fa++ ) {
      m_sendPtrs[ l_el*i_nElFas + l_fa ] = nullptr;
      m_recvPtrs[ l_el*i_nElFas + l_fa ] = nullptr;
    }
  }

  // init message structure and communication data
  m_sendMsgs = (t_msg*) io_dynMem.allocate( m_nChs * sizeof(t_msg) );
  m_recvMsgs = (t_msg*) io_dynMem.allocate( m_nChs * sizeof(t_msg) );

  std::size_t l_offSend = 0;
  std::size_t l_offRecv = 0;
  std::size_t l_first = 0;
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    // unpack channel: local time group, adjacent rank, adjacent time group, #comm faces
    std::size_t l_tg    = i_commStruct[1 + l_ch*4 + 0];
    std::size_t l_raAd  = i_commStruct[1 + l_ch*4 + 1];
    std::size_t l_tgAd  = i_commStruct[1 + l_ch*4 + 2];
    std::size_t l_nSeRe = i_commStruct[1 + l_ch*4 + 3];

    l_sizeSend = (l_tg > l_tgAd) ? l_nSeRe * i_nByFa * 2 : l_nSeRe * i_nByFa;
    l_sizeRecv = (l_tg < l_tgAd) ? l_nSeRe * i_nByFa * 2 : l_nSeRe * i_nByFa;

    // assign
    m_sendMsgs[l_ch].lt      = (l_tg < l_tgAd);
    m_sendMsgs[l_ch].tg      = l_tg;
    m_sendMsgs[l_ch].rank    = l_raAd;
    m_sendMsgs[l_ch].tag     = l_tg*i_nTgs + l_tgAd;
    m_sendMsgs[l_ch].size    = l_sizeSend;
    m_sendMsgs[l_ch].ptr     = m_sendBuffer + l_offSend;
    m_sendMsgs[l_ch].test    = 0;
    m_sendMsgs[l_ch].request = MPI_REQUEST_NULL;

    m_recvMsgs[l_ch].lt      = (l_tg < l_tgAd);
    m_recvMsgs[l_ch].tg      = l_tg;
    m_recvMsgs[l_ch].rank    = l_raAd;
    m_recvMsgs[l_ch].tag     = l_tgAd*i_nTgs + l_tg;
    m_recvMsgs[l_ch].size    = l_sizeRecv;
    m_recvMsgs[l_ch].ptr     = m_recvBuffer + l_offRecv;
    m_recvMsgs[l_ch].test    = 0;
    m_recvMsgs[l_ch].request = MPI_REQUEST_NULL;

    for( std::size_t l_co = 0; l_co < l_nSeRe; l_co++ ) {
      std::size_t l_seFa = i_sendFa[l_first+l_co];
      std::size_t l_seEl = i_sendEl[l_first+l_co];
      std::size_t l_reFa = i_recvFa[l_first+l_co];
      std::size_t l_reEl = i_recvEl[l_first+l_co];

      m_sendPtrs[ l_seEl*i_nElFas + l_seFa ] = m_sendBuffer+l_offSend;
      m_recvPtrs[ l_reEl*i_nElFas + l_reFa ] = m_recvBuffer+l_offRecv;

      l_offSend += (l_tg > l_tgAd) ? i_nByFa * 2 : i_nByFa;
      l_offRecv += (l_tg < l_tgAd) ? i_nByFa * 2 : i_nByFa;
    }
    l_first += l_nSeRe;
  }
#endif
}

edge::parallel::MpiRemix::MpiRemix( int    i_argc,
                                    char * i_argv[] ) {
      // set default values for non-mpi runs
      g_nRanks = 1;
      g_rank = 0;
      g_rankStr = std::to_string(0);
#ifdef PP_USE_MPI
      m_comm = MPI_COMM_WORLD;

      // initialize MPI, get size and rank
      if( g_nThreads == 1 ) {
        MPI_Init( &i_argc,
                  &i_argv );
      }
      else {
        int l_tdSu;
        MPI_Init_thread( &i_argc,
                         &i_argv,
                         MPI_THREAD_SERIALIZED,
                         &l_tdSu );
        // ensure the required threading support of MPI
        EDGE_CHECK( l_tdSu == MPI_THREAD_SERIALIZED );
      }
      MPI_Comm_size ( m_comm, &g_nRanks );
      MPI_Comm_rank( m_comm, &g_rank );
      MPI_Get_version( m_verStd, m_verStd+1 );
      g_rankStr = std::to_string( g_rank );
#endif
}

std::string edge::parallel::MpiRemix::getVerStr() {
  return std::to_string( m_verStd[0] ) + "." + std::to_string( m_verStd[1] );
}

void edge::parallel::MpiRemix::fin() {
#ifdef PP_USE_MPI
      MPI_Barrier(m_comm);
      MPI_Finalize();
#endif
}

void edge::parallel::MpiRemix::beginSends( bool           i_lt,
                                           unsigned short i_tg ) {
#ifdef PP_USE_MPI
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    // only send if the message's time group matches
    bool l_match = checkSendTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // get the message on the way
    if( l_match ) {
      int l_err = MPI_Isend( m_sendMsgs[l_ch].ptr,
                             m_sendMsgs[l_ch].size,
                             MPI_BYTE,
                             m_sendMsgs[l_ch].rank,
                             m_sendMsgs[l_ch].tag,
                             m_comm,
                             &m_sendMsgs[l_ch].request );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      m_sendMsgs[l_ch].test = 0;
    }
  }
#endif
}

void edge::parallel::MpiRemix::beginRecvs( bool           i_lt,
                                           unsigned short i_tg ) {
#ifdef PP_USE_MPI
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // get the message on the way
    if( l_match ) {
      int l_err = MPI_Irecv( m_recvMsgs[l_ch].ptr,
                             m_recvMsgs[l_ch].size,
                             MPI_BYTE,
                             m_recvMsgs[l_ch].rank,
                             m_recvMsgs[l_ch].tag,
                             m_comm,
                             &m_recvMsgs[l_ch].request );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      m_recvMsgs[l_ch].test = 0;
    }
  }
#endif
}

void edge::parallel::MpiRemix::comm() {
#ifdef PP_USE_MPI
  for( std::size_t l_it = 0; l_it < m_nIterComm; l_it++ ) {
    std::size_t l_nFinSend = 0;
    std::size_t l_nFinRecv = 0;

    for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
      // test on send
      int l_test = -1;
      int l_err = MPI_Test( &m_sendMsgs[l_ch].request,
                            &l_test,
                            MPI_STATUS_IGNORE );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      EDGE_CHECK_NE( l_test, -1 );
      if( l_test == 1 ) {
        l_nFinSend++;
        m_sendMsgs[l_ch].request = MPI_REQUEST_NULL;
        m_sendMsgs[l_ch].test = l_test;
      }

      // test on receive
      l_test = -1;
      l_err = MPI_Test( &m_recvMsgs[l_ch].request,
                        &l_test,
                        MPI_STATUS_IGNORE );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      EDGE_CHECK_NE( l_test, -1 );
      if( l_test == 1 ) {
        l_nFinRecv++;
        m_recvMsgs[l_ch].request = MPI_REQUEST_NULL;
        m_recvMsgs[l_ch].test = l_test;
      }
    }

    // abort early if everything is finished already
    if( l_nFinSend == m_nChs && l_nFinRecv == m_nChs ) break;
  }
#endif
}

bool edge::parallel::MpiRemix::finSends( bool           i_lt,
                                         unsigned short i_tg ) const {
#ifdef PP_USE_MPI
  // iterate over send messages of the time group
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkSendTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    if( l_match && m_sendMsgs[l_ch].test != 1 ) return false;
  }
#endif

  return true;
}

bool edge::parallel::MpiRemix::finRecvs( bool           i_lt,
                                         unsigned short i_tg ) const {
#ifdef PP_USE_MPI
  // iterate over recv messages of the time group
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    if( l_match && m_recvMsgs[l_ch].test != 1 ) return false;
  }
#endif

  return true;
}

void edge::parallel::MpiRemix::syncData( std::size_t           i_nByFa,
                                         unsigned char const * i_sendData,
                                         unsigned char       * o_recvData ) {
#ifdef PP_USE_MPI
  MPI_Request * l_sendReqs = new MPI_Request[ m_nChs ];
  MPI_Request * l_recvReqs = new MPI_Request[ m_nChs ];

  unsigned char const * l_sendPtr = i_sendData;
  unsigned char       * l_recvPtr = o_recvData;

  // issue communication
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    std::size_t l_size = m_nSeRe[l_ch] * i_nByFa;

    int l_err = MPI_Irecv( l_recvPtr,
                           l_size,
                           MPI_BYTE,
                           m_recvMsgs[l_ch].rank,
                           m_recvMsgs[l_ch].tag,
                           m_comm,
                           l_recvReqs+l_ch );
    EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

    l_err = MPI_Isend( l_sendPtr,
                       l_size,
                       MPI_BYTE,
                       m_sendMsgs[l_ch].rank,
                       m_sendMsgs[l_ch].tag,
                       m_comm,
                       l_sendReqs+l_ch );
    EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

    l_sendPtr += l_size;
    l_recvPtr += l_size;
  }

  // wait for communication to finish
  int l_err = MPI_Waitall( m_nChs,
                           l_recvReqs,
                           MPI_STATUSES_IGNORE );
  EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
  l_err = MPI_Waitall( m_nChs,
                       l_sendReqs,
                       MPI_STATUSES_IGNORE );
  EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

  delete[] l_sendReqs;
  delete[] l_recvReqs;
#endif
}
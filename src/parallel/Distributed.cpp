/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Basic interface for distributed memory parallelizations.
 **/
#include "Distributed.h"
#ifdef PP_USE_MPI
#include "mpi_wrapper.inc"
#endif

void edge::parallel::Distributed::sendCommBuffers2( unsigned short   i_tg,
                                                    unsigned short & o_cbL,
                                                    unsigned short & o_cbLtR,
                                                    unsigned short & o_cbGeR ) const {
  o_cbL = m_nSendsSync[i_tg]%2;
  o_cbLtR = std::numeric_limits< unsigned short >::max();
  if( m_nSendsSync[i_tg]%2 == 1 ) o_cbLtR = (m_nSendsSync[i_tg] / 2)%2;
  o_cbGeR = m_nSendsSync[i_tg]%2;
}

void edge::parallel::Distributed::recvCommBuffers2( unsigned short   i_tg,
                                                    unsigned short & o_cbLtL,
                                                    unsigned short & o_cbGeL ) const {
  o_cbLtL = (m_nRecvsSync[i_tg]/2)%2;
  o_cbGeL = m_nRecvsSync[i_tg]%2;
}

edge::parallel::Distributed::Distributed( int    i_argc,
                                          char * i_argv[] ) {
  // set default values for non-mpi runs
  g_nRanks = 1;
  g_rank = 0;
  g_rankStr = std::to_string(0);
#ifdef PP_USE_MPI
  // only check if MPI is already initialized
  int l_initialized = 0;
  int l_err = MPI_Initialized( &l_initialized );
  EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

  // init if not already done
  if( !l_initialized ) {
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
  }

  MPI_Comm_size ( MPI_COMM_WORLD, &g_nRanks );
  MPI_Comm_rank( MPI_COMM_WORLD, &g_rank );
  MPI_Get_version( m_verStd, m_verStd+1 );
  g_rankStr = std::to_string( g_rank );
#endif
}

void edge::parallel::Distributed::fin() {
#ifdef PP_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}

std::string edge::parallel::Distributed::getVerStr() {
  return std::to_string( m_verStd[0] ) + "." + std::to_string( m_verStd[1] );
}

bool edge::parallel::Distributed::checkSendTgLt( std::size_t    i_ch,
                                                 bool           i_lt,
                                                 unsigned short i_tg ) const {
  // message's time group has to match
  bool l_match = (m_sendMsgs[i_ch].tgL == i_tg);
  // either message is greater-equal or less-than was requested
  if( m_sendMsgs[i_ch].tgL < m_sendMsgs[i_ch].tgR ) {
    l_match = l_match && i_lt;
  }
  return l_match;
}

bool edge::parallel::Distributed::checkRecvTgLt( std::size_t    i_ch,
                                                 bool           i_lt,
                                                 unsigned short i_tg ) const {
  // message's time group has to match
  bool l_match = (m_recvMsgs[i_ch].tgL == i_tg);
  // either message is greater-equal or less-than was requested
  if( m_recvMsgs[i_ch].tgL < m_recvMsgs[i_ch].tgR ) {
    l_match = l_match && i_lt;
  }
  return l_match;
}

void edge::parallel::Distributed::init( unsigned short         i_nTgs,
                                        unsigned short         i_nElFas,
                                        std::size_t            i_nEls,
                                        std::size_t            i_nByFa,
                                        std::size_t    const * i_commStruct,
                                        unsigned short const * i_sendFa,
                                        std::size_t    const * i_sendEl,
                                        unsigned short const * i_recvFa,
                                        std::size_t    const * i_recvEl,
                                        unsigned short         i_nCommBuffers,
                                        data::Dynamic        & io_dynMem ) {
  m_nTgs = i_nTgs;
  m_nCommBuffers = i_nCommBuffers;
  EDGE_CHECK( m_nCommBuffers == 1 || m_nCommBuffers == 2 );

  m_nSendsSync = (std::size_t *) io_dynMem.allocate( m_nTgs * sizeof(std::size_t) );
  m_nRecvsSync = (std::size_t *) io_dynMem.allocate( m_nTgs * sizeof(std::size_t) );

  // derive the number of communication channels and communicating faces
  if( i_commStruct != nullptr ) {
    m_nChs = i_commStruct[0];
  }
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
  m_sendBufferSize = l_sizeSend;
  m_sendBuffers = (unsigned char*) io_dynMem.allocate( l_sizeSend*m_nCommBuffers );
  m_recvBufferSize = l_sizeRecv;
  m_recvBuffers = (unsigned char*) io_dynMem.allocate( l_sizeRecv*m_nCommBuffers );

  unsigned short l_nPtrsSend = m_nCommBuffers;
  unsigned short l_nPtrsRecv = (m_nCommBuffers == 2) ? 4 : 1;

  // allocate pointer data structure
  std::size_t l_sizePtrs = i_nEls * i_nElFas * sizeof(unsigned char*);
  unsigned char** l_sendPtrs = (unsigned char** ) io_dynMem.allocate( l_sizePtrs * l_nPtrsSend );
  unsigned char** l_recvPtrs = (unsigned char** ) io_dynMem.allocate( l_sizePtrs * l_nPtrsRecv );

  for( unsigned short l_po = 0; l_po < l_nPtrsSend; l_po++ )
    m_sendPtrs[l_po] = l_sendPtrs + l_po * i_nEls * i_nElFas;
  for( unsigned short l_po = 0; l_po < l_nPtrsRecv; l_po++ )
    m_recvPtrs[l_po] = l_recvPtrs + l_po * i_nEls * i_nElFas;

  // init with null pointers (no communication)
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < i_nElFas; l_fa++ ) {
      for( unsigned short l_po = 0; l_po < l_nPtrsSend; l_po++ )
        m_sendPtrs[l_po][l_el*i_nElFas + l_fa] = nullptr;
      for( unsigned short l_po = 0; l_po < l_nPtrsRecv; l_po++ )
        m_recvPtrs[l_po][l_el*i_nElFas + l_fa] = nullptr;
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
    m_sendMsgs[l_ch].tgL     = l_tg;
    m_sendMsgs[l_ch].tgR     = l_tgAd;
    m_sendMsgs[l_ch].rank    = l_raAd;
    m_sendMsgs[l_ch].tag     = l_tg*i_nTgs + l_tgAd;
    m_sendMsgs[l_ch].size    = l_sizeSend;
    m_sendMsgs[l_ch].offL    = l_offSend;

    m_recvMsgs[l_ch].tgL     = l_tg;
    m_recvMsgs[l_ch].tgR     = l_tgAd;
    m_recvMsgs[l_ch].rank    = l_raAd;
    m_recvMsgs[l_ch].tag     = l_tgAd*i_nTgs + l_tg;
    m_recvMsgs[l_ch].size    = l_sizeRecv;
    m_recvMsgs[l_ch].offL    = l_offRecv;

    for( std::size_t l_co = 0; l_co < l_nSeRe; l_co++ ) {
      std::size_t l_seFa = i_sendFa[l_first+l_co];
      std::size_t l_seEl = i_sendEl[l_first+l_co];
      std::size_t l_reFa = i_recvFa[l_first+l_co];
      std::size_t l_reEl = i_recvEl[l_first+l_co];

      for( unsigned short l_po = 0; l_po < l_nPtrsSend; l_po++ )
        m_sendPtrs[l_po][l_seEl*i_nElFas + l_seFa] = m_sendBuffers + (l_po%2)*m_sendBufferSize + l_offSend;

      for( unsigned short l_po = 0; l_po < l_nPtrsRecv; l_po++ ) {
        if( l_tg >= l_tgAd ) {
          m_recvPtrs[l_po][l_reEl*i_nElFas + l_reFa] = m_recvBuffers + (l_po%2)*m_recvBufferSize + l_offRecv;
        }
        else {
          m_recvPtrs[l_po][l_reEl*i_nElFas + l_reFa] = m_recvBuffers + ( (l_po/2)%2 )*m_recvBufferSize + l_offRecv;
        }
      }

      l_offSend += (l_tg > l_tgAd) ? i_nByFa * 2 : i_nByFa;
      l_offRecv += (l_tg < l_tgAd) ? i_nByFa * 2 : i_nByFa;
    }
    l_first += l_nSeRe;
  }

  reset();
}

void edge::parallel::Distributed::reset() {
  for( unsigned short l_tg = 0; l_tg < m_nTgs; l_tg++ ) {
    m_nRecvsSync[l_tg] = std::numeric_limits< std::size_t >::max();
    m_nSendsSync[l_tg] = std::numeric_limits< std::size_t >::max();
  }
}

void edge::parallel::Distributed::min( std::size_t      i_nVals,
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

void edge::parallel::Distributed::syncData( std::size_t           i_nByCh,
                                            std::size_t           i_nByFa,
                                            unsigned char const * i_sendData,
                                            unsigned char       * o_recvData ) {
#ifdef PP_USE_MPI
  MPI_Request * l_sendReqs = new MPI_Request[ m_nChs ];
  MPI_Request * l_recvReqs = new MPI_Request[ m_nChs ];

  unsigned char const * l_sendPtr = i_sendData;
  unsigned char       * l_recvPtr = o_recvData;

  // issue communication
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    std::size_t l_size = m_nSeRe[l_ch] * i_nByFa + i_nByCh;

    int l_err = MPI_Irecv( l_recvPtr,
                           l_size,
                           MPI_BYTE,
                           m_recvMsgs[l_ch].rank,
                           m_recvMsgs[l_ch].tag,
                           MPI_COMM_WORLD,
                           l_recvReqs+l_ch );
    EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

    l_err = MPI_Isend( l_sendPtr,
                       l_size,
                       MPI_BYTE,
                       m_sendMsgs[l_ch].rank,
                       m_sendMsgs[l_ch].tag,
                       MPI_COMM_WORLD,
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
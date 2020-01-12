/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * GASPI interface.
 **/
#include "Gaspi.h"
#include "io/logging.h"

void edge::parallel::Gaspi::ringWait( unsigned short     i_nReqs,
                                      gaspi_queue_id_t & io_queue ) {
  // get the maximum number of requests per queue
  gaspi_number_t l_nQueueReqsMax;
  gaspi_return_t l_err = gaspi_queue_size_max( &l_nQueueReqsMax );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

  // get the current number of requests for the given queue
  gaspi_number_t l_nQueueReqsCur;
  l_err = gaspi_queue_size(  io_queue,
                            &l_nQueueReqsCur );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

  // get the number of available queues
  gaspi_number_t l_nQueues;
  l_err = gaspi_queue_num( &l_nQueues );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

  if( l_nQueueReqsCur + i_nReqs > l_nQueueReqsMax ) {
    io_queue = (io_queue + 1) % l_nQueues;

    // flush the queue if its full
    l_err = gaspi_wait( io_queue,
                        GASPI_BLOCK );
    EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
  }                                  
}

edge::parallel::Gaspi::Gaspi( int    i_argc,
                              char * i_argv[] ): Distributed( i_argc, i_argv ) {
  // init
  gaspi_return_t l_err = gaspi_proc_init( GASPI_BLOCK );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

  // check for MPI-consistency
  gaspi_rank_t l_rank, l_nRanks;
  l_err = gaspi_proc_rank( &l_rank );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
  EDGE_CHECK_EQ( l_rank, g_rank );

  l_err = gaspi_proc_num( &l_nRanks );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
  EDGE_CHECK_EQ( l_nRanks, g_nRanks );

  // get meta-information
  l_err = gaspi_version( &m_verStd );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
}

void edge::parallel::Gaspi::init( unsigned short         i_nTgs,
                                  unsigned short         i_nElFas,
                                  std::size_t            i_nEls,
                                  std::size_t            i_nByFa,
                                  std::size_t    const * i_commStruct,
                                  unsigned short const * i_sendFa,
                                  std::size_t    const * i_sendEl,
                                  unsigned short const * i_recvFa,
                                  std::size_t    const * i_recvEl,
                                  data::Dynamic        & io_dynMem ) {
  // init parent
  Distributed::init( i_nTgs,
                     i_nElFas,
                     i_nEls,
                     i_nByFa,
                     i_commStruct,
                     i_sendFa,
                     i_sendEl,
                     i_recvFa,
                     i_recvEl,
                     2,
                     io_dynMem );

  // init receive counters
  std::size_t l_nRecvsSize = m_nTgs * sizeof(std::size_t);
  m_nRecvsOn = (std::size_t *) io_dynMem.allocate( l_nRecvsSize );
  for( unsigned short l_tg = 0; l_tg < m_nTgs; l_tg++ ) m_nRecvsOn[l_tg] = 0;

  // get remote offsets
  std::size_t l_offSize = m_nChs * sizeof(std::size_t);
  m_sendOffR = (std::size_t *) io_dynMem.allocate( l_offSize );
  std::size_t *l_recvOffL = new std::size_t[ m_nChs ];

  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    l_recvOffL[l_ch] = m_recvMsgs[l_ch].offL;
  }

  syncData( sizeof(std::size_t),
            0,
            (unsigned char *) l_recvOffL,
            (unsigned char *) m_sendOffR );

  delete[] l_recvOffL;

  // register segments
  for( unsigned short l_cb = 0; l_cb < 2; l_cb++ ) {
    gaspi_return_t l_err = gaspi_segment_use( m_sendSegs[l_cb],
                                              m_sendBuffers+(m_sendBufferSize*l_cb),
                                              m_sendBufferSize,
                                              GASPI_GROUP_ALL,
                                              GASPI_BLOCK,
                                              0 );
    EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

    l_err = gaspi_segment_use( m_recvSegs[l_cb],
                               m_recvBuffers+(m_recvBufferSize*l_cb),
                               m_recvBufferSize,
                               GASPI_GROUP_ALL,
                               GASPI_BLOCK,
                               0 );
    EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
  }

  // ensure availability of notification slots
  gaspi_number_t l_nNosMax;
  gaspi_return_t l_err = gaspi_notification_num( &l_nNosMax );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
  EDGE_CHECK_LE( m_nChs, l_nNosMax );

  // get remote notification ids
  std::size_t *l_notL = new std::size_t[ m_nChs ];
  std::size_t l_notSize = m_nChs * sizeof(std::size_t);
  m_notR = (std::size_t*) io_dynMem.allocate( l_notSize );
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ )
    l_notL[l_ch] = l_ch;

  syncData( sizeof(std::size_t),
            0,
            (unsigned char *) l_notL,
            (unsigned char *) m_notR );
  delete[] l_notL;
}

void edge::parallel::Gaspi::fin() {
  gaspi_return_t l_err = gaspi_proc_term( GASPI_BLOCK );
  EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );

  Distributed::fin();
}

void edge::parallel::Gaspi::beginSends( bool           i_lt,
                                        unsigned short i_tg ) {
  // increase or init send counter
  if( m_nSendsSync[i_tg] != std::numeric_limits< std::size_t >::max() ) {
    m_nSendsSync[i_tg]++;
  }
  else {
    m_nSendsSync[i_tg] = 0;
  }

  // derive local and remote comm buffers
  unsigned short l_cbL, l_cbLtR, l_cbGeR;
  sendCommBuffers2( i_tg,
                    l_cbL,
                    l_cbLtR,
                    l_cbGeR );

  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    // only send if the message's time group matches
    bool l_match = checkSendTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // get the message on the way
    if( l_match ) {
      // get a queue for the request
      ringWait( 1,
                m_queue );

      // derive remote communication buffer
      unsigned short l_cbR = l_cbGeR;
      if( m_sendMsgs[l_ch].tgL < m_sendMsgs[l_ch].tgR ) {
        l_cbR = l_cbLtR;
      }

      gaspi_return_t l_err = gaspi_write_notify( m_sendSegs[l_cbL],
                                                 m_sendMsgs[l_ch].offL,
                                                 m_sendMsgs[l_ch].rank,
                                                 m_recvSegs[l_cbR],
                                                 m_sendOffR[l_ch],
                                                 m_sendMsgs[l_ch].size,
                                                 m_notR[l_ch],
                                                 1,
                                                 m_queue,
                                                 GASPI_BLOCK );
      EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
    }
  }
}

void edge::parallel::Gaspi::beginRecvs( bool           i_lt,
                                        unsigned short i_tg ) {
  // increase or init recv counter
  if( m_nRecvsSync[i_tg] != std::numeric_limits< std::size_t >::max() ) {
    m_nRecvsSync[i_tg]++;
  }
  else {
    m_nRecvsSync[i_tg] = 0;
  }

  // set the number of ongoing receives
  EDGE_CHECK_EQ( m_nRecvsOn[i_tg], 0 );
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );
    if( l_match ) m_nRecvsOn[i_tg]++;
  }
}

bool edge::parallel::Gaspi::finSends( bool,
                                      unsigned short ) const {
  // current send buffer can be overwritten if the previous recvs are complete
  if( m_nSendsSync <= m_nRecvsSync+1 ) return true;

  return false;
}

bool edge::parallel::Gaspi::finRecvs( bool           i_lt,
                                      unsigned short i_tg ) const {
  // return immediately if nothing is ongoing
  if( m_nRecvsOn[i_tg] == 0 ) return true;

  // update ongoing counter otherwise
  unsigned short l_cbLtL, l_cbGeL;
  recvCommBuffers2( i_tg,
                    l_cbLtL,
                    l_cbGeL );

  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    bool l_match = checkRecvTgLt( l_ch,
                                  i_lt,
                                  i_tg );

    // check for notification
    if( l_match ) {
      unsigned short l_cbL = l_cbGeL;
      if( m_recvMsgs[l_ch].tgL < m_recvMsgs[l_ch].tgR )
        l_cbL = l_cbLtL;

      gaspi_notification_id_t l_not;
      gaspi_return_t l_err = gaspi_notify_waitsome( m_recvSegs[l_cbL],
                                                    l_ch,
                                                    1,
                                                    &l_not,
                                                    GASPI_TEST );
      EDGE_CHECK_NE( l_err, GASPI_ERROR );

      if( l_err == GASPI_SUCCESS ) {
        EDGE_CHECK_EQ( l_not, l_ch );
        gaspi_notification_t l_notVal;
        l_err = gaspi_notify_reset( m_recvSegs[l_cbL],
                                    l_ch,
                                    &l_notVal );
        EDGE_CHECK_EQ( l_err, GASPI_SUCCESS );
        EDGE_CHECK_EQ( l_notVal, 1 );

        EDGE_CHECK_NE( m_nRecvsOn[i_tg], 0 );
        m_nRecvsOn[i_tg]--;
        if( m_nRecvsOn[i_tg] == 0 ) return true;
      }
    }
  }

  return false;
}
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
 * MPI parallelization.
 **/

#include "Mpi.h"
#include "io/logging.h"
#include <set>

void edge::parallel::Mpi::start( int i_argc, char *i_argv[] ) {
      // set default values for non-mpi runs
      g_nRanks = 1;
      g_rank = 0;
      g_rankStr = std::to_string(0);
      m_verStd[0] = 0; m_verStd[1] = 0;
#ifdef PP_USE_MPI
      // initialize MPI, get size and rank
      if( g_nThreads == 1 ) {
        MPI_Init( &i_argc, &i_argv );
      }
      else {
        int l_tdSu;
        MPI_Init_thread( &i_argc, &i_argv, MPI_THREAD_FUNNELED, &l_tdSu );
        // ensure the required threading support of MPI
        EDGE_CHECK( l_tdSu == MPI_THREAD_FUNNELED );
      }
      m_mpiComm = MPI_COMM_WORLD;
      MPI_Comm_size (   m_mpiComm, &g_nRanks   );
      MPI_Comm_rank(    m_mpiComm, &g_rank     );
      MPI_Get_version(  m_verStd,   m_verStd+1 );
      g_rankStr = std::to_string(   g_rank );
#endif
}

void edge::parallel::Mpi::initLayout( const t_enLayout   &i_enLayout,
                                      const void         *i_data,
                                            std::size_t   i_bytesPerEntry,
                                            int_tg        i_tgFirst,
                                            int_tg        i_nTgGlo,
                                            unsigned int  i_iter ) {
#ifdef PP_USE_MPI
  EDGE_LOG_INFO << "  initializing entity specific MPI-layout";

  m_nMsgs = 0;
  m_nIterPerCheck = i_iter;
  m_cmmTdLast = -1;

  // prepare the messages
  m_send.resize( i_enLayout.timeGroups.size() );
  m_recv.resize( i_enLayout.timeGroups.size() );

  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    m_send[l_tg].resize( i_enLayout.timeGroups[l_tg].neRanks.size() );
    m_recv[l_tg].resize( i_enLayout.timeGroups[l_tg].neRanks.size() );

    m_nMsgs += m_send[l_tg].size() + m_recv[l_tg].size();
  }

  // assign the neighboring ranks and tags
  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK( i_enLayout.timeGroups[l_tg].neRanks.size() ==
                i_enLayout.timeGroups[l_tg].neTgs.size() );

    for( unsigned int l_ne = 0; l_ne < i_enLayout.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      m_send[l_tg][l_ne].rank = i_enLayout.timeGroups[l_tg].neRanks[l_ne];
      m_recv[l_tg][l_ne].rank = i_enLayout.timeGroups[l_tg].neRanks[l_ne];

      m_send[l_tg][l_ne].tag  =  (i_tgFirst+l_tg) * i_nTgGlo
                                + i_enLayout.timeGroups[l_tg].neTgs[l_ne];

      m_recv[l_tg][l_ne].tag  =   i_enLayout.timeGroups[l_tg].neTgs[l_ne] * i_nTgGlo
                                + i_tgFirst+l_tg;
    }
  }

  // init requests and test flags
  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    for( unsigned int l_ne = 0; l_ne < i_enLayout.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      m_send[l_tg][l_ne].test = 0;
      m_recv[l_tg][l_ne].test = 0;

      m_send[l_tg][l_ne].request = MPI_REQUEST_NULL;
      m_recv[l_tg][l_ne].request = MPI_REQUEST_NULL;
    }
  }

  // initialize pointers and size
  static_assert( sizeof(char) == 1, "char is assumed 1 byte in size" );
  char *l_data = (char*) i_data;

  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    l_data += i_enLayout.timeGroups[l_tg].inner.size * i_bytesPerEntry;
    for( unsigned int l_ne = 0; l_ne < m_send[l_tg].size(); l_ne++ ) {
      // check that we don't send empty messages
      EDGE_CHECK_NE( i_enLayout.timeGroups[l_tg].send[l_ne].size, 0 ) << l_tg << " " << l_ne;

      m_send[l_tg][l_ne].ptr  = l_data;
      m_send[l_tg][l_ne].size = i_enLayout.timeGroups[l_tg].send[l_ne].size * i_bytesPerEntry;
      l_data += m_send[l_tg][l_ne].size;
    }
    for( unsigned int l_ne = 0; l_ne < m_recv[l_tg].size(); l_ne++ ) {
      // check that we don't receive empty messages
      EDGE_CHECK_NE( i_enLayout.timeGroups[l_tg].receive[l_ne].size, 0 ) << l_tg << " " << l_ne;

      m_recv[l_tg][l_ne].ptr  = l_data;
      m_recv[l_tg][l_ne].size = i_enLayout.timeGroups[l_tg].receive[l_ne].size * i_bytesPerEntry;
      l_data += m_recv[l_tg][l_ne].size;
    }
  }

  // initiaize comm thread
  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    for( unsigned int l_ne = 0; l_ne < i_enLayout.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      m_send[l_tg][l_ne].cmmTd = -2;
      m_recv[l_tg][l_ne].cmmTd = -2;
    }
  }
#else
  // check that nothing is communicated for non-mpi settings
  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK( i_enLayout.timeGroups[l_tg].send.size()    == 0 );
    EDGE_CHECK( i_enLayout.timeGroups[l_tg].receive.size() == 0 );
  }
#endif
}

void edge::parallel::Mpi::comm(                bool  i_return,
                                const volatile bool &i_finished,
                                               bool  i_isLead ) {
#ifdef PP_USE_MPI
  std::set< volatile t_mpiMsg* > l_msgs;

  while( i_finished == false ) {
    // collect messages this thread is repsosible for
    for( int_tg l_tg = 0; l_tg < m_send.size(); l_tg++ ) {
      for( unsigned int l_ne = 0; l_ne < m_send[l_tg].size(); l_ne++ ) {
        // get send message
        volatile t_mpiMsg* l_send = &m_send[l_tg][l_ne];

        // add message if requrired
        if( (l_send->cmmTd == g_thread || l_send->cmmTd == -1) && l_send->test == 0 ) {
          l_msgs.insert( l_send );
        }

        // get recv message
        volatile t_mpiMsg* l_recv = &m_recv[l_tg][l_ne];

        // add message if required
        if( (l_recv->cmmTd == g_thread || l_recv->cmmTd == -1) && l_recv->test == 0 ) {
          l_msgs.insert( l_recv );
        }
      }
    }

    // progress communication
    for( unsigned int l_it = 0; l_it < m_nIterPerCheck; l_it++ ) {
      for ( std::set< volatile t_mpiMsg * >::iterator l_msg = l_msgs.begin(); l_msg != l_msgs.end(); ) {
        MPI_Request l_requ = (*l_msg)->request;
        int         l_test;

        MPI_Test( &l_requ,
                  &l_test,
                  MPI_STATUS_IGNORE );

        // in the case of more than one thread per message,
        // only the comm lead is allowed to signal results with the scheduling thread
        if( ( (*l_msg)->cmmTd == -1 && i_isLead ) ||
            ( (*l_msg)->cmmTd == g_thread       )    ) {
          // remove message if finished
          if( l_test == 1 ) {
            (*l_msg)->cmmTd = -2;
            (*l_msg)->test  =  1;
          }
        }

        // delete message form thread-local set in any case
        if( l_test == 1 ) {
          l_msg = l_msgs.erase(l_msg);
          continue;
        }

        l_msg++;
      }

      // exit if there's no communication
      if( l_msgs.size() == 0 ) break;
    }

    if( i_return == true ) break;
  }
#endif
}

void edge::parallel::Mpi::beginSends( int_tg i_tg ) {
#ifdef PP_USE_MPI
  for( std::size_t l_msg = 0; l_msg < m_send[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_mpiMsg *l_send = &m_send[i_tg][l_msg];

    MPI_Request l_req;
    int l_error = MPI_Isend(  l_send->ptr,
                              l_send->size,
                              MPI_BYTE,
                              l_send->rank,
                              l_send->tag,
                              m_mpiComm,
                             &l_req );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    // update message info
    l_send->request = l_req;
    l_send->test = 0;
    l_send->cmmTd = -1;
  }
#endif
}

void edge::parallel::Mpi::beginRecvs( int_tg i_tg ) {
#ifdef PP_USE_MPI
  for( std::size_t l_msg = 0; l_msg < m_recv[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_mpiMsg *l_recv = &m_recv[i_tg][l_msg];

    MPI_Request l_req;
    int l_error = MPI_Irecv(   l_recv->ptr,
                               l_recv->size,
                               MPI_BYTE,
                               l_recv->rank,
                               l_recv->tag,
                               m_mpiComm,
                              &l_req );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    // update message info
    l_recv->request = l_req;
    l_recv->test = 0;
    l_recv->cmmTd = -1;
  }
#endif
}

bool edge::parallel::Mpi::finSends( int_tg i_tg ) {
#ifdef PP_USE_MPI
  EDGE_CHECK( i_tg < m_send.size() );

  // iterate over send messages of the time group
  for( std::size_t l_msg = 0; l_msg < m_send[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_mpiMsg *l_send = &m_send[i_tg][l_msg];

    if( l_send->test == 0 ) return false;
  }
#endif

  return true;
}

bool edge::parallel::Mpi::finRecvs( int_tg i_tg ) {
#ifdef PP_USE_MPI
  EDGE_CHECK( i_tg < m_recv.size() );

  // iterate over send messages of the time group
  for( std::size_t l_msg = 0; l_msg < m_recv[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_mpiMsg *l_recv = &m_recv[i_tg][l_msg];

    if( l_recv->test == 0 ) return false;
  }
#endif

  return true;
}

#ifdef PP_USE_MPI
void edge::parallel::Mpi::iSendTgRg( char         const * i_buff,
                                     std::size_t          i_nBytesPerRg,
                                     unsigned int         i_nRgns,
                                     int          const * i_neRanks,
                                     MPI_Request        * o_requests,
                                     int                  i_tag ) {
  // check correct size of char type
  static_assert( sizeof(char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  char const * l_buffRgn = i_buff;

  // iterate over the send-regions
  for( unsigned int l_rg = 0; l_rg < i_nRgns; l_rg++ ) {
    int l_error = MPI_Isend( l_buffRgn,
                             i_nBytesPerRg,
                             MPI_BYTE,
                             i_neRanks[l_rg],
                             i_tag,
                             MPI_COMM_WORLD,
                             o_requests+l_rg );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    l_buffRgn += i_nBytesPerRg;
  }

}
#endif

#ifdef PP_USE_MPI
void edge::parallel::Mpi::iRecvTgRg( char               * o_buff,
                                     std::size_t          i_nBytesPerRg,
                                     unsigned int         i_nRgns,
                                     int          const * i_neRanks,
                                     MPI_Request        * o_requests,
                                     int                  i_tag ) {
  // check correct size of char type
  static_assert( sizeof(char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  char * l_buffRgn = o_buff;

  // iterate over the send-regions
  for( unsigned int l_rg = 0; l_rg < i_nRgns; l_rg++ ) {
    int l_error = MPI_Irecv( l_buffRgn,
                             i_nBytesPerRg,
                             MPI_BYTE,
                             i_neRanks[l_rg],
                             i_tag,
                             MPI_COMM_WORLD,
                             o_requests+l_rg );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    l_buffRgn += i_nBytesPerRg;
  }

}
#endif

#ifdef PP_USE_MPI
void edge::parallel::Mpi::iSendTgEn( char         const * i_buff,
                                     std::size_t          i_nBytesPerEn,
                                     unsigned int         i_nRgns,
                                     t_timeRegion const * i_sendRgns,
                                     int          const * i_neRanks,
                                     MPI_Request        * o_requests,
                                     int                  i_tag ) {
  // check correct size of char type
  static_assert( sizeof(char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  char const * l_buffRgn = i_buff;

  // iterate over the send-regions
  for( unsigned int l_rg = 0; l_rg < i_nRgns; l_rg++ ) {
    unsigned int l_size = i_nBytesPerEn * i_sendRgns[l_rg].size;

    // check that the message fits in int-type
    EDGE_CHECK_LT( l_size, std::numeric_limits< int >::max() );

    int l_error = MPI_Isend( l_buffRgn,
                             l_size,
                             MPI_BYTE,
                             i_neRanks[l_rg],
                             i_tag,
                             MPI_COMM_WORLD,
                             o_requests+l_rg );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    l_buffRgn += l_size;
  }
}
#endif

#ifdef PP_USE_MPI
void edge::parallel::Mpi::waitAll( unsigned int   i_nRequests,
                                   MPI_Request  * i_requests ) {

  int l_error = MPI_Waitall( i_nRequests,
                             i_requests,
                             MPI_STATUS_IGNORE );
  EDGE_CHECK_EQ( l_error, MPI_SUCCESS );
}
#endif

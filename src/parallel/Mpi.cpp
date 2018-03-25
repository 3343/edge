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
      m_comm = MPI_COMM_WORLD;
      MPI_Comm_size (   m_comm, &g_nRanks   );
      MPI_Comm_rank(    m_comm, &g_rank     );
      MPI_Get_version(  m_verStd,   m_verStd+1 );
      g_rankStr = std::to_string(   g_rank );
#endif
}

#ifdef PP_USE_MPI
void edge::parallel::Mpi::init( int_tg             i_tgFirst,
                                int_tg             i_nTgGlo,
                                std::uintptr_t     i_id,
                                unsigned short     i_grpId,
                                t_enLayout const & i_enLay,
                                t_grp            & o_grp ) {
  // set identifier
  o_grp.id = i_id;

  // prepare the messages
  o_grp.send.resize( i_enLay.timeGroups.size() );
  o_grp.recv.resize( i_enLay.timeGroups.size() );

  for( int_tg l_tg = 0; l_tg < i_enLay.timeGroups.size(); l_tg++ ) {
    o_grp.send[l_tg].resize( i_enLay.timeGroups[l_tg].neRanks.size() );
    o_grp.recv[l_tg].resize( i_enLay.timeGroups[l_tg].neRanks.size() );
  }

  // assign the neighboring ranks and tags
  for( int_tg l_tg = 0; l_tg < i_enLay.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK_EQ( i_enLay.timeGroups[l_tg].neRanks.size(),
                   i_enLay.timeGroups[l_tg].neTgs.size() );

    for( unsigned int l_ne = 0; l_ne < i_enLay.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      o_grp.send[l_tg][l_ne].rank = i_enLay.timeGroups[l_tg].neRanks[l_ne];
      o_grp.recv[l_tg][l_ne].rank = i_enLay.timeGroups[l_tg].neRanks[l_ne];
      o_grp.send[l_tg][l_ne].tag  =  i_grpId * i_nTgGlo * i_nTgGlo;
      o_grp.send[l_tg][l_ne].tag += (i_tgFirst+l_tg) * i_nTgGlo;
      o_grp.send[l_tg][l_ne].tag +=  i_enLay.timeGroups[l_tg].neTgs[l_ne];

      o_grp.recv[l_tg][l_ne].tag  =  i_grpId * i_nTgGlo * i_nTgGlo;
      o_grp.recv[l_tg][l_ne].tag +=  i_enLay.timeGroups[l_tg].neTgs[l_ne] * i_nTgGlo;
      o_grp.recv[l_tg][l_ne].tag +=  i_tgFirst+l_tg;
    }
  }

  // init requests and test flags
  for( int_tg l_tg = 0; l_tg < i_enLay.timeGroups.size(); l_tg++ ) {
    for( unsigned int l_ne = 0; l_ne < i_enLay.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      o_grp.send[l_tg][l_ne].test = 0;
      o_grp.recv[l_tg][l_ne].test = 0;

      o_grp.send[l_tg][l_ne].request = MPI_REQUEST_NULL;
      o_grp.recv[l_tg][l_ne].request = MPI_REQUEST_NULL;
    }
  }

  // initiaize comm thread
  for( int_tg l_tg = 0; l_tg < i_enLay.timeGroups.size(); l_tg++ ) {
    for( unsigned int l_ne = 0; l_ne < i_enLay.timeGroups[l_tg].neRanks.size(); l_ne++ ) {
      o_grp.send[l_tg][l_ne].cmmTd = -2;
      o_grp.recv[l_tg][l_ne].cmmTd = -2;
    }
  }
}
#endif

unsigned short edge::parallel::Mpi::addDefault( const t_enLayout     &i_enLayout,
                                                const void           *i_data,
                                                      std::size_t     i_bytesPerEntry,
                                                      int_tg          i_tgFirst,
                                                      int_tg          i_nTgGlo,
                                                      std::uintptr_t  i_id ) {
#ifdef PP_USE_MPI
  EDGE_LOG_INFO << "  adding MPI-group #" << m_grps.size() << " (default), id: " << i_id;

  // add a new group
  m_grps.resize( m_grps.size()+1 );

  init(  i_tgFirst,
         i_nTgGlo,
         i_id,
        (m_grps.size()-1),
         i_enLayout,
         m_grps.back() );

  // initialize pointers and size
  static_assert( sizeof(unsigned char) == 1, "unsigned char is assumed 1 byte in size" );
  unsigned char *l_data = (unsigned char*) i_data;

  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    l_data += i_enLayout.timeGroups[l_tg].inner.size * i_bytesPerEntry;
    for( unsigned int l_ne = 0; l_ne < m_grps.back().send[l_tg].size(); l_ne++ ) {
      m_grps.back().send[l_tg][l_ne].ptr  = l_data;
      m_grps.back().send[l_tg][l_ne].size = i_enLayout.timeGroups[l_tg].send[l_ne].size * i_bytesPerEntry;
      l_data += m_grps.back().send[l_tg][l_ne].size;
    }
    for( unsigned int l_ne = 0; l_ne < m_grps.back().recv[l_tg].size(); l_ne++ ) {
      m_grps.back().recv[l_tg][l_ne].ptr  = l_data;
      m_grps.back().recv[l_tg][l_ne].size = i_enLayout.timeGroups[l_tg].receive[l_ne].size * i_bytesPerEntry;
      l_data += m_grps.back().recv[l_tg][l_ne].size;
    }
  }

  return (m_grps.size() == 0 ) ? std::numeric_limits< unsigned short >::max() : m_grps.size()-1;
#else
  // check that nothing is communicated for non-mpi settings
  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK_EQ( i_enLayout.timeGroups[l_tg].send.size(),    0 );
    EDGE_CHECK_EQ( i_enLayout.timeGroups[l_tg].receive.size(), 0 );
  }

  return std::numeric_limits< unsigned short >::max();
#endif

}

unsigned short  edge::parallel::Mpi::addCustom( t_enLayout            const &i_enLayout,
                                                std::vector<
                                                    std::vector<
                                                      unsigned char *
                                                    >
                                                >                     const &i_sendData,
                                                std::vector<
                                                    std::vector<
                                                      unsigned char *
                                                    >
                                                >                     const &i_recvData,
                                                int_tg                       i_tgFirst,
                                                int_tg                       i_nTgGlo,
                                                std::uintptr_t               i_id ) {
#ifdef PP_USE_MPI
  EDGE_LOG_INFO << "  adding MPI-group #" << m_grps.size() << " (custom), id: " << i_id;

  // add a new group
  m_grps.resize( m_grps.size()+1 );

  init(  i_tgFirst,
         i_nTgGlo,
         i_id,
        (m_grps.size()-1),
         i_enLayout,
         m_grps.back() );

  // initialize pointers and size
  static_assert( sizeof(unsigned char) == 1, "char is assumed 1 byte in size" );

  // perform some initial checks on the data
  EDGE_CHECK_EQ( i_sendData.size(), i_recvData.size() );
  EDGE_CHECK_EQ( i_sendData.size(), i_enLayout.timeGroups.size() );

  // iterate over time groups in the custom regions
  for( std::size_t l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK_EQ( i_sendData[l_tg].size(),
                   i_recvData[l_tg].size() );
    EDGE_CHECK_EQ( i_enLayout.timeGroups[l_tg].send.size() + 1,
                   i_sendData[l_tg].size() );

    // iterate over MPI-regions
    for( std::size_t l_mr = 0; l_mr < i_enLayout.timeGroups[l_tg].send.size(); l_mr++ ) {
      m_grps.back().send[l_tg][l_mr].ptr  = i_sendData[l_tg][l_mr];
      m_grps.back().send[l_tg][l_mr].size = i_sendData[l_tg][l_mr+1] -  i_sendData[l_tg][l_mr];

      m_grps.back().recv[l_tg][l_mr].ptr  = i_recvData[l_tg][l_mr];
      m_grps.back().recv[l_tg][l_mr].size = i_recvData[l_tg][l_mr+1] -  i_recvData[l_tg][l_mr];
    }
  }

  return (m_grps.size() == 0 ) ? std::numeric_limits< unsigned short >::max() : m_grps.size()-1;
#else
  // check that nothing is communicated for non-mpi settings
  EDGE_CHECK_EQ( i_sendData.size(), 0 );
  EDGE_CHECK_EQ( i_recvData.size(), 0 );

  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    EDGE_CHECK_EQ( i_enLayout.timeGroups[l_tg].send.size(),    0 );
    EDGE_CHECK_EQ( i_enLayout.timeGroups[l_tg].receive.size(), 0 );
  }

  return std::numeric_limits< unsigned short >::max();
#endif
}

void edge::parallel::Mpi::comm(                bool  i_return,
                                const volatile bool &i_finished,
                                               bool  i_isLead ) {
#ifdef PP_USE_MPI
  std::set< volatile t_msg* > l_msgs;

  while( i_finished == false ) {
    // collect messages this thread is repsosible for
    for( unsigned short l_mg = 0; l_mg < m_grps.size(); l_mg++ ) {
      for( int_tg l_tg = 0; l_tg < m_grps[l_mg].send.size(); l_tg++ ) {
        for( unsigned int l_ne = 0; l_ne < m_grps[l_mg].send[l_tg].size(); l_ne++ ) {
          // get send message
          volatile t_msg* l_send = &m_grps[l_mg].send[l_tg][l_ne];

          // add message if requrired
          if( (l_send->cmmTd == g_thread || l_send->cmmTd == -1) && l_send->test == 0 ) {
            l_msgs.insert( l_send );
          }

          // get recv message
          volatile t_msg* l_recv = &m_grps[l_mg].recv[l_tg][l_ne];

          // add message if required
          if( (l_recv->cmmTd == g_thread || l_recv->cmmTd == -1) && l_recv->test == 0 ) {
            l_msgs.insert( l_recv );
          }
        }
      }
    }

    // progress communication
    for( unsigned int l_it = 0; l_it < m_nIterPerCheck; l_it++ ) {
      for ( std::set< volatile t_msg * >::iterator l_msg = l_msgs.begin(); l_msg != l_msgs.end(); ) {
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

unsigned short edge::parallel::Mpi::getMg( std::uintptr_t i_id ) {
  unsigned short l_mg = std::numeric_limits< unsigned short >::max();

#ifdef PP_USE_MPI
  // nullpointer ids are invalid, abort
  if( i_id == 0 ) return l_mg;

  for( std::size_t l_gr = 0; l_gr < m_grps.size(); l_gr++ )
    if( m_grps[l_gr].id == i_id ) {
      // check that there's only one group
      EDGE_CHECK_EQ( l_mg, std::numeric_limits< unsigned short >::max() );

      l_mg = l_gr;
    }
#endif

  return l_mg;
}

void edge::parallel::Mpi::beginSends( int_tg         i_tg,
                                      unsigned short i_mg ) {
  // ignore command, if group does not exist
  if( i_mg == std::numeric_limits< unsigned short >::max() ) return;

#ifdef PP_USE_MPI
  for( std::size_t l_msg = 0; l_msg < m_grps[i_mg].send[i_tg].size(); l_msg++ ) {
    // only issue non-empty messages
    if( m_grps[i_mg].send[i_tg][l_msg].size == 0 ) continue;

    // get message
    volatile t_msg *l_send = &m_grps[i_mg].send[i_tg][l_msg];

    MPI_Request l_req;
    int l_error = MPI_Isend(  l_send->ptr,
                              l_send->size,
                              MPI_BYTE,
                              l_send->rank,
                              l_send->tag,
                              m_comm,
                             &l_req );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    // update message info
    l_send->request = l_req;
    l_send->test = 0;
    l_send->cmmTd = -1;
  }
#endif
}

void edge::parallel::Mpi::beginRecvs( int_tg         i_tg,
                                      unsigned short i_mg ) {
  // ignore command, if group does not exist
  if( i_mg == std::numeric_limits< unsigned short >::max() ) return;

#ifdef PP_USE_MPI
  for( std::size_t l_msg = 0; l_msg < m_grps[i_mg].recv[i_tg].size(); l_msg++ ) {
    // only issue non-empty messages
    if( m_grps[i_mg].recv[i_tg][l_msg].size == 0 ) continue;

    // get message
    volatile t_msg *l_recv = &m_grps[i_mg].recv[i_tg][l_msg];

    MPI_Request l_req;
    int l_error = MPI_Irecv(   l_recv->ptr,
                               l_recv->size,
                               MPI_BYTE,
                               l_recv->rank,
                               l_recv->tag,
                               m_comm,
                              &l_req );
    EDGE_CHECK( l_error == MPI_SUCCESS );

    // update message info
    l_recv->request = l_req;
    l_recv->test = 0;
    l_recv->cmmTd = -1;
  }
#endif
}

bool edge::parallel::Mpi::finSends( int_tg         i_tg,
                                    unsigned short i_mg ) {
  // ignore command, if group does not exist
  if( i_mg == std::numeric_limits< unsigned short >::max() ) return true;

#ifdef PP_USE_MPI
  EDGE_CHECK_LT( i_mg, m_grps.size() );
  EDGE_CHECK_LT( i_tg, m_grps[i_mg].send.size() );

  // iterate over send messages of the time group
  for( std::size_t l_msg = 0; l_msg < m_grps[i_mg].send[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_msg *l_send = &m_grps[i_mg].send[i_tg][l_msg];

    if( l_send->test == 0 ) return false;
  }
#endif

  return true;
}

bool edge::parallel::Mpi::finRecvs( int_tg         i_tg,
                                    unsigned short i_mg ) {
  // ignore command, if group does not exist
  if( i_mg == std::numeric_limits< unsigned short >::max() ) return true;

#ifdef PP_USE_MPI
  EDGE_CHECK_LT( i_mg, m_grps.size() );
  EDGE_CHECK_LT( i_tg, m_grps[i_mg].recv.size() );

  // iterate over send messages of the time group
  for( std::size_t l_msg = 0; l_msg < m_grps[i_mg].recv[i_tg].size(); l_msg++ ) {
    // get message
    volatile t_msg *l_recv = &m_grps[i_mg].recv[i_tg][l_msg];

    if( l_recv->test == 0 ) return false;
  }
#endif

  return true;
}

#ifdef PP_USE_MPI
void edge::parallel::Mpi::iSendTgRg( unsigned char const * i_buff,
                                     std::size_t           i_nBytesPerRg,
                                     unsigned int          i_nRgns,
                                     int           const * i_neRanks,
                                     MPI_Request         * o_requests,
                                     int                   i_tag ) {
  // check correct size of char type
  static_assert( sizeof(unsigned char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  unsigned char const * l_buffRgn = i_buff;

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
void edge::parallel::Mpi::iRecvTgRg( unsigned char      * o_buff,
                                     std::size_t          i_nBytesPerRg,
                                     unsigned int         i_nRgns,
                                     int          const * i_neRanks,
                                     MPI_Request        * o_requests,
                                     int                  i_tag ) {
  // check correct size of char type
  static_assert( sizeof(unsigned char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  unsigned char * l_buffRgn = o_buff;

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
void edge::parallel::Mpi::iSendTgEn( unsigned char const * i_buff,
                                     std::size_t           i_nBytesPerEn,
                                     unsigned int          i_nRgns,
                                     t_timeRegion  const * i_sendRgns,
                                     int           const * i_neRanks,
                                     MPI_Request         * o_requests,
                                     int                   i_tag ) {
  // check correct size of char type
  static_assert( sizeof(unsigned char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  unsigned char const * l_buffRgn = i_buff;

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
void edge::parallel::Mpi::iRecvTgEn( unsigned char       * o_buff,
                                     std::size_t           i_nBytesPerEn,
                                     unsigned int          i_nRgns,
                                     t_timeRegion  const * i_recvRgns,
                                     int           const * i_neRanks,
                                     MPI_Request         * o_requests,
                                     int                   i_tag ) {
  // check correct size of char type
  static_assert( sizeof(unsigned char) == 1, "size of char not 1 byte" );

  // start of the current region's buffer
  unsigned char * l_buffRgn = o_buff;

  // iterate over the receive-regions
  for( unsigned int l_rg = 0; l_rg < i_nRgns; l_rg++ ) {
    unsigned int l_size = i_nBytesPerEn * i_recvRgns[l_rg].size;

    // check that the message fits in int-type
    EDGE_CHECK_LT( l_size, std::numeric_limits< int >::max() );

    int l_error = MPI_Irecv( l_buffRgn,
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
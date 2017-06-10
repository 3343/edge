/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
#ifndef MPI_H_
#define MPI_H_

#ifdef PP_USE_MPI
#include "mpi_wrapper.inc"
#endif
#include <string>
#include "data/layout.hpp"

#include "parallel/global.h"

namespace edge {
  namespace parallel {
    class Mpi;
  }
}

class edge::parallel::Mpi {
  private:
#ifdef PP_USE_MPI
  //! communicator
  MPI_Comm m_mpiComm;

  // MPI messages
  typedef struct {
    //! neighboring rank
    int         rank;
    //! mpi tag
    int         tag;
    //! test flag
    int         test;
    //! MPI request
    MPI_Request request;
    //! pointer to start of message
    void*       ptr;
    //! size of the message in bytes
    std::size_t size;
    //! responsible communication thread; -2 is inactive, -1 is all, 0+ is thread id
    int         cmmTd;
  } t_mpiMsg;

  //! outgoing messages [*][]: time region, [][*]: neigh rank
  std::vector< std::vector< t_mpiMsg > > m_send;

  //! incoming messages [*][]: time region, [][*]: neigh rank
  std::vector< std::vector< t_mpiMsg > > m_recv;

  //! last assigned cmm thread
  int m_cmmTdLast;

  //! number of messages present
  unsigned int m_nMsgs;

  //! number of iterations over comm list until a check for new work is performed
  unsigned int m_nIterPerCheck;
#endif

  public:
    //! max. version of the supported mpi-standard, 0: major, 1: minor
    static int m_verStd[2];

    /**
     * Initializes MPI.
     *
     * @param i_argc number of command line parameters.
     * @param i_argv values of command line parameters.
     **/
    void start( int i_argc, char *i_argv[] );

    /**
     * Initializes the communication layout.
     *
     * @param i_enLayout data layout of the entities which are communicated.
     * @param i_data pointer for first data point of the entity layout.
     * @param i_bytesPerEntry number of bytes per entry in the data layout.
     * @param i_tgFirst global time group associated with local time group 0.
     * @param i_nTgGl number of global time groups.
     * @param i_iter number iterations used, before comm. thread check for new assigned messages.
     **/
    void initLayout( const t_enLayout  &i_enLayout,
                     const void        *i_data,
                           std::size_t  i_bytesPerEntry,
                           int_tg       i_tgFirst,
                           int_tg       i_nTgGlo,
                           unsigned int i_iter=100 );

    /**
     * Progresses communication using the calling thread.
     *
     * @i_return if true, the function returns after the number of iterations specified in init;
     *           if false until finished is true.
     * @i_finished abort criterion if i_return is false.
     * @i_isLead true if the thread is the lead of the comm threads, signaling with the scheduling thread.
     *
     **/
    void comm(                bool  i_return,
               const volatile bool &i_finished,
                              bool  i_isLead=false );

    /**
     * Begins the send-operations for the specified time group.
     *
     * @param i_tg time group for which send-operations are issued.
     **/
    void beginSends( int_tg i_tg );

    /**
     * Begins the receive-operations for the specified time group.
     *
     * @param i_tg time group for which receive-operations are issued.
     **/
    void beginRecvs( int_tg i_tg );

    /**
     * Checks if all sends for the specified time group are finished.
     *
     * @param i_tg time group which is checked.
     * @return true if all sends are finished, false if sends are ongoing.
     **/
    bool finSends( int_tg i_tg );

    /**
     * Checks if all receives for the specified time group are finished.
     *
     * @param i_tg time group which is checked.
     * @return true if all receives are finished, false if receives are ongoing.
     **/
    bool finRecvs( int_tg i_tg );

    /**
     * Sends data for all send-regions of a time group.
     *
     * @param i_buff send buffer containing the data.
     * @param i_nBytesPerRg number of bytes every per-region datum occupies.
     * @param i_nRgns number of send-regions in the time group.
     * @param i_neRanks neighboring ranks.
     * @param o_requests will be set to the MPI_REQUESTs (one per send-region).
     * @param i_tag used mpi tag.
     **/
#ifdef PP_USE_MPI
    static void iSendTgRg( char         const * i_buff,
                           std::size_t          i_nBytesPerRg,
                           unsigned int         i_nRgns,
                           int          const * i_neRanks,
                           MPI_Request        * o_requests,
                           int                  i_tag=0 );
#endif

    /**
     * Receives data for all recv-regions of a time group.
     *
     * @param i_buff recv buffer containing the data.
     * @param i_nBytesPerRg number of bytes every per-region datum occupies.
     * @param i_nRgns number of receive-regions in the time group.
     * @param i_neRanks neighboring ranks.
     * @param o_requests will be set to the MPI_REQUESTs (one per send-region).
     * @param i_tag used mpi tag.
     **/
#ifdef PP_USE_MPI
    static void iRecvTgRg( char               * o_buff,
                           std::size_t          i_nBytesPerRg,
                           unsigned int         i_nRgns,
                           int          const * i_neRanks,
                           MPI_Request        * o_requests,
                           int                  i_tag=0 );
#endif

    /**
     * Sends data for all send-entities of a time group.
     *
     * @param i_buff send buffer containing the data.
     * @param i_nBytesPerEn number of bytes every per-entity datum occupies.
     * @param i_nRgns number of send-regions in the time group.
     * @param i_sendRgns entity layout of the send-regions.
     * @param i_neRanks neighboring ranks.
     * @param o_requests will be set to the MPI_REQUESTs (one per send-region).
     * @param i_tag used mpi tag.
     **/
#ifdef PP_USE_MPI
    static void iSendTgEn( char         const * i_buff,
                           std::size_t          i_nBytesPerEn,
                           unsigned int         i_nRgns,
                           t_timeRegion const * i_sendRgns,
                           int          const * i_neRanks,
                           MPI_Request        * o_requests,
                           int                  i_tag=0 );
#endif

#ifdef PP_USE_MPI
    /**
     * Waits until all MPI-requests are complete.
     *
     * @param i_nRequests number of MPI-requests.
     * @param i_requests MPI-requests.
     **/
    static void waitAll( unsigned int   i_nRequests,
                         MPI_Request  * i_requests );
#endif

    /**
     * Finalizes MPI.
     **/
    void fin() {
#ifdef PP_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
#endif
    }
};

#endif

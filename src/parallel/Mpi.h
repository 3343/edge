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
#ifndef EDGE_PARALLEL_MPI_H_
#define EDGE_PARALLEL_MPI_H_

#ifdef PP_USE_MPI
#include "mpi_wrapper.inc"
#endif
#include <string>
#include <cstdint>
#include <limits>
#include "data/EntityLayout.type"

#include "parallel/global.h"

namespace edge {
  namespace parallel {
    class Mpi;
  }
}

class edge::parallel::Mpi {
  private:
    //! number of iterations over comm list until a check for new work is performed
    unsigned int m_nIterPerCheck;

#ifdef PP_USE_MPI
    //! communicator
    MPI_Comm m_comm;

    //! MPI messages
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
      unsigned char* ptr;
      //! size of the message in bytes
      std::size_t size;
      //! responsible communication thread; -2 is inactive, -1 is all, 0+ is thread id
      int         cmmTd;
    } t_msg;

    //! MPI group, consisting of a collection of messages
    typedef struct {
      //! outgoing messages [*][]: time region, [][*]: neigh rank
      std::vector< std::vector< t_msg > > send;

      //! incoming messages [*][]: time region, [][*]: neigh rank
      std::vector< std::vector< t_msg > > recv;

      //! pointer-sized id of the group, defaults to 0 if not set
      std::uintptr_t id;
    } t_grp;

    //! mpi groups
    std::vector< t_grp > m_grps;

    /**
     * Initializes an MPI-group.
     *
     * @param i_tgFirst global time group associated with local time group 0.
     * @param i_nTgGlo number of global time groups.
     * @param i_id optional identifier of the group.
     * @param i_grpId local id of the group.
     * @param i_enLay entity layout for the group.
     * @param o_grp MPI group, which is initialized.
     **/
    void init( int_tg             i_tgFirst,
               int_tg             i_nTgGlo,
               std::uintptr_t     i_id,
               unsigned short     i_grpId,
               t_enLayout const & i_enLay,
               t_grp            & o_grp );
#endif

  public:
    //! max. version of the supported mpi-standard, 0: major, 1: minor
    static int m_verStd[2];

    /**
     * Constructor.
     *
     * @param i_iter number of iterations in message progressions.
     **/
    Mpi( unsigned int i_iter=100 ): m_nIterPerCheck( i_iter ){};

    /**
     * Initializes MPI.
     *
     * @param i_argc number of command line parameters.
     * @param i_argv values of command line parameters.
     **/
    void start( int i_argc, char *i_argv[] );

    /**
     * Adds the given default communication data.
     *
     * @param i_enLayout layout of the entities for which data is communicated.
     * @param i_data pointer for first data point of the entity layout.
     * @param i_bytesPerEntry number of bytes per entry in the data layout.
     * @param i_tgFirst global time group associated with local time group 0.
     * @param i_nTgGlo number of global time groups.
     * @param i_id identifier of this group (0 if not required).
     * 
     * @return internal id of the mpi group added MPI group (not the identifier). numeric_limits< unsigned short >::max() if none.
     **/
    unsigned short addDefault( const t_enLayout     &i_enLayout,
                               const void           *i_data,
                                     std::size_t     i_bytesPerEntry,
                                     int_tg          i_tgFirst,
                                     int_tg          i_nTgGlo,
                                     std::uintptr_t  i_id = 0 );

    /**
     * Adds the given custom communication data.
     *
     * @param i_enLayout layout of the entities for which data is communicated.
     * @param i_sendData pointers to the send region data of all time groups (one additional pointer per time group for size computations).
     * @param i_recvData pointers to the recv region data of all time groups (one additional pointer per time group for size computations).
     * @param i_tgFirst first considered time group.
     * @param i_tgFirst global time group associated with local time group 0.
     * @param i_nTgGlo number of global time groups.
     * @param i_id identifier of the group (0 if not required).
     *
     * @return internal id of mpi group added MPI group (not the identifier). numeric_limits< unsigned short >::max() if none.
     **/
    unsigned short addCustom( t_enLayout          const &i_enLayout,
                              std::vector<
                                std::vector<
                                  unsigned char *
                                >
                              >                   const &i_sendData,
                              std::vector<
                                std::vector<
                                  unsigned char *
                                >
                              >                   const &i_recvData,
                              int_tg                     i_tgFirst,
                              int_tg                     i_nTgGlo,
                              std::uintptr_t             i_id = 0 );

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
     * Determines the message group for the given identifier.
     *
     * @param i_id id of the MPI group.
     * @return message group, numeric_limits< unsigned short >::max() if not found.
     **/
    unsigned short getMg( std::uintptr_t i_id );

    /**
     * Begins the send-operations for the specified time group.
     * Empty messages are directly set to finished.
     *
     * @param i_tg time group for which send-operations are issued.
     * @param i_mg mpi group for which send-operations are issued.
     **/
    void beginSends( int_tg         i_tg,
                     unsigned short i_mg );

    /**
     * Begins the receive-operations for the specified time group.
     * Empty messages are directly set to finished.
     *
     * @param i_tg time group for which receive-operations are issued.
     * @param i_mg mpi group for which receive-operations are issued.
     **/
    void beginRecvs( int_tg         i_tg,
                     unsigned short i_mg );

    /**
     * Checks if all sends for the specified time group are finished.
     *
     * @param i_tg time group which is checked.
     * @param i_mg mpi group, which is checked.
     * @return true if all sends are finished, false if sends are ongoing.
     **/
    bool finSends( int_tg         i_tg,
                   unsigned short i_mg );

    /**
     * Checks if all receives for the specified time group are finished.
     *
     * @param i_tg time group which is checked.
     * @param i_mg message group, which is checked.
     * @return true if all receives are finished, false if receives are ongoing.
     **/
    bool finRecvs( int_tg         i_tg,
                   unsigned short i_mg );

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
    static void iSendTgRg( unsigned char const * i_buff,
                           std::size_t           i_nBytesPerRg,
                           unsigned int          i_nRgns,
                           int           const * i_neRanks,
                           MPI_Request         * o_requests,
                           int                   i_tag=0 );
#endif

    /**
     * Receives data for all recv-regions of a time group.
     *
     * @param o_buff recv buffer containing the data.
     * @param i_nBytesPerRg number of bytes every per-region datum occupies.
     * @param i_nRgns number of receive-regions in the time group.
     * @param i_neRanks neighboring ranks.
     * @param o_requests will be set to the MPI_REQUESTs (one per send-region).
     * @param i_tag used mpi tag.
     **/
#ifdef PP_USE_MPI
    static void iRecvTgRg( unsigned char       * o_buff,
                           std::size_t           i_nBytesPerRg,
                           unsigned int          i_nRgns,
                           int           const * i_neRanks,
                           MPI_Request         * o_requests,
                           int                   i_tag=0 );
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
    static void iSendTgEn( unsigned char const * i_buff,
                           std::size_t           i_nBytesPerEn,
                           unsigned int          i_nRgns,
                           t_timeRegion  const * i_sendRgns,
                           int           const * i_neRanks,
                           MPI_Request         * o_requests,
                           int                   i_tag=0 );
#endif

    /**
     * Receives data for all recv-entities of a time group.
     *
     * @param i_buff send buffer containing the data.
     * @param i_nBytesPerEn number of bytes every per-entity datum occupies.
     * @param i_nRgns number of send-regions in the time group.
     * @param i_recvRgns entity layout of the receive-regions.
     * @param i_neRanks neighboring ranks.
     * @param o_requests will be set to the MPI_REQUESTs (one per recv-region).
     * @param i_tag used mpi tag.
     **/
#ifdef PP_USE_MPI
    static void iRecvTgEn( unsigned char       * i_buff,
                           std::size_t           i_nBytesPerEn,
                           unsigned int          i_nRgns,
                           t_timeRegion  const * i_recvRgns,
                           int           const * i_neRanks,
                           MPI_Request         * o_requests,
                           int                   i_tag=0 );
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
     * Determines, if the calling rank holds the minimum for the values
     *
     * @param i_nVals number of values.
     * @param i_vals values.
     * @param o_min 1 if rank is holds minimum, 0 otherwise.
     */
    static void min( std::size_t      i_nVals,
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

    /*
     * Synchronizes the sparse types.
     *
     * @param i_enLayout entity layout.
     * @param io_chars entity characteristics.
     *
     * @paramt TL_T_EN_CHARS type of entity characteristics, offering member .spType.
     */
    template< typename TL_T_EN_CHARS >
    static void syncSpTypes( t_enLayout const &i_enLayout,
                             TL_T_EN_CHARS    *io_chars ) {
#ifdef PP_USE_MPI
      // iterate over time groups
      for( unsigned short l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
        std::size_t l_nRgns = i_enLayout.timeGroups[l_tg].neRanks.size();
        if( l_nRgns == 0 ) continue;

        int_el l_sendFirst = i_enLayout.timeGroups[l_tg].inner.first + i_enLayout.timeGroups[l_tg].inner.size;
        int_el l_recvFirst = i_enLayout.timeGroups[l_tg].inner.first + i_enLayout.timeGroups[l_tg].nEntsOwn;
        int_el l_nSend = i_enLayout.timeGroups[l_tg].nEntsOwn - i_enLayout.timeGroups[l_tg].inner.size;
        int_el l_nRecv = i_enLayout.timeGroups[l_tg].nEntsNotOwn;

        std::vector< int_spType  > l_buffSend( l_nSend );
        std::vector< int_spType  > l_buffRecv( l_nRecv );
        std::vector< MPI_Request > l_reqSend(  l_nRgns );
        std::vector< MPI_Request > l_reqRecv(  l_nRgns );

        // copy send data to buffer
        for( int_el l_se = 0; l_se < l_nSend; l_se++ )
          l_buffSend[l_se] = io_chars[ l_sendFirst + l_se ].spType;

        // exchange the data
        iSendTgEn( (unsigned char *) &l_buffSend[0],
                    sizeof(int_spType),
                    l_nRgns,
                    &i_enLayout.timeGroups[l_tg].send[0],
                    &i_enLayout.timeGroups[l_tg].neRanks[0],
                    &l_reqSend[0] );
        iRecvTgEn( (unsigned char *) &l_buffRecv[0],
                    sizeof(int_spType),
                    l_nRgns,
                    &i_enLayout.timeGroups[l_tg].receive[0],
                    &i_enLayout.timeGroups[l_tg].neRanks[0],
                    &l_reqRecv[0] );

        waitAll(   l_nRgns,
                 &l_reqSend[0] );
        waitAll(   l_nRgns,
                 &l_reqRecv[0] );

        // store the data
        for( int_el l_re = 0; l_re < l_nRecv; l_re++ )
          io_chars[ l_recvFirst + l_re ].spType = l_buffRecv[l_re];
      }
#endif
    }

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

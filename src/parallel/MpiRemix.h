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
#ifndef EDGE_PARALLEL_MPI_REMIX_H
#define EDGE_PARALLEL_MPI_REMIX_H

#ifdef PP_USE_MPI
#include <mpi.h>
#endif

#include <cstddef>
#include "data/Dynamic.h"

namespace edge {
  namespace parallel {
    class MpiRemix;
  }
}

/**
 * MPI interface.
 **/
class edge::parallel::MpiRemix {
  private:
#ifdef PP_USE_MPI
    //! communicator
    MPI_Comm m_comm;
#endif

    //! max. version of the supported mpi-standard, 0: major, 1: minor
    int m_verStd[2] = {0, 0};

    //! number of channels
    std::size_t m_nChs = 0;

    //! number of send and receiver element-face pairs
    std::size_t m_nSeRe = 0;

    //! send buffer
    unsigned char * m_sendBuffer = nullptr;

    //! receive buffer
    unsigned char * m_recvBuffer = nullptr;

    //! pointers into the send buffer
    unsigned char ** m_sendPtrs = nullptr;

    //! pointers into the receiver buffer
    unsigned char ** m_recvPtrs = nullptr;

    //! message structure
    typedef struct {
      //! true if message is part of a less-than LTS relation
      bool lt;

      //! local time group of this message
      unsigned short tg;

      //! neighboring rank
      int rank;

      //! mpi tag
      int tag;

      //! test flag
      volatile int test;

#ifdef PP_USE_MPI
      //! MPI request
      MPI_Request request;
#endif

      //! pointer to start of message
      unsigned char* ptr;

      //! size of the message in bytes
      std::size_t size;
    } t_msg;

    //! send messages
    t_msg * m_sendMsgs = nullptr;

    //! receive messages
    t_msg * m_recvMsgs = nullptr;

    //! number of iterations used in the message progression
    std::size_t m_nIterComm = 0;

    /**
     * Checks if the send message of the given channel matches the time group and less-than requirement.
     *
     * @param i_ch communication channel.
     * @param i_lt less-than relationship. If the message is less-than, than this argument has also be true for true return.
     * @param i_tg time group.
     * @return true if the conditions match, false if not.
     **/
    bool checkSendTgLt( std::size_t    i_ch,
                        bool           i_lt,
                        unsigned short i_tg ) const;

    /**
     * Checks if the recv message of the given channel matches the time group and less-than requirement.
     *
     * @param i_ch communication channel.
     * @param i_lt less-than relationship. If the message is less-than, than this argument has also be true for true return.
     * @param i_tg time group.
     * @return true if the conditions match, false if not.
     **/
    bool checkRecvTgLt( std::size_t    i_ch,
                        bool           i_lt,
                        unsigned short i_tg ) const;

  public:
    /**
     * Constructor.
     *
     * @param i_nTgs number of time groups.
     * @param i_nElFas number of faces per element.
     * @param i_nEls number of elements in the mesh.
     * @param i_nByFa number of bytes for every communicating face (excludes LTS).
     * @param i_commStruct communication structure as specified in EDGE-V.
     * @param i_sendFa send faces.
     * @param i_sendEl send elements.
     * @param i_recvFa receive faces.
     * @param i_recvEl receive elements.
     * @param io_dynMem dynamic memory allocations.
     * @param i_nIterComm number of iterations in the message-progresson calls.
     **/
    void init( unsigned short         i_nTgs,
               unsigned short         i_nElFas,
               std::size_t            i_nEls,
               std::size_t            i_nByFa,
               std::size_t    const * i_commStruct,
               unsigned short const * i_sendFa,
               std::size_t    const * i_sendEl,
               unsigned short const * i_recvFa,
               std::size_t    const * i_recvEl,
               data::Dynamic        & io_dynMem,
               std::size_t            i_nIterComm = 100 );

    /**
     * Destructor.
     **/
    ~MpiRemix(){};

    /**
     * Initializes MPI.
     *
     * @param i_argc number of command line parameters.
     * @param i_argv values of command line parameters.
     **/
    MpiRemix( int    i_argc,
              char * i_argv[] );

    /**
     * Finalizes MPI.
     **/
    void fin();

    /**
     * Calls MPI to initiate the sends for the given time group.
     *
     * @param i_lt if true sends are also issued for less-than LTS relations.
     * @param i_tg time group for which data is send.
     **/
    void beginSends( bool           i_lt,
                     unsigned short i_tg );

    /**
     * Calls MPI to initiate the receives for the given time group.
     *
     * @param i_lt if true receives are also issued for less-than LTS relations.
     * @param i_tg time group for which data is received.
     **/
    void beginRecvs( bool           i_lt,
                     unsigned short i_tg );

    /**
     * Progresses MPI communication.
     **/
    void comm();

    /**
     * Checks if all sends for the specified time group are finished.
     *
     * @param i_lt if true if sends for regions in less-than LTS relations should also be checked.
     * @param i_tg time group for which the sends are checked.
     * @return true if all sends are finished, false if sends are ongoing.
     **/
    bool finSends( bool           i_lt,
                   unsigned short i_tg ) const;

    /**
     * Checks if all receives for the specified time group are finished.
     *
     * @param i_lt if true if receives for regions in less-than LTS relations should also be checked.
     * @param i_tg time group for which the receives are checked.
     * @return true if all receives are finished, false if receives are ongoing.
     **/
    bool finRecvs( bool           i_lt,
                   unsigned short i_tg ) const;

    /**
     * Gets the pointers to the send buffer.
     *
     * @return send pointers.
     **/
    unsigned char ** getSendPtrs() { return m_sendPtrs; }

    /**
     * Gets the pointers to the receive buffer.
     *
     * @return receive pointers.
     **/
    unsigned char ** getRecvPtrs() { return m_recvPtrs; }
};

#endif
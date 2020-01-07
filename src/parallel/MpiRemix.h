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
#include "mpi_wrapper.inc"
#endif

#include "Distributed.h"

namespace edge {
  namespace parallel {
    class MpiRemix;
  }
}

/**
 * MPI interface.
 **/
class edge::parallel::MpiRemix: public Distributed {
  protected:
    //! send test flags
    volatile int *m_sendTests = nullptr;

    //! receive test flags
    volatile int *m_recvTests = nullptr;

    //! send requests
    MPI_Request *m_sendReqs = nullptr;

    //! receive requests
    MPI_Request *m_recvReqs = nullptr;

    //! number of iterations used in the message progression
    std::size_t m_nIterComm = 0;

  public:
    /**
     * Constructor.
     **/
    MpiRemix( int    i_argc,
              char * i_argv[] ): Distributed( i_argc,
                                              i_argv ){};

    /**
     * Initializes the MPI communication structure.
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
     * Calls MPI to initiate the sends for the given time group.
     *
     * @param i_lt if true sends are also issued for less-than LTS relations.
     * @param i_tg time group for which data is send.
     **/
    void beginSends( bool           i_lt,
                     unsigned short i_tg,
                     unsigned short );

    /**
     * Calls MPI to initiate the receives for the given time group.
     *
     * @param i_lt if true receives are also issued for less-than LTS relations.
     * @param i_tg time group for which data is received.
     **/
    void beginRecvs( bool           i_lt,
                     unsigned short i_tg,
                     unsigned short );

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
                   unsigned short i_tg,
                   unsigned short ) const;

    /**
     * Checks if all receives for the specified time group are finished.
     *
     * @param i_lt if true if receives for regions in less-than LTS relations should also be checked.
     * @param i_tg time group for which the receives are checked.
     * @return true if all receives are finished, false if receives are ongoing.
     **/
    bool finRecvs( bool           i_lt,
                   unsigned short i_tg,
                   unsigned short ) const;

    /**
     * Dummy reset, returns immediately.
     **/
    void reset(){};
};

#endif
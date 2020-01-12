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
#ifndef EDGE_PARALLEL_GASPI_H
#define EDGE_PARALLEL_GASPI_H

#include "Distributed.h"
#include <GASPI.h>

namespace edge {
  namespace parallel {
    class Gaspi;
  }
}

/**
 * GASPI interface.
 **/
class edge::parallel::Gaspi: public Distributed {
  private:
    //! version of the GASPI implementation
    float m_verStd = 0;

    //! number of ongoing receive "messages"
    std::size_t *m_nRecvsOn;

    //! segment id of the send buffers
    const gaspi_segment_id_t m_sendSegs[2] = {0, 1};

    //! segment id of the recv buffers
    const gaspi_segment_id_t m_recvSegs[2] = {2, 3};

    //! remote offsets of the sends
    std::size_t *m_sendOffR = nullptr;

    //! remote offsets of the recvs
    std::size_t *m_recvOffR = nullptr;

    //! remote notification ids
    std::size_t *m_notR = nullptr;

    //! current communication queue
    gaspi_queue_id_t m_queue = 0;

    /**
     * Gets a queue with a sufficient number of available communication requests.
     *
     * @param i_nReqs number communication requests which should be available.
     * @param io_queue current queue which is updated with the next one if it does not have a sufficient number of requests.
     **/
    static void ringWait( unsigned short     i_nReqs,
                          gaspi_queue_id_t & io_queue );

  public:
    /**
     * Initializes GASPI.
     *
     * @param i_argc number of command line parameters.
     * @param i_argv values of command line parameters.
     **/
    Gaspi( int    i_argc,
           char * i_argv[] );

    /**
     * Initializes the GASPI communication structure.
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
     * @param io_dynMem will be used for dynamic memory allocations.
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
               data::Dynamic        & io_dynMem );

    /**
     * Finalizes GASPI (and MPI) if initialized.
     **/
    void fin();

    /**
     * Uses RMA-writes to initiate the sends for the given time group.
     *
     * @param i_lt if true sends are also issued for less-than LTS relations.
     * @param i_tg time group for which data is send.
     **/
    void beginSends( bool           i_lt,
                     unsigned short i_tg );

    /**
     * Begins the "receives" by initializing the number of expected messages.
     *
     * @param i_lt if true receives are also issued for less-than LTS relations.
     * @param i_tg time group for which data is received.
     **/
    void beginRecvs( bool           i_lt,
                     unsigned short i_tg );

    /**
     * Checks if the current send buffer for the specified time group can be overwritten.
     * This means that the previous (double buffer) receives were successful. 
     *
     * @return true if all sends are finished, false if sends are ongoing.
     **/
    bool finSends( bool,
                   unsigned short ) const;

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
     * Dummy function, returns immediately.
     **/
    void comm(){};

    /**
     * Resets the internal status of the dual send/recv buffer.
     **/
    void reset();
};

#endif
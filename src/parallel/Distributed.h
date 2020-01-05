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
#ifndef EDGE_PARALLEL_DISTRIBUTED_H
#define EDGE_PARALLEL_DISTRIBUTED_H

#include <cstddef>
#include "data/Dynamic.h"

namespace edge {
  namespace parallel {
    class Distributed;
  }
}

/**
 * Distributed memory interface.
 **/
class edge::parallel::Distributed {
  protected:
    //! number of time groups
    std::size_t m_nTgs = 0;

    //! number of channels
    std::size_t m_nChs = 0;

    //! number of send-receive faces
    std::size_t *m_nSeRe = nullptr;

    //! size of the send buffer
    std::size_t m_sendBufferSize = 0;

    //! send buffer
    unsigned char * m_sendBuffer = nullptr;

    //! size of the receive buffer
    std::size_t m_recvBufferSize = 0;

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

      //! size of the message in bytes
      std::size_t size;

      //! offset of the message in local memory
      std::size_t offL;
    } t_msg;

    //! send messages
    t_msg * m_sendMsgs = nullptr;

    //! receive messages
    t_msg * m_recvMsgs = nullptr;

    /**
     * Initializes the communication structure.
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
};

#endif
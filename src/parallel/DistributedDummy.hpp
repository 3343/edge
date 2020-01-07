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
 * Dummy interface to the distributed memory parallelization.
 **/
#ifndef EDGE_PARALLEL_DISTRIBUTED_DUMMY_H
#define EDGE_PARALLEL_DISTRIBUTED_DUMMY_H

#include "Distributed.h"

namespace edge {
  namespace parallel {
    class DistributedDummy;
  }
}

/**
 * Dummy implementation for shared memory runs.
 **/
class edge::parallel::DistributedDummy: public Distributed {

  public:
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
               data::Dynamic        & io_dynMem ) {
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
    }

    /**
     * Constructor.
     **/
    DistributedDummy( int    i_argc,
                      char * i_argv[] ): Distributed( i_argc,
                                                      i_argv ){}; 

    /**
     * Dummy sends, returns immediately.
     **/
    void beginSends( bool,
                     unsigned short,
                     unsigned short ){};

    /**
     * Dummy receives, returns immediately.
     **/
    void beginRecvs( bool,
                     unsigned short,
                     unsigned short ){};

    /**
     * Dummy progress, returns immediately.
     **/
    void comm(){};

    /**
     * Dummy check, always returns true.
     **/
    bool finSends( bool,
                   unsigned short,
                   unsigned short ) const { return true; };

    /**
     * Dummy check, always returns true.
     **/
    bool finRecvs( bool,
                   unsigned short,
                   unsigned short ) const { return true; };

    /**
     * Dummy reset, returns immediately.
     **/
    void reset(){};
};

#endif
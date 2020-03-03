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
 * Writes background meshes in msh4 format.
 **/
#ifndef EDGE_V_IO_BGMESHMSH4_H
#define EDGE_V_IO_BGMESHMSH4_H

#include <iostream>
#include "../constants.h"

namespace edge_v {
  namespace io {
    class BgMeshMsh4;
  }
}

class edge_v::io::BgMeshMsh4 {
  public:
    /**
     * Writes the string representation of the msh4 background mesh.
     *
     * @param i_elTy element type.
     * @param i_nVes number of vertices.
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds vertex coordinates.
     * @param i_targetLengths target lengths at the vertices.
     * @param io_stream stream which will be written.
     **/
    static void write( t_entityType                i_elTy,
                       t_idx                       i_nVes,
                       t_idx                       i_nEls,
                       t_idx              const  * i_elVe,
                       double             const (* i_veCrds)[3],
                       float              const  * i_targetLengths,
                       std::ostream              & io_stream );
};

#endif
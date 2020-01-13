/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2018, Regents of the University of California
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
 * Writes a GMSH view, based on the given coordinates and values.
 **/
#ifndef EDGE_V_IO_GMSH_VIEW_HPP
#define EDGE_V_IO_GMSH_VIEW_HPP

#include <string>
#include "constants.h"

namespace edge_v {
  namespace io {
    class GmshView;
  }
}

/**
 * Inteface to Gmsh's view.
 **/
class edge_v::io::GmshView {
  public:
    /**
     * Writes the Gmsh-view.
     * 
     * @param i_path path of the output file.
     * @param i_elType element type: 'tria3' or 'tet4'
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds coordinates of the vertices.
     * @param i_values values, which are assigned to the nodes.
     **/
    static void write( std::string const  & i_path,
                       t_entityType         i_elType,
                       std::size_t          i_nEls,
                       std::size_t const (* i_elVe),
                       double      const (* i_veCrds)[3],
                       float       const  * i_values );
};

#endif
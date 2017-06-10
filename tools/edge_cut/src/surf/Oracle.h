/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Oracle for the surface mesh.
 **/

#ifndef ORACLE_H_
#define ORACLE_H_

#include <string>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include "constants.hpp"
#include "Topo.h"

namespace edge_cut {
  namespace surf {
    class Oracle;
  }
}

class edge_cut::surf::Oracle {
  private:
    //! topography
    Topo const m_topo;

    //! surrounding box
    double m_box[5];

  public:
    /**
     * Oracle for a surface mesh with a surrounding box (top open) and topography.
     * For the top dimension the surface topography is considered as boundary,
     * for the remainder the box.
     *
     * Sketch:
     *
     *             __________________________
     *            |                          |
     *            |          outside         |
     *            |                          |
     *            |            *****         |   *
     *            |  **********     ******   | *   ***** <-- topography
     *            | *                     ***|*
     *          **|*                         |
     *  outside   |          inside          |  outside
     *            |                          |
     *            |__________________________|
     *
     *                       outside
     *
     * @param i_box boundaries of the box [box[0], box[1]] x [box[2], box[3]] x [box[4], infinity (topo)]
     * @param i_topoFile location of the topgraphic data.
     **/
    Oracle( double              i_box[5],
            std::string const & i_topoFile );

    CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT operator() (
      CGAL::Surface_mesh_default_triangulation_3::Geom_traits::Point_3 i_pt
    ) const;
};

#endif

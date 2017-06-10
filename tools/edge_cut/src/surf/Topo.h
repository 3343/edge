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
 * Meshing of topographic data.
 **/

#ifndef TOPO_H_
#define TOPO_H_

#include <string>

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>

#include "constants.hpp"

namespace edge_cut {
  namespace surf {
    class Topo;
  }
}

class edge_cut::surf::Topo {
  private:
    // 2.5D delaunay triangulation of the topographic data
    CGAL::Delaunay_triangulation_2 <
      CGAL::Projection_traits_xy_3 <
        CGAL::Cartesian<double>
      >
    > m_delTria;

    /*
     * Computes the 2.5D delaynay triangulation from the given topographic data.
     *
     * @param i_topoFile location of the topographic data.
     */
    void computeDelaunay( std::string const & i_topoFile );

  public:
    /**
     * Constructor: Derives the delaunay triangulation of the topo data.
     *
     * @param i_topoFile location of the topographic data.
     **/
    Topo( std::string const & i_topoFile );

    /**
     * Check if the given vertical ray intersects with the topo mesh.
     *
     * @param i_pt point where the ray starts.
     * @return true if an intersection was found, false otherwise.
     **/
    bool interRay( CGAL::Point_3< CGAL::Cartesian<double> > const & i_pt,
                   bool                                             i_positive=true ) const;

    /**
     * Intersects the given segment with the topo mesh.
     *
     * @param i_segPt1 first point of the segment.
     * @param i_segPt2 second point of the segment.
     * @param o_inters will be set to intersections of the segment and the topomesh. If more than one intersection occurs, multiple point will be returned.
     **/
    unsigned short interSeg( CGAL::Point_3< CGAL::Cartesian<double> > const & i_segPt1,
                             CGAL::Point_3< CGAL::Cartesian<double> > const & i_segPt2,
                             CGAL::Point_3< CGAL::Cartesian<double> >         o_inters[C_MAX_SURF_INTER] ) const;
};

#endif

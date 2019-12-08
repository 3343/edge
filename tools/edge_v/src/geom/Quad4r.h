/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Geometry computations for quad4r elements.
 **/
#ifndef EDGE_V_GEOM_QUAD4R_H
#define EDGE_V_GEOM_QUAD4R_H

#include <cstddef>

namespace edge_v {
  namespace geom {
    class Quad4r;
  }
}

class edge_v::geom::Quad4r {
  private:
    /**
     * Computes the maximum per-dimension distance of the vertices.
     *
     * @param i_veCrds vertex coordinates.
     * @param o_distMax will be set to maximum per-dimension distance.
     **/
    static void distMax( double const (*i_veCrds)[3],
                         double         o_distMax[2] );
  public:
    /**
     * Computes the area.
     *
     * @param i_veCrds vertex coordinates of the rectangle.
     * @return area of the rectangle.
     **/
    static double area( double const (*i_veCrds)[3] );

    /**
     * Computes the incircle-diameter.
     *
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameter( double const (*i_veCrds)[3] );

    /**
     * Normalizes the order of the vertices and faces.
     * The input is assumed to be ascending (w.r.t. the vertex ids) for all adjacency info.
     *
     * @param i_veCrds coordinates of the element's vertices.
     * @param io_elVe vertices adjacent to the element (ordered ascending by the ids).
     * @param io_elFa faces adjacent to the elements (ordered ascending by vertex ids).
     * @param io_elFaEl elements adjacent to elements, faces as bridge (ordered ascending by vertex ids of faces).
     **/
    static void normVesFas( double      const (* i_veCrds)[3],
                            std::size_t        * io_elVe,
                            std::size_t        * io_elFa,
                            std::size_t        * io_elFaEl );
};

#endif

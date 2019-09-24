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
 * Geometry computations for 3-node triangle elements.
 **/
#ifndef EDGE_V_MESH_GEOM_TRIA3_H
#define EDGE_V_MESH_GEOM_TRIA3_H

#include <cstdlib>

namespace edge_v {
  namespace mesh {
    class GeomTria3;
  }
}

class edge_v::mesh::GeomTria3 {
  public:
    /**
     * Computes the area.
     *
     * @param i_veCrds vertex coordinates of the triangle.
     * @return area of the triangle.
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
     * Example for the mapping of the initial, vertex-ordered face assignment to an orientation-based assignment:
     *
     *                 15                             15
     *                 *                              *
     *               *    *             map         *    *
     *             * 0    1  *         ---->      * 0     2 *
     *           *      2       *               *       1      *
     *       35 ******************** 80     35 ******************** 80
     *
     * Example for a valid mapping (0->1, 1->2, 2-> 0 is counter-clockwise):
     *
     *                0                   2 *
     *                *              Ref.   * *
     *             *      *         ---->   *   *
     *           *             *            *     *
     *       1 ******************** 2     0 ********* 1
     *
     * In this case we only check that the first face has vertices 0 and 1, the second face has vertices 1 and 2 and the
     * third face has vertices 2 and 0. No changes are performed.
     *
     * Example for an invalid mapping (0->1, 1->2, 2-> 0 is clockwise):
     *
     *                0                   1 *
     *                *              Ref.   * *
     *             *      *         ---->   *   *
     *           *             *            *     *
     *       2 ******************** 1     0 ********* 2
     *
     * In this case we exchange the local position of vertices 1 and 2.
     * Additional we change the ordering of the faces to match this.
     *
     * @param i_veCrds coordinates of the element's vertices.
     * @param io_elVe vertices adjacent to the element (ordered ascending by the ids).
     * @param io_elFa faces adjacent to the element (ordered ascending by the vertex ids of the faces).
     * @param io_elFaEl elements adjacent to the elements (ordered ascending by the vertex ids of the faces).
     **/
    static void normVesFas( double      const (* i_veCrds)[3],
                            std::size_t        * io_elVe,
                            std::size_t        * io_elFa,
                            std::size_t        * io_elFaEl );
};

#endif
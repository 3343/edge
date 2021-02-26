/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Geometry computations for 4-node tetrahedral elements.
 **/
#ifndef EDGE_V_GEOM_TET4_H
#define EDGE_V_GEOM_TET4_H

#include <cstdlib>
#include "../constants.h"

namespace edge_v {
  namespace geom {
    class Tet4;
  }
}

class edge_v::geom::Tet4 {
  private:
    //! tolerance used for zero-comparisons
    static constexpr double m_tol = 1E-6;

  public:
    /**
     * Computes the volume.
     *
     * @param i_veCrds vertex coordinates of the tetrahedron.
     * @return volume of the tetrahedron.
     **/
    static double volume(  double const (*i_veCrds)[3] );

    /**
     * Computes the insphere-diameter.
     *
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameter( double const (*i_veCrds)[3] );

    /**
     * Normalizes the order for faces and vertices.
     *
     * Given the tetrahedral vertices v0-v3. We keep the assignment of vertices v0 and v1.
     * Then we consider outer pointing normal of face v1-v2-v3.
     * Looking from the outside to this face (opposite direction of the normal),
     * we enforce counter-clockwise storage of vertices v2 and v3 w.r.t. to face 3.
     *
     * While exchanging v2 and v3 leaves the nodes (not the ordering) of face 2 (0-3-2)
     * and face 3 (1-2-3) untouched, we have to change the exchange the positions of face 0 and 1:
     *   face 0 (0-2-1) gets face 1 (0-1-3).
     *   face 1 (0-1-3) gets face 0 (0-2-2).
     *
     *
     * Example:
     *
     *                 v2                                          v3
     *                 *                                           *
     *               * . *                                       * . *
     *             *    .  *                                   *    .  *
     *           *       .   *                               *       .   *
     *         *          .    *                           *          .    *
     *       *  .            .   *              ->       *  .            .   *
     *  v1 *                    .  *                v1 *                    .  *
     *            *                . *                        *                . *
     *                    *           .*                              *           .*
     *                            *      *                                    *      *
     *                                     *                                            *
     *                                       v3                                          v2
     *
     * @param i_veCrds coordinates of the element's vertices.
     * @param io_elVe vertices adjacent to the element (ordered ascending by the ids).
     * @param io_elFa faces adjacent to the element (ordered ascending by the vertex ids of the faces).
     * @param io_elFaEl elements adjacent to the elements (ordered ascending by the vertex ids of the faces).
     **/
    static void normVesFas( double const (* i_veCrds)[3],
                            t_idx         * io_elVe,
                            t_idx         * io_elFa,
                            t_idx         * io_elFaEl );

    /**
     * Gets the adjacent elements' face-vertex ids for the given element-faces.
     * This is the local id of the vertex, which matches the faces' first vertex.
     *
     * In the case of periodic boundaries conditions, the vertex ids of the elements don't match for the respective face.
     *   In that case, one can provide the vertex coordinates through i_veCrds to trigger a
     *   coordinate-based search after a failed id-search.
     *
     * @param i_nFas number of element-faces to get ids for.
     * @param i_elOff element offset (added to element-ids in queries).
     * @param i_el elements to which the faces belong.
     * @param i_fa local face ids w.r.t. the elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_elFalEl elements adjacent to elements (faces as bridge).
     * @param i_veCrds coordinates of the vertices which will be used if id-search failed; ignored if nullptr is passed.
     * @param o_veIdsAd will be set to vertex ids w.r.t. to the adjacent elements.
     **/
    static void getVeIdsAd( t_idx                   i_nFas,
                            t_idx                   i_elOff,
                            t_idx          const  * i_el,
                            unsigned short const  * i_fa,
                            t_idx          const  * i_elVe,
                            t_idx          const  * i_elFaEl,
                            double         const (* i_veCrds)[3],
                            unsigned short        * o_veIdsAd );
};

#endif
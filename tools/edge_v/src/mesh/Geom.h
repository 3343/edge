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
 * Geometry computations for the mesh.
 **/
#ifndef EDGE_V_MESH_GEOM_H
#define EDGE_V_MESH_GEOM_H

#include "../constants.h"
#include <cmath>

namespace edge_v {
  namespace mesh {
    class Geom;
  }
}

class edge_v::mesh::Geom {
  private:
    /**
     * Computes the length of a line.
     *
     * @param i_veCrds vertex coordinates of the line.
     * @return length of the line.
     **/
    static double lengthLine( double const (*i_veCrds)[3] );

    /**
     * Computes the area of a 3-node triangle.
     *
     * @param i_veCrds vertex coordinates of the triangle.
     * @return area of the triangle.
     **/
    static double areaTria3( double const (*i_veCrds)[3] );

    /**
     * Computes the volume of a 4-node tetrahedron.
     *
     * @param i_veCrds vertex coordinates of the tetrahedron.
     * @return volume of the tetrahedron.
     **/
    static double volumeTet4(  double const (*i_veCrds)[3] );

    /**
     * Computes the incircle-diameter of 3-node triangles.
     *
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameterTria3( double const (*i_veCrds)[3] );

    /**
     * Computes the insphere-diameter of 4-node tetrahedrons.
     *
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameterTet4( double const (*i_veCrds)[3] );

    /**
     * Computes the directed normal for a line.
     *
     * @param i_veCrds vertices of the line.
     * @param i_nPt point lying on the non-normal side of the line (normal will point away from this point).
     * @param o_normal will be set to normal.
     **/
    static void normalLine( double const (*i_veCrds)[3],
                            double const   i_nPt[3],
                            double         o_normal[3] );

    /**
     * Computes the tangent for a line (right-handed coordinate system with normal).
     *
     * @param i_veCrds vertices of the line.
     * @param i_nPt point lying on the non-normal side of the line (normal is assumed to point away from this point).
     * @param o_tangent will be set to tangent.
     **/
    static void tangentLine( double const (*i_veCrds)[3],
                             double const   i_nPt[3],
                             double         o_tangent[3] );

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
    static void normVesFasTria3( double      const (* i_veCrds)[3],
                                 std::size_t        * io_elVe,
                                 std::size_t        * io_elFa,
                                 std::size_t        * io_elFaEl );

    /**
     * Normalizes the order for faces and vertices of tetrahedral elements.
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
    static void normVesFasTet4( double      const (* i_veCrds)[3],
                                std::size_t        * io_elVe,
                                std::size_t        * io_elFa,
                                std::size_t        * io_elFaEl );

  public:
    /**
     * Computes the volume of the given entity.
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @return volume of the entity.
     **/
    static double volume( t_entityType         i_enTy,
                          double       const (*i_veCrds)[3] );

    /**
     * Computes the normal for the given entity
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @param i_nPt normal point on one side of the entity (normal will point to the other direction).
     * @param o_normal will be set to normal.
     **/
    static void normal( t_entityType         i_enTy,
                        double       const (*i_veCrds)[3],
                        double       const   i_nPt[3],
                        double               o_normal[3] );

    /**
     * Computes the length (1d), incircle (2d) or insphere diameter (3d).
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameter( t_entityType         i_enTy,
                              double       const (*i_veCrds)[3] );

    /**
     * Normalizes the order of the vertices and faces.
     * The input is assumed to be ascending (w.r.t. the vertex ids) for all adjacency info.
     *
     * @param i_elTy element type.
     * @param i_veCrds coordinates of the element's vertices.
     * @param io_elVe vertices adjacent to the element (ordered ascending by the ids).
     * @param io_elFa faces adjacent to the element (ordered ascending by the vertex ids of the faces).
     * @param io_elFaEl elements adjacent to the elements (ordered ascending by the vertex ids of the faces).
     **/
    static void normVesFas( t_entityType         i_elTy,
                            double      const (* i_veCrds)[3],
                            std::size_t        * io_elVe,
                            std::size_t        * io_elFa,
                            std::size_t        * io_elFaEl );
};

#endif
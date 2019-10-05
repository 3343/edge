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
#ifndef EDGE_V_GEOM_GEOM_H
#define EDGE_V_GEOM_GEOM_H

#include "../constants.h"
#include <cmath>

namespace edge_v {
  namespace geom {
    class Geom;
  }
}

class edge_v::geom::Geom {
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
     * Computes the tangents for the given entity
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @param i_nPt normal point on one side of the entity (normal is assumed to point in the other direction).
     * @param o_tangents will be set to tangents.
     **/
    static void tangents( t_entityType         i_enTy,
                          double       const (*i_veCrds)[3],
                          double       const   i_nPt[3],
                          double               o_tangents[2][3] );

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
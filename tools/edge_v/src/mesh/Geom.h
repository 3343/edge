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

#include "constants.h"
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
     * Computes the length (1d), incircle (2d) or insphere diameter (3d).
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @return diameter.
     **/
    static double inDiameter( t_entityType         i_enTy,
                              double       const (*i_veCrds)[3] );
};

#endif
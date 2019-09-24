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
 * Geometry computations for line elements.
 **/
#ifndef EDGE_V_MESH_GEOM_LINE_H
#define EDGE_V_MESH_GEOM_LINE_H

namespace edge_v {
  namespace mesh {
    class GeomLine;
  }
}

class edge_v::mesh::GeomLine {
  public:
    /**
     * Computes the length.
     *
     * @param i_veCrds vertex coordinates of the line.
     * @return length of the line.
     **/
    static double length( double const (*i_veCrds)[3] );

    /**
     * Computes the directed normal.
     *
     * @param i_veCrds vertices of the line.
     * @param i_nPt point lying on the non-normal side of the line (normal will point away from this point).
     * @param o_normal will be set to normal.
     **/
    static void normal( double const (*i_veCrds)[3],
                        double const   i_nPt[3],
                        double         o_normal[3] );

    /**
     * Computes the tangent for (right-handed coordinate system with normal).
     *
     * @param i_veCrds vertices of the line.
     * @param i_nPt point lying on the non-normal side of the line (normal is assumed to point away from this point).
     * @param o_tangent will be set to tangent.
     **/
    static void tangent( double const (*i_veCrds)[3],
                         double const   i_nPt[3],
                         double         o_tangent[3] );
};

#endif
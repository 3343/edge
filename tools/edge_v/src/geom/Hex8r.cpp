/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Geometry computations for hex8r elements.
 **/
#include "Hex8r.h"
#include "Generic.h"
#include <cmath>
#include "../io/logging.h"

void edge_v::geom::Hex8r::distMax( double const (*i_veCrds)[3],
                                   double         o_distMax[3] ) {
  o_distMax[0] = o_distMax[1] = o_distMax[2] = 0;

  for( unsigned short l_v0 = 0; l_v0 < 8; l_v0++ ) {
    for( unsigned short l_v1 = 0; l_v1 < 8; l_v1++ ) {
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        double l_dist = std::abs( i_veCrds[l_v0][l_di] - i_veCrds[l_v1][l_di] );
        o_distMax[l_di] = std::max( o_distMax[l_di], l_dist );
      }
    }
  }
  EDGE_V_CHECK_NE( o_distMax[0], 0.0 );
  EDGE_V_CHECK_NE( o_distMax[1], 0.0 );
  EDGE_V_CHECK_NE( o_distMax[2], 0.0 );
}

double edge_v::geom::Hex8r::volume( double const (*i_veCrds)[3] ) {
  double l_distMax[3] = {0};
  distMax( i_veCrds,
           l_distMax );
  double l_volume = l_distMax[0] * l_distMax[1] * l_distMax[2];

  return l_volume;
}

double edge_v::geom::Hex8r::inDiameter( double const (*i_veCrds)[3] ) {
  double l_distMax[3] = {0};
  distMax( i_veCrds,
           l_distMax );
  double l_dia = std::min( l_distMax[0], l_distMax[1] );
  l_dia = std::min( l_dia, l_distMax[2] );

  return l_dia;
}

void edge_v::geom::Hex8r::normVesFas( double const (* i_veCrds)[3],
                                      t_idx         * io_elVe,
                                      t_idx         * io_elFa,
                                      t_idx         * io_elFaEl ) {
  // lambda which compares two vertices
  auto l_veLess = [ i_veCrds ]( unsigned short i_ve0,
                                unsigned short i_ve1 ) {
    if( i_veCrds[i_ve1][2] - i_veCrds[i_ve0][2] > m_tol ) {
      return true;
    }
    else if(    std::abs( i_veCrds[i_ve1][2] - i_veCrds[i_ve0][2] ) < m_tol
             && i_veCrds[i_ve1][1] - i_veCrds[i_ve0][1] > m_tol ) {
      return true;
    }
    else if(    std::abs( i_veCrds[i_ve1][2] - i_veCrds[i_ve0][2] ) < m_tol
             && std::abs( i_veCrds[i_ve1][1] - i_veCrds[i_ve0][1] ) < m_tol
             && i_veCrds[i_ve1][0] - i_veCrds[i_ve0][0] > m_tol ) {
      return true;
    }
    return false;
  };

  unsigned short l_sorted[8] = {0};
  for( unsigned short l_ve = 0; l_ve < 8; l_ve++ ) {
    l_sorted[l_ve] = l_ve;
  }

  // sort the vertices
  std::sort( l_sorted,
             l_sorted+8,
             l_veLess );

  /* EDGE's hex8r reference element:
   *
   *   face 0: 0-3-2-1
   *   face 1: 0-1-5-4
   *   face 2: 1-2-6-5
   *   face 3: 3-7-6-2
   *   face 4: 0-4-7-3
   *   face 5: 4-5-6-7
   * we are forcing counter-clockwise storage of face vertices w.r.t. to opposite face.
   *
   *
   *           7 x*******************x 6
   *            **                  **
   *           * *                 * *
   *          *  *                *  *
   *         *   *               *   *
   *        *    *              *    *
   *     4 x*******************x 5   *
   *       *     *             *     *
   *       *   3 x************ * ****x 2
   *       *    *              *    *
   *       *   *               *   *
   *       |  /                *  *
   *  zeta | / eta             * *
   *       |/                  **
   *       x---****************x
   *     0   xi                 1
   */

  // adjust sorting to match that of EDGE
  unsigned short l_tmp = l_sorted[2];
  l_sorted[2] = l_sorted[3];
  l_sorted[3] = l_tmp;

  l_tmp = l_sorted[6];
  l_sorted[6] = l_sorted[7];
  l_sorted[7] = l_tmp;

  // perform the element-vertex reordering
  t_idx l_elVeTmp[8] = {0};
  for( unsigned short l_ve = 0; l_ve < 8; l_ve++ ) {
    l_elVeTmp[l_ve] = io_elVe[l_ve];
  }
  for( unsigned short l_ve = 0; l_ve < 8; l_ve++ ) {
    io_elVe[l_ve] = l_elVeTmp[ l_sorted[l_ve] ];
  }

  // get vertices of faces
  t_idx l_faVe[6][4] = { {0, 3, 2, 1},
                         {0, 1, 5, 4},
                         {1, 2, 6, 5},
                         {3, 7, 6, 2},
                         {0, 4, 7, 3},
                         {4, 5, 6, 7} };
  for( unsigned short l_fa = 0; l_fa < 6; l_fa++ ) {
    for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
      l_faVe[l_fa][l_ve] = io_elVe[ l_faVe[l_fa][l_ve] ];
    }
    std::sort( l_faVe[l_fa], l_faVe[l_fa]+4 );
  }

  // derive mapping from current vertex-ordering of faces to that of the input data
  t_idx l_mapping[6] = {0};
  Generic::sortLex( 6,
                    4,
                    l_faVe[0],
                    l_mapping );

  // reorder face-data structures
  t_idx l_elFaTmp[6] = {0};
  t_idx l_elFaElTmp[6] = {0};
  for( unsigned short l_fa = 0; l_fa < 6; l_fa++ ) {
    l_elFaTmp[l_fa] = io_elFa[l_fa];
    l_elFaElTmp[l_fa] = io_elFaEl[l_fa];
  }
  for( unsigned short l_fa = 0; l_fa < 6; l_fa++ ) {
    t_idx l_id = l_mapping[l_fa];
    io_elFa[l_id] = l_elFaTmp[l_fa];
    io_elFaEl[l_id] = l_elFaElTmp[l_fa];
  }
}
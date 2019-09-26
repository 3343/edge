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
#include <limits>
#include <Eigen/Dense>
#include "io/logging.h"
#include "GeomTria3.h"

double edge_v::mesh::GeomTria3::area( double const (*i_veCrds)[3] ) {
  // get edges
  Eigen::Vector3d l_a, l_b;
  l_a[0] = i_veCrds[1][0] - i_veCrds[0][0];
  l_a[1] = i_veCrds[1][1] - i_veCrds[0][1];
  l_a[2] = i_veCrds[1][2] - i_veCrds[0][2];

  l_b[0] = i_veCrds[2][0] - i_veCrds[0][0];
  l_b[1] = i_veCrds[2][1] - i_veCrds[0][1];
  l_b[2] = i_veCrds[2][2] - i_veCrds[0][2];

  // compute cross product
  Eigen::Vector3d l_cross = l_a.cross( l_b );

  // compute volume
  double l_vol = l_cross.norm() / 2;

  return l_vol;
}

void edge_v::mesh::GeomTria3::tangents( double const (*i_veCrds)[3],
                                        double const   i_nPt[3],
                                        double         o_tangents[2][3] ) {
  /*
   * we can derive the second (normal) tangent from the first tangent
   *
   *            x <- d1
   * d0 -> x   *
   *       *  *
   *       * *
   *       x******x <- n: that's the vector we want
   *
   *  <d0, n> = 0
   *  <d0, d1 + alpha * d0> = 0
   *  <d0, d1> + <d0, alpha * d0> = 0
   *  => alpha = - <d0,d1> / <d0,d0>
   *  => n = d1 + alpha * d0
   */
  Eigen::Vector3d l_v0( i_veCrds[0] );
  Eigen::Vector3d l_v1( i_veCrds[1] );
  Eigen::Vector3d l_v2( i_veCrds[2] );

  // directions
  Eigen::Vector3d l_d0 = l_v1 - l_v0;
  Eigen::Vector3d l_d1 = l_v2 - l_v0;

  // compute alpha
  double l_alpha = -l_d0.dot( l_d1 ) / l_d0.dot( l_d0 );

  // compute in-plane normal to d0 and normalize
  Eigen::Vector3d l_n0 = l_d1 + l_alpha * l_d0;
  l_d0.normalize();
  l_n0.normalize();

  // compute normal to d0 and n0
  Eigen::Vector3d l_n1 = l_d0.cross( l_n0 );

  // check if we have to change the orientation
  // which implies a change of the tangents right-handed coordinate system
  Eigen::Vector3d l_np( i_nPt );
  Eigen::Vector3d l_dn = l_np - l_v0;
  double l_dp = l_n1.dot( l_dn );
  if( l_dp > 0 ) {
    Eigen::Vector3d l_tmp = l_d0;
    l_d0 = l_n0;
    l_n0 = l_tmp;
  }

  // store results
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    o_tangents[0][l_di] = l_d0[l_di];
    o_tangents[1][l_di] = l_n0[l_di];
  }
}

void edge_v::mesh::GeomTria3::normal( double const (*i_veCrds)[3],
                                      double const   i_nPt[3],
                                      double         o_normal[3] ) {
  // compute tangents
  double l_ts[2][3] = {0};
  tangents( i_veCrds,
            i_nPt,
            l_ts );

  Eigen::Vector3d l_t0( l_ts[0] );
  Eigen::Vector3d l_t1( l_ts[1] );

  // derive normal as cross-product
  Eigen::Vector3d l_n = l_t0.cross( l_t1 );

  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    o_normal[l_di] = l_n[l_di];
}

double edge_v::mesh::GeomTria3::inDiameter( double const (*i_veCrds)[3] ) {
  // check zero last dimension
  EDGE_V_CHECK_EQ( i_veCrds[0][2], 0 );
  EDGE_V_CHECK_EQ( i_veCrds[1][2], 0 );
  EDGE_V_CHECK_EQ( i_veCrds[2][2], 0 );

  // extract points
  Eigen::Vector3d l_a( i_veCrds[0] );
  Eigen::Vector3d l_b( i_veCrds[1] );
  Eigen::Vector3d l_c( i_veCrds[2] );

  // compute lengths
  double l_lengths[3];
  l_lengths[0] = ( l_a - l_b ).norm();
  l_lengths[1] = ( l_b - l_c ).norm();
  l_lengths[2] = ( l_c - l_a ).norm();

  // compute the area of the triangle
  double l_area = area( i_veCrds );

  // compute the diameter
  double l_dia = 4 * l_area;
         l_dia /= l_lengths[0] + l_lengths[1] + l_lengths[2];

  return l_dia;
}

void edge_v::mesh::GeomTria3::normVesFas( double      const (* i_veCrds)[3],
                                          std::size_t        * io_elVe,
                                          std::size_t        * io_elFa,
                                          std::size_t        * io_elFaEl ) {
  /*
   * reorder face-information in ascending order of vertex ids:
   * id0--fa0-->id1--fa1-->id2--fa2-->
   */
  std::size_t l_faTmp = io_elFa[2];
  io_elFa[2] = io_elFa[1];
  io_elFa[1] = l_faTmp;

  std::size_t l_elTmp = io_elFaEl[2];
  io_elFaEl[2] = io_elFaEl[1];
  io_elFaEl[1] = l_elTmp;

  // get vectors point from 0->1 and 0->2
  Eigen::Matrix2d l_m;
  l_m(0, 0) = i_veCrds[1][0] - i_veCrds[0][0];
  l_m(1, 0) = i_veCrds[1][1] - i_veCrds[0][1];

  l_m(0, 1) = i_veCrds[2][0] - i_veCrds[0][0];
  l_m(1, 1) = i_veCrds[2][1] - i_veCrds[0][1];

  // check if we have to reorder based on the determinant
  double l_det = l_m.determinant();

  // negative determinant -> clockwise -> change pos of 2nd and 3rd vertex
  if( l_det < 0 ) {
    std::size_t l_veTmp = io_elVe[1];
    io_elVe[1] = io_elVe[2];
    io_elVe[2] = l_veTmp;

    /*
     * change position of face-information accordingly
     *
     *              v0                                v0
     *               *                                *
     *          f2 *     * f0      ---->        f0  *     * f2
     *           *    f1    *                    *      f1    *
     *      v2 ****************** v1         v1 ****************** v2
     *
     * Swapping positions of vertices v1 and v2 means that we have to swap positions
     * of faces f0 and f2 also.
     */
    l_faTmp = io_elFa[0];
    io_elFa[0] = io_elFa[2];
    io_elFa[2] = l_faTmp;

    l_elTmp = io_elFaEl[0];
    io_elFaEl[0] = io_elFaEl[2];
    io_elFaEl[2] = l_elTmp;
  }
}
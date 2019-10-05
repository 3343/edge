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
 * Geometry computations for 4-node tetrahedral elements.
 **/
#include "Tet4.h"

#include <limits>
#include "io/logging.h"
#include <Eigen/Dense>

double edge_v::geom::Tet4::volume( double const (*i_veCrds)[3] ) {

  // assembly matrix
  Eigen::Matrix4d l_mat;

  for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      l_mat(l_ve, l_di) = i_veCrds[l_ve][l_di];
    }
    l_mat(l_ve, 3) = 1.0;
  }

  // compute volume
  double l_vol = std::abs( l_mat.determinant() ) / 6;

  return l_vol;
}

double edge_v::geom::Tet4::inDiameter( double const (*i_veCrds)[3] ) {
  /*
   * Reference: John Burkardt
   *            The Inscribed Sphere of a Tetrahedron
   *            Computational Geometry Lab: TETRAHEDRONS
   *            2010
   */
  // define lambda for computation of normal vector (faces) length
  auto l_faNorm = []( Eigen::Vector3d const i_a,  Eigen::Vector3d const i_b, Eigen::Vector3d const i_c ) {
    Eigen::Vector3d l_cross = (i_b-i_a).cross( (i_c-i_a) );
    return l_cross.norm();
  };

  // construct vectors of the vertices
  Eigen::Vector3d l_a( i_veCrds[0] );
  Eigen::Vector3d l_b( i_veCrds[1] );
  Eigen::Vector3d l_c( i_veCrds[2] );
  Eigen::Vector3d l_d( i_veCrds[3] );

  double l_div  = l_faNorm( l_a, l_b, l_c );
         l_div += l_faNorm( l_a, l_b, l_d );
         l_div += l_faNorm( l_a, l_c, l_d );
         l_div += l_faNorm( l_b, l_c, l_d );

  double l_vol = volume( i_veCrds );

  double l_dia = (12 * l_vol) / l_div;

  return l_dia;
}

void edge_v::geom::Tet4::normVesFas( double      const (* i_veCrds)[3],
                                     std::size_t        * io_elVe,
                                     std::size_t        * io_elFa,
                                     std::size_t        * io_elFaEl ) {
  // get vectors point from 0->1, 0->2 and 0->3
  Eigen::Matrix3d l_m;
  for( unsigned short l_d0 = 0; l_d0 < 3; l_d0++ )
    for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ )
      l_m(l_d0, l_d1) = i_veCrds[l_d1+1][l_d0] - i_veCrds[0][l_d0];

  // check if we have to reorder based on the determinant
  double l_det = l_m.determinant();

  // assert non-planar vertices
  EDGE_V_CHECK_GT( std::abs(l_det), 1E-5 );

  // negative determinant -> clockwise -> exchange vertices 2,3 and faces 0,1
  if( l_det < 0 ) {
    // exchange vertices
    std::size_t l_tmpVe = io_elVe[2];
    io_elVe[2] = io_elVe[3];
    io_elVe[3] = l_tmpVe;

    // exchange faces
    std::size_t l_tmpFa = io_elFa[0];
    io_elFa[0] = io_elFa[1];
    io_elFa[1] = l_tmpFa;

    std::size_t l_tmpEl = io_elFaEl[0];
    io_elFaEl[0] = io_elFaEl[1];
    io_elFaEl[1] = l_tmpEl;
  }
}
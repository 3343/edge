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
#include "Geom.h"

#include <limits>
#include "io/logging.h"
#include <Eigen/Dense>

double edge_v::mesh::Geom::lengthLine( double const (*i_veCrds)[3] ) {
  double l_len = 0;

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    double l_diff = i_veCrds[1][l_di] - i_veCrds[0][l_di];
    l_len += l_diff * l_diff;
  }
  l_len = std::sqrt( l_len );

  return l_len;
}

double edge_v::mesh::Geom::areaTria3( double const (*i_veCrds)[3] ) {
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

double edge_v::mesh::Geom::volumeTet4( double const (*i_veCrds)[3] ) {

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

double edge_v::mesh::Geom::volume( t_entityType         i_enTy,
                                   double       const (*i_veCrds)[3] ) {
  double l_vol = std::numeric_limits< double >::max();

  if( i_enTy == LINE ) {
    l_vol = lengthLine( i_veCrds );
  }
  else if( i_enTy == TRIA3 ) {
    l_vol = areaTria3( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_vol = volumeTet4( i_veCrds );
  }
  else EDGE_V_LOG_FATAL;

  return l_vol;
}

double edge_v::mesh::Geom::inDiameterTria3( double const (*i_veCrds)[3] ) {
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
  double l_area = areaTria3( i_veCrds );

  // compute the diameter
  double l_dia = 4 * l_area;
         l_dia /= l_lengths[0] + l_lengths[1] + l_lengths[2];

  return l_dia;
}

double edge_v::mesh::Geom::inDiameterTet4( double const (*i_veCrds)[3] ) {
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

  double l_vol = volumeTet4( i_veCrds );

  double l_dia = (12 * l_vol) / l_div;

  return l_dia;
}

double edge_v::mesh::Geom::inDiameter( t_entityType         i_enTy,
                                       double       const (*i_veCrds)[3] ) {
  double l_dia = std::numeric_limits< double >::max();
  if( i_enTy == TRIA3 ) {
    l_dia = inDiameterTria3( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_dia = inDiameterTet4( i_veCrds );
  }
  else EDGE_V_LOG_FATAL;

  return l_dia;
}
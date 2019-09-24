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
#include "GeomLine.h"

#include <limits>
#include "io/logging.h"
#include <Eigen/Dense>

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
    l_vol = GeomLine::length( i_veCrds );
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

void edge_v::mesh::Geom::normVesFasTria3( double      const (* i_veCrds)[3],
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

void edge_v::mesh::Geom::normVesFasTet4( double      const (* i_veCrds)[3],
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

void edge_v::mesh::Geom::normal( t_entityType         i_enTy,
                                 double       const (*i_veCrds)[3],
                                 double       const   i_nPt[3],
                                 double               o_normal[3] ) {
  if( i_enTy == LINE ) {
    GeomLine::normal( i_veCrds,
                      i_nPt,
                      o_normal );
  }
  else EDGE_V_LOG_FATAL;
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

void edge_v::mesh::Geom::normVesFas( t_entityType         i_elTy,
                                     double      const (* i_veCrds)[3],
                                     std::size_t        * io_elVe,
                                     std::size_t        * io_elFa,
                                     std::size_t        * io_elFaEl ) {
  if( i_elTy == TRIA3 ) {
    normVesFasTria3( i_veCrds,
                     io_elVe,
                     io_elFa,
                     io_elFaEl );
  }
  else if( i_elTy == TET4 ) {
    normVesFasTet4( i_veCrds,
                    io_elVe,
                    io_elFa,
                    io_elFaEl );
  }
  else EDGE_V_LOG_FATAL;
}
/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
#include "Generic.h"

#include <limits>
#include "../io/logging.h"
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

void edge_v::geom::Tet4::normVesFas( double const (* i_veCrds)[3],
                                     t_idx         * io_elVe,
                                     t_idx         * io_elFa,
                                     t_idx         * io_elFaEl ) {
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
    t_idx l_tmpVe = io_elVe[2];
    io_elVe[2] = io_elVe[3];
    io_elVe[3] = l_tmpVe;

    // exchange faces
    t_idx l_tmpFa = io_elFa[0];
    io_elFa[0] = io_elFa[1];
    io_elFa[1] = l_tmpFa;

    t_idx l_tmpEl = io_elFaEl[0];
    io_elFaEl[0] = io_elFaEl[1];
    io_elFaEl[1] = l_tmpEl;
  }
}

void edge_v::geom::Tet4::getVeIdsAd( t_idx                   i_nFas,
                                     t_idx                   i_elOff,
                                     t_idx          const  * i_el,
                                     unsigned short const  * i_fa,
                                     t_idx          const  * i_elVe,
                                     t_idx          const  * i_elFaEl,
                                     double         const (* i_veCrds)[3],
                                     unsigned short        * o_veIdsAd ) {
  for( t_idx l_id = 0; l_id < i_nFas; l_id++ ) {
    // get element and face id
    t_idx l_el = i_el[l_id] + i_elOff;
    unsigned short l_fa = i_fa[l_id];

    // adjacent element
    t_idx l_elAd = i_elFaEl[l_el*4 + l_fa];

    // vertex at position zero of the face
    t_idx l_ve0 = (l_fa <= 2) ? i_elVe[l_el*4+0] : i_elVe[l_el*4+1];

    // get the respective position in the adjacent element
    t_idx l_veAd = std::numeric_limits< t_idx >::max();
    for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
      if( i_elVe[l_elAd*4+l_ve] == l_ve0 ) l_veAd = l_ve;
    }
    // go through coordinates if available (periodic boundaries)
    if( i_veCrds != nullptr && l_veAd == std::numeric_limits< t_idx >::max() ) {
      for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
        t_idx l_veId = i_elVe[l_elAd*4+l_ve];

        unsigned short l_nEq = 0;
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          double l_diff = i_veCrds[l_ve0][l_di] - i_veCrds[l_veId][l_di];
          if( std::abs(l_diff) < 1E-5 ) l_nEq++;
        }
        if( l_nEq == 2 ) {
          EDGE_V_CHECK_EQ( l_veAd, std::numeric_limits< t_idx >::max() );
          l_veAd = l_ve;
        }
      }
    }
    EDGE_V_CHECK_LT( l_veAd, 4 );

    // get the face id
    unsigned short l_faId = std::numeric_limits< unsigned short >::max();
    Generic::getFaIdsAd( TET4,
                         1,
                         0,
                         &l_el,
                         &l_fa,
                         i_elFaEl,
                         &l_faId );
    EDGE_V_CHECK_LT( l_faId, 4 );

    /*
      * decide, depending on the neighboring face, what vertex combination we have
      *
      * Example:
      *
      *      face 0                   face 1
      *         3                       2
      *         *                         *
      *       *    *         neighbors    *  *
      *     *        *                    *     *
      * 0  ************* 1              0 ********** 3 <-- dominant 0 goes here
      *
      * -> We have vertex combi 1 out of possible 0-2.
      */
      if( l_faId == 0 ) {
        if(      l_veAd == 0 ) o_veIdsAd[l_id] = 0;
        else if( l_veAd == 2 ) o_veIdsAd[l_id] = 1;
        else if( l_veAd == 1 ) o_veIdsAd[l_id] = 2;
        else EDGE_V_LOG_FATAL;
      }
      else if( l_faId == 1 ) {
        if(      l_veAd == 0 ) o_veIdsAd[l_id] = 0;
        else if( l_veAd == 1 ) o_veIdsAd[l_id] = 1;
        else if( l_veAd == 3 ) o_veIdsAd[l_id] = 2;
        else EDGE_V_LOG_FATAL;
      }
      else if( l_faId == 2 ) {
        if(      l_veAd == 0 ) o_veIdsAd[l_id] = 0;
        else if( l_veAd == 3 ) o_veIdsAd[l_id] = 1;
        else if( l_veAd == 2 ) o_veIdsAd[l_id] = 2;
        else EDGE_V_LOG_FATAL;
      }
      else if( l_faId == 3 ) {
        if(      l_veAd == 1 ) o_veIdsAd[l_id] = 0;
        else if( l_veAd == 2 ) o_veIdsAd[l_id] = 1;
        else if( l_veAd == 3 ) o_veIdsAd[l_id] = 2;
        else EDGE_V_LOG_FATAL;
      }
      else EDGE_V_LOG_FATAL;
  }
}
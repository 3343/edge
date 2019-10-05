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
#include "GeomLine.h"

#include <Eigen/Dense>

double edge_v::geom::Line::length( double const (*i_veCrds)[3] ) {
  double l_len = 0;

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    double l_diff = i_veCrds[1][l_di] - i_veCrds[0][l_di];
    l_len += l_diff * l_diff;
  }
  l_len = std::sqrt( l_len );

  return l_len;
}

void edge_v::geom::Line::normal( double const (*i_veCrds)[3],
                                 double const   i_nPt[3],
                                 double         o_normal[3] ) {
  // unit vector: v0 -> v1
  Eigen::Vector2d l_v0( i_veCrds[0] );
  Eigen::Vector2d l_v1( i_veCrds[1] );

  Eigen::Vector2d l_d = l_v1 - l_v0;
  l_d.normalize();

  // derive normal
  Eigen::Vector2d l_n;
  l_n[0] = -l_d[1];
  l_n[1] =  l_d[0];

  /*
   * compute dot product of vector pointing from a vertex to the normal point and the normal
   *
   *            x * * * * * o <-- normal point
   *          *  *     .   o
   *        *     *  .  <-----example normal
   *      x     1  *     o <-- vector pointing from a vertex
   *         *      *   o
   *             *   * o
   *                * x
   *  translates to:
   *  |        o
   *  |       o \
   *  |      o   \
   *  |     o   a . <--- projection (a=90deg) by the dot product > 0:
   *  |    o    .        we are on the wrong side here.
   *  |   o   .
   *  |  o  .
   *  | o .
   *  |o.
   *  |_______
   */
  Eigen::Vector2d l_nPt( i_nPt );
  Eigen::Vector2d l_dn = l_nPt - l_v0;

  double l_dp = l_n.dot( l_dn );

  // if the dot product is positive, the angle is below 90deg:
  // we want the normal to point in the other direction; therefore we have to change the sign
  if( l_dp > 0 ) {
    l_n *= -1;
  }

  o_normal[0] = l_n[0];
  o_normal[1] = l_n[1];
  o_normal[2] = 0;
}

void edge_v::geom::Line::tangent( double const (*i_veCrds)[3],
                                  double const   i_nPt[3],
                                  double         o_tangent[3] ) {
  double l_n[3] = {0};
  normal( i_veCrds, i_nPt, l_n );

  o_tangent[0] = -l_n[1];
  o_tangent[1] =  l_n[0];
  o_tangent[2] =  0;
}
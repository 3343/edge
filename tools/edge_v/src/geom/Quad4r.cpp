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
 * Geometry computations for quad4r elements.
 **/
#include "Quad4r.h"
#include "Tria3.h"
#include <cmath>
#include <limits>
#include "../io/logging.h"

void edge_v::geom::Quad4r::distMax( double const (*i_veCrds)[3],
                                    double         o_distMax[2] ) {
  o_distMax[0] = o_distMax[1] = 0;

  for( unsigned short l_v0 = 0; l_v0 < 4; l_v0++ ) {
    for( unsigned short l_v1 = 0; l_v1 < 4; l_v1++ ) {
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        double l_dist = std::abs( i_veCrds[l_v0][l_di] - i_veCrds[l_v1][l_di] );
        o_distMax[l_di] = std::max( o_distMax[l_di], l_dist );
      }
    }
  }
  EDGE_V_CHECK_NE( o_distMax[0], 0.0 );
  EDGE_V_CHECK_NE( o_distMax[1], 0.0 );
}

double edge_v::geom::Quad4r::area( double const (*i_veCrds)[3] ) {
  // split quad4r in equal-area triangles
  double l_area = Tria3::area( i_veCrds );
  l_area *= 2;

  return l_area;
}

void edge_v::geom::Quad4r::tangents( double const (*i_veCrds)[3],
                                     double const   i_nPt[3],
                                     double         o_tangents[2][3] ) {
  // quad4r elements are planar -> ignore a vertex
  Tria3::tangents( i_veCrds,
                    i_nPt,
                    o_tangents );
}

void edge_v::geom::Quad4r::normal( double const (*i_veCrds)[3],
                                   double const   i_nPt[3],
                                   double         o_normal[3] ) {
  // quad4r elements are planar -> ignore a vertex
  Tria3::normal( i_veCrds,
                 i_nPt,
                 o_normal );
}

double edge_v::geom::Quad4r::inDiameter( double const (*i_veCrds)[3] ) {
  double l_distMax[2] = {0};
  distMax( i_veCrds,
           l_distMax );
  double l_dia = std::min( l_distMax[0], l_distMax[1] );

  return l_dia;
}

void edge_v::geom::Quad4r::normVesFas( double const (* i_veCrds)[3],
                                       t_idx         * io_elVe,
                                       t_idx         * io_elFa,
                                       t_idx         * io_elFaEl ) {
  // store coordinates of element's vertices for reordering
  double l_veCrds[4][2];
  for( unsigned short l_ve = 0; l_ve < 4; l_ve++ )
    for( unsigned short l_di = 0; l_di < 2; l_di++ )
      l_veCrds[l_ve][l_di] = i_veCrds[l_ve][l_di];

  // local reordering
  unsigned short l_elVe[4] = {0, 1, 2, 3};

  // sort by y-coord
  for( unsigned short l_v0 = 0; l_v0 < 4; l_v0++ ) {
    for( unsigned short l_v1 = l_v0; l_v1 < 4; l_v1++ ) {
      if( l_veCrds[l_v1][1] < l_veCrds[l_v0][1] ) {
        t_idx l_veTmp = l_elVe[l_v1];
        l_elVe[l_v1] = l_elVe[l_v0];
        l_elVe[l_v0] = l_veTmp;

        for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
          double l_crdTmp = l_veCrds[l_v1][l_di];
          l_veCrds[l_v1][l_di] = l_veCrds[l_v0][l_di];
          l_veCrds[l_v0][l_di] = l_crdTmp; 
        }
      }
    }
  }

  // sort by x
  if( l_veCrds[0][0] > l_veCrds[1][0] ) {
    t_idx l_veTmp = l_elVe[0];
    l_elVe[0] = l_elVe[1];
    l_elVe[1] = l_veTmp;
  }

  if( l_veCrds[3][0] > l_veCrds[2][0] ) {
    t_idx l_veTmp = l_elVe[3];
    l_elVe[3] = l_elVe[2];
    l_elVe[2] = l_veTmp;
  }

  // derive location of faces in old lexicographical order
  unsigned short l_faOr[4] = { 3, 3, 3, 3 };
  unsigned short l_vePos[2] = { std::numeric_limits< unsigned short >::max(),
                                std::numeric_limits< unsigned short >::max() };
  for( unsigned short l_ve = 0; l_ve < 4; l_ve++ )
    if( l_elVe[l_ve] == 0 ) l_vePos[0] = l_ve;
  for( unsigned short l_ve = 0; l_ve < 4; l_ve++ )
    if( l_elVe[l_ve] == 1 ) l_vePos[1] = l_ve;

  if( l_elVe[ (l_vePos[0]+3)%4 ] > l_elVe[ (l_vePos[0]+1)%4 ] ) {
    l_faOr[ l_vePos[0]       ] = 0;
    l_faOr[ (l_vePos[0]+3)%4 ] = 1;
  }
  else {
    l_faOr[ l_vePos[0]       ] = 1;
    l_faOr[ (l_vePos[0]+3)%4 ] = 0;
  }

  if(    (    l_elVe[ (l_vePos[1]+3)%4 ] > l_elVe[ (l_vePos[1]+1)%4 ]
           && l_faOr[ l_vePos[1] ] == 3 )
      || l_faOr[ (l_vePos[1]+3)%4 ] != 3 ) {
    l_faOr[ l_vePos[1]       ] = 2;
  }
  else {
    l_faOr[ (l_vePos[1]+3)%4 ] = 2;
  }

  // perform the reordering
  t_idx l_elVeRe[4];
  t_idx l_elFaRe[4];
  t_idx l_elFaElRe[4];
  for( unsigned short l_en = 0; l_en < 4; l_en++ ) {
    l_elVeRe[l_en]   = io_elVe[ l_elVe[l_en] ];
    l_elFaRe[l_en]   = io_elFa[ l_faOr[l_en] ];
    l_elFaElRe[l_en] = io_elFaEl[ l_faOr[l_en] ];
  }
  for( unsigned short l_en = 0; l_en < 4; l_en++ ) {
    io_elVe[l_en]   = l_elVeRe[l_en];
    io_elFa[l_en]   = l_elFaRe[l_en];
    io_elFaEl[l_en] = l_elFaElRe[l_en];
  }
}
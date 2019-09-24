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
#include "GeomTria3.h"
#include "GeomTet4.h"

#include <limits>
#include "io/logging.h"
#include <Eigen/Dense>

double edge_v::mesh::Geom::volume( t_entityType         i_enTy,
                                   double       const (*i_veCrds)[3] ) {
  double l_vol = std::numeric_limits< double >::max();

  if( i_enTy == LINE ) {
    l_vol = GeomLine::length( i_veCrds );
  }
  else if( i_enTy == TRIA3 ) {
    l_vol = GeomTria3::area( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_vol = GeomTet4::volume( i_veCrds );
  }
  else EDGE_V_LOG_FATAL;

  return l_vol;
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

void edge_v::mesh::Geom::tangents( t_entityType         i_enTy,
                                   double       const (*i_veCrds)[3],
                                   double       const   i_nPt[3],
                                   double               o_tangents[2][3] ) {
  // init to zero
  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    o_tangents[0][l_di] = o_tangents[1][l_di] = 0;

  if( i_enTy == LINE ) {
    GeomLine::tangent( i_veCrds,
                       i_nPt,
                       o_tangents[0] );
  }
  else EDGE_V_LOG_FATAL;
}

double edge_v::mesh::Geom::inDiameter( t_entityType         i_enTy,
                                       double       const (*i_veCrds)[3] ) {
  double l_dia = std::numeric_limits< double >::max();
  if( i_enTy == TRIA3 ) {
    l_dia = GeomTria3::inDiameter( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_dia = GeomTet4::inDiameter( i_veCrds );
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
    GeomTria3::normVesFas( i_veCrds,
                           io_elVe,
                           io_elFa,
                           io_elFaEl );
  }
  else if( i_elTy == TET4 ) {
    GeomTet4::normVesFas( i_veCrds,
                          io_elVe,
                          io_elFa,
                          io_elFaEl );
  }
  else EDGE_V_LOG_FATAL;
}
/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
#include "Generic.h"
#include "Line.h"
#include "Quad4r.h"
#include "Tria3.h"
#include "Tet4.h"

#include <limits>
#include "io/logging.h"
#include <Eigen/Dense>

double edge_v::geom::Geom::volume( t_entityType         i_enTy,
                                   double       const (*i_veCrds)[3] ) {
  double l_vol = std::numeric_limits< double >::max();

  if( i_enTy == LINE ) {
    l_vol = Line::length( i_veCrds );
  }
  else if( i_enTy == QUAD4R ) {
    l_vol = Quad4r::area( i_veCrds );
  }
  else if( i_enTy == TRIA3 ) {
    l_vol = Tria3::area( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_vol = Tet4::volume( i_veCrds );
  }
  else EDGE_V_LOG_FATAL;

  return l_vol;
}

void edge_v::geom::Geom::normal( t_entityType         i_enTy,
                                 double       const (*i_veCrds)[3],
                                 double       const   i_nPt[3],
                                 double               o_normal[3] ) {
  if( i_enTy == LINE ) {
    Line::normal( i_veCrds,
                  i_nPt,
                  o_normal );
  }
  else if( i_enTy == TRIA3 ) {
    Tria3::normal( i_veCrds,
                   i_nPt,
                   o_normal );
  }
  else EDGE_V_LOG_FATAL;
}

void edge_v::geom::Geom::tangents( t_entityType         i_enTy,
                                   double       const (*i_veCrds)[3],
                                   double       const   i_nPt[3],
                                   double               o_tangents[2][3] ) {
  // init to zero
  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    o_tangents[0][l_di] = o_tangents[1][l_di] = 0;

  if( i_enTy == LINE ) {
    Line::tangent( i_veCrds,
                   i_nPt,
                   o_tangents[0] );
  }
  else if( i_enTy == TRIA3 ) {
    Tria3::tangents( i_veCrds,
                     i_nPt,
                     o_tangents );
  }
  else EDGE_V_LOG_FATAL;
}

double edge_v::geom::Geom::inDiameter( t_entityType         i_enTy,
                                       double       const (*i_veCrds)[3] ) {
  double l_dia = std::numeric_limits< double >::max();
  if( i_enTy == QUAD4R ) {
    l_dia = Quad4r::inDiameter( i_veCrds );
  }
  else if( i_enTy == TRIA3 ) {
    l_dia = Tria3::inDiameter( i_veCrds );
  }
  else if( i_enTy == TET4 ) {
    l_dia = Tet4::inDiameter( i_veCrds );
  }
  else EDGE_V_LOG_FATAL;

  return l_dia;
}

void edge_v::geom::Geom::normVesFas( t_entityType         i_elTy,
                                     double      const (* i_veCrds)[3],
                                     t_idx              * io_elVe,
                                     t_idx              * io_elFa,
                                     t_idx              * io_elFaEl ) {
  if( i_elTy == QUAD4R ) {
    Quad4r::normVesFas( i_veCrds,
                        io_elVe,
                        io_elFa,
                        io_elFaEl );
  }
  else if( i_elTy == TRIA3 ) {
    Tria3::normVesFas( i_veCrds,
                       io_elVe,
                       io_elFa,
                       io_elFaEl );
  }
  else if( i_elTy == TET4 ) {
    Tet4::normVesFas( i_veCrds,
                      io_elVe,
                      io_elFa,
                      io_elFaEl );
  }
  else EDGE_V_LOG_FATAL;
}

void edge_v::geom::Geom::getVeIdsAd( t_entityType           i_elTy,
                                     t_idx                  i_nFas,
                                     t_idx                  i_elOff,
                                     t_idx          const * i_el,
                                     unsigned short const * i_fa,
                                     t_idx          const * i_elVe,
                                     t_idx          const * i_elFaEl,
                                     unsigned short       * o_veIdsAd ) {
  if( i_elTy == TET4 ) {
    Tet4::getVeIdsAd( i_nFas,
                      i_elOff,
                      i_el,
                      i_fa,
                      i_elVe,
                      i_elFaEl,
                      o_veIdsAd );
  }
  else {
    for( t_idx l_id = 0; l_id < i_nFas; l_id++ ) {
      o_veIdsAd[l_id] = 0;
    }
  }
}


void edge_v::geom::Geom::getFaIdsAd( t_entityType           i_elTy,
                                     t_idx                  i_nFas,
                                     t_idx                  i_elOff,
                                     t_idx          const * i_el,
                                     unsigned short const * i_fa,
                                     t_idx          const * i_elFaEl,
                                     unsigned short       * o_faIdsAd ) {
  Generic::getFaIdsAd( i_elTy,
                       i_nFas,
                       i_elOff,
                       i_el,
                       i_fa,
                       i_elFaEl,
                       o_faIdsAd );
}
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
 * Generic geometry computations.
 **/
#include "Generic.h"
#include "../io/logging.h"
#include <limits>

void edge_v::geom::Generic::sortLex( t_idx   i_n0,
                                     t_idx   i_n1,
                                     t_idx * io_data,
                                     t_idx * o_sorted ) {
  // duplicate data
  t_idx * l_data = new t_idx[ i_n0 * i_n1 ];
  for( t_idx l_id = 0; l_id < i_n0*i_n1; l_id++ ) {
    l_data[l_id] = io_data[l_id];
  }

  // alloc and init sorted aray
  t_idx * l_sorted = new t_idx[ i_n0 ];
  for( t_idx l_i0 = 0; l_i0 < i_n0; l_i0++ ) {
    l_sorted[l_i0] = l_i0;
  }

  // lambda which compares two entries
  auto l_entryLess = [ i_n1, l_data ]( t_idx i_en0,
                                       t_idx i_en1 ) {
    for( t_idx l_i1 = 0; l_i1 < i_n1; l_i1++ ) {
      if( l_data[i_en0*i_n1 + l_i1] < l_data[i_en1*i_n1 + l_i1] ) {
        return true;
      }
      else if( l_data[i_en0*i_n1 + l_i1] > l_data[i_en1*i_n1 + l_i1] ) {
        return false;
      }
    }
    return false;
  };

  // sort the entries
  std::sort( l_sorted,
             l_sorted+i_n0,
             l_entryLess );

  // store the result
  for( t_idx l_i0 = 0; l_i0 < i_n0; l_i0++ ) {
    t_idx l_id = l_sorted[l_i0];
    if( o_sorted != nullptr ) o_sorted[l_i0] = l_id;

    for( t_idx l_i1 = 0; l_i1 < i_n1; l_i1++ ) {
      io_data[l_i0*i_n1 + l_i1] = l_data[l_id*i_n1 + l_i1];
    }
  }

  delete[] l_sorted;
  delete[] l_data;
}

void edge_v::geom::Generic::getFaIdsAd( t_entityType           i_elTy,
                                        t_idx                  i_nFas,
                                        t_idx                  i_elOff,
                                        t_idx    const       * i_el,
                                        unsigned short const * i_fa,
                                        t_idx    const       * i_elFaEl,
                                        unsigned short       * o_faIdsAd ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );

  for( t_idx l_id = 0; l_id < i_nFas; l_id++ ) {
    t_idx l_el = i_el[l_id] + i_elOff;
    t_idx l_fa = i_fa[l_id];

    t_idx l_elAd = i_elFaEl[l_el*l_nElFas + l_fa];
    EDGE_V_CHECK_NE( l_elAd, std::numeric_limits< t_idx >::max() );

    for( unsigned short l_ad = 0; l_ad < l_nElFas; l_ad++ ) {
      if( i_elFaEl[l_elAd*l_nElFas + l_ad] == l_el ) {
        o_faIdsAd[l_id] = l_ad;
        break;
      }
      EDGE_V_CHECK_NE( l_ad+1, l_nElFas );
    }
  }
}
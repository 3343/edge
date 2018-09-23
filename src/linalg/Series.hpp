/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Operations on series.
 **/

#ifndef EDGE_LINALG_SERIES_HPP
#define EDGE_LINALG_SERIES_HPP

#include "constants.hpp"

namespace edge {
  namespace linalg {
      template< unsigned short TL_N_SERIES >
    class Series;
  }
}

/**
 * Operations on series data, typically time series.
 * All operations might operate on multiple series in parallel.
 *
 * @paramt TL_N_TS number of series considered. Only the values are assumed to be different, not the meta-data.
 **/
template< unsigned short TL_N_SERIES >
class edge::linalg::Series {
  // private:
  public:
    /**
     * Integrates the given function in [x1, x2].
     * Remark: Linear interpolation between points is used.
     *
     * @param i_dx distance of the sampled points.
     * @param i_xStart temporal of first sampled point.
     * @param i_nPts number of sampled points.
     * @param i_vals values at the points. These might be cover multiple series.
     * @param i_x1 lower integration point.
     * @param i_x2 upper integration point.
     * @param o_int will be set to the integrated values for the series.
     * @param i_na assumed value for everything not covered by the series.
     **/
    template< typename TL_T_REAL >
    static void integrate( TL_T_REAL         i_dx,
                           TL_T_REAL         i_xStart,
                           std::size_t       i_nPts,
                           TL_T_REAL const (*i_vals)[TL_N_SERIES],
                           TL_T_REAL         i_x1,
                           TL_T_REAL         i_x2,
                           TL_T_REAL         o_int[TL_N_SERIES],
                           TL_T_REAL         i_na=0 ) {
      TL_T_REAL l_dxH = i_dx * (TL_T_REAL) 0.5;

      // reset integration to zero
      for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) o_int[l_se] = 0;

      // derive points covered by [x1, x2]
      std::size_t l_first = std::max( (i_x1 - i_xStart + i_dx)/i_dx, (TL_T_REAL) 0 );
      std::size_t l_last  = std::max( (i_x2 - i_xStart       )/i_dx, (TL_T_REAL) 0 );
                  l_last  = std::min( l_last,  i_nPts-1 );
                  l_first = std::min( l_first, l_last );


      // add contribution of the covered intervals
      for( std::size_t l_pt = l_first+1; l_pt < l_last+1; l_pt++ ) {
        for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) {
          // area of triangle
          o_int[l_se] += l_dxH*( i_vals[l_pt][l_se] - i_vals[l_pt-1][l_se] );
          // constant contribution
          o_int[l_se] += i_vals[l_pt-1][l_se]*i_dx;
        }
      }

      // determine the case we are in
      // false: x1/x2 is outside the defined series
      // true:  x1/x2 it outside the defined series
      bool l_bL = false; bool l_bR = false;
      TL_T_REAL l_tEnd = i_xStart + (i_nPts-1)*i_dx;
      if( i_x1 + TOL.LINALG > i_xStart && i_x1 - TOL.LINALG < l_tEnd ) l_bL = true;
      if( i_x2 + TOL.LINALG > i_xStart && i_x2 - TOL.LINALG < l_tEnd ) l_bR = true;

      // add remaining contribution before
      TL_T_REAL l_diff = i_xStart+l_first*i_dx - i_x1;

      if( l_first > 0 && l_bL ) {
        for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) {
          // area of triangle
          o_int[l_se] +=   (TL_T_REAL) 0.5 * l_diff * (l_diff / i_dx)
                         * (i_vals[l_first-1][l_se] - i_vals[l_first][l_se]);
          // add const contribution
          o_int[l_se] += i_vals[l_first][l_se]*l_diff;
        }
      }
      else {
        for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) {
          o_int[l_se] += i_na * l_diff;
        }
      }


      // add remaining contribution after
      l_diff = i_x2 - (i_xStart+l_last*i_dx);

      if( l_last < i_nPts-1 && l_bR ) {
        for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) {
          // area of triangle
          o_int[l_se] +=   (TL_T_REAL) 0.5 * l_diff
                         * (l_diff / i_dx) * (i_vals[l_last+1][l_se] - i_vals[l_last][l_se]);
          // constant contribution
          o_int[l_se] += i_vals[l_last][l_se] * l_diff;
        }
      }
      else {
        for( unsigned short l_se = 0; l_se < TL_N_SERIES; l_se++ ) {
          o_int[l_se] += i_na * l_diff;
        }
      }
    }
};

#endif

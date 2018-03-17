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
 * Detection criteria for the admissibility of candidate solutions.
 **/
#ifndef EDGE_SC_DETECTIONS_HPP
#define EDGE_SC_DETECTIONS_HPP

#include "constants.hpp"

namespace edge {
  namespace sc {
    template< t_entityType   TL_T_EL,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    class Detections;
  }
}

/**
 * Detection criteria for the admissibility of candidate solutions.
 *
 * @paramt TL_N_QTS number of quantities.
 * @paramt TL_N_CRS number of fused runs.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_N_QTS,
          unsigned short TL_N_CRS >
class edge::sc::Detections {
  private:
    //! epsilon determing the maximum allowed delta of new extrema in the candidate solution w.r.t. to the current solution
    constexpr static const double m_eps  = 0.001;

    //! minimum allowed delta of new extrema in the candidate solution
    constexpr static const double m_eps0 = 0.0001;

  //! number of faces per element
  static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

  public:
    /**
     * Numerical admissibility through the discrete maximum principle in the Voronoi neighborhood of an element.
     *
     * @param i_extCur extrema of the element's current solution, [0][][] is min, [1][][] is max.
     * @param i_extCurAd extrema of the adjacent elements' current solutions. Offset determined by i_elVeEl, [0][][] is min, [1][][] is max.
     * @param i_extCan extrema of the candidate solution, [0][][] is min, [1][][] is max.
     * @param i_nElVeEl number of elements in the Voronoi neighborhood.
     * @param i_elVeEl elements in the Voronoi neighborhood (adjacent through vertices) of the element for which admissiblity is determined.
     * @param o_adm will be set to true if the element is admissible for the respective fused run, false otherwise.
     *
     * @paramt TL_T_REAL floating point precision.
     * @paramt TL_T_INT_LID integer type of the elment ids.
     * 
     **/
    template < typename TL_T_REAL, typename TL_T_INT_LID >
    static void dmpVe( TL_T_REAL      const   i_extCur[2][TL_N_QTS][TL_N_CRS],
                       TL_T_REAL      const (*i_extCurAd)[2][TL_N_QTS][TL_N_CRS],
                       TL_T_REAL      const   i_extCan[2][TL_N_QTS][TL_N_CRS],
                       unsigned short         i_nElVeEl,
                       TL_T_INT_LID   const (*i_elVeEl),
                       bool                   o_adm[TL_N_CRS] ) {
      // init admissibility
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_adm[l_cr] = true;

      // reduced extrema of current solution
      TL_T_REAL l_extCurRed[2][TL_N_QTS][TL_N_CRS];
      // init with element itself
      for( unsigned short l_ex = 0; l_ex < 2; l_ex++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_extCurRed[l_ex][l_qt][l_cr] = i_extCur[l_ex][l_qt][l_cr];
          }
        }
      }
      // reduce over Voronoi neighborhood
      for( unsigned short l_ne = 0; l_ne < i_nElVeEl; l_ne++ ) {
        TL_T_INT_LID l_id = i_elVeEl[l_ne];
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_extCurRed[0][l_qt][l_cr] = std::min( l_extCurRed[0][l_qt][l_cr], i_extCurAd[l_id][0][l_qt][l_cr] );
            l_extCurRed[1][l_qt][l_cr] = std::max( l_extCurRed[1][l_qt][l_cr], i_extCurAd[l_id][1][l_qt][l_cr] );
          }
        }
      }

      // determine parameter delta
      TL_T_REAL l_delta[TL_N_QTS][TL_N_CRS];
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          l_delta[l_qt][l_cr] = (TL_T_REAL) m_eps * ( l_extCurRed[1][l_qt][l_cr] - l_extCurRed[0][l_qt][l_cr] );
          l_delta[l_qt][l_cr] = std::max( (TL_T_REAL) m_eps0, l_delta[l_qt][l_cr] );
        }
      }

      // check discrete maximum principle
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          if( l_extCurRed[0][l_qt][l_cr] - l_delta[l_qt][l_cr] > i_extCan[0][l_qt][l_cr] ||
              i_extCan[1][l_qt][l_cr]                          > l_extCurRed[1][l_qt][l_cr] + l_delta[l_qt][l_cr] ) {
            o_adm[l_cr] = false;
          }
        }
      }

    }

    /**
     * Numerical admissibility through the discrete maximum principle in the face neighborhood of an element.
     *
     * @param i_extCur extrema of the element's current solution, [0][][] is min, [1][][] is max.
     * @param i_extCurAd extrema of the adjacent elements' current solutions.
     * @param i_extCan extrema of the candidate solution, [0][][] is min, [1][][] is max.
     * @param i_elFaEl elements adjacent to the considered element (faces as bridge).
     * @param o_adm will be set to true if the element is admissible for the respective fused run, false otherwise.
     *
     * @paramt TL_T_REAL floating point precision.
     * @paramt TL_T_INT_LID integer type of the elment ids.
     **/
    template < typename TL_T_REAL, typename TL_T_INT_LID >
    static void dmpFa( TL_T_REAL      const   i_extCur[2][TL_N_QTS][TL_N_CRS],
                       TL_T_REAL      const (*i_extCurAd)[2][TL_N_QTS][TL_N_CRS],
                       TL_T_REAL      const   i_extCan[2][TL_N_QTS][TL_N_CRS],
                       TL_T_INT_LID   const   i_elFaEl[TL_N_FAS],
                       bool                   o_adm[TL_N_CRS] ) {
      // init admissibility
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_adm[l_cr] = true;

      // reduced extrema of current solution
      TL_T_REAL l_extCurRed[2][TL_N_QTS][TL_N_CRS];
      // init with element itself
      for( unsigned short l_ex = 0; l_ex < 2; l_ex++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_extCurRed[l_ex][l_qt][l_cr] = i_extCur[l_ex][l_qt][l_cr];
          }
        }
      }
      // reduce over face neighborhood
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        TL_T_INT_LID l_elAd = i_elFaEl[l_fa];
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_extCurRed[0][l_qt][l_cr] = std::min( l_extCurRed[0][l_qt][l_cr], i_extCurAd[l_elAd][0][l_qt][l_cr] );
            l_extCurRed[1][l_qt][l_cr] = std::max( l_extCurRed[1][l_qt][l_cr], i_extCurAd[l_elAd][1][l_qt][l_cr] );
          }
        }
      }

      // determine parameter delta
      TL_T_REAL l_delta[TL_N_QTS][TL_N_CRS];
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          l_delta[l_qt][l_cr] = (TL_T_REAL) m_eps * ( l_extCurRed[1][l_qt][l_cr] - l_extCurRed[0][l_qt][l_cr] );
          l_delta[l_qt][l_cr] = std::max( (TL_T_REAL) m_eps0, l_delta[l_qt][l_cr] );
        }
      }

      // check discrete maximum principle
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          if( l_extCurRed[0][l_qt][l_cr] - l_delta[l_qt][l_cr] > i_extCan[0][l_qt][l_cr] ||
              i_extCan[1][l_qt][l_cr]                          > l_extCurRed[1][l_qt][l_cr] + l_delta[l_qt][l_cr] ) {
            o_adm[l_cr] = false;
          }
        }
      }
    }
};
#endif

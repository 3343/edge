/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Finite Volume solver for the shallow water equations.
 **/
#ifndef EDGE_SWE_FINITE_VOLUME_HPP
#define EDGE_SWE_FINITE_VOLUME_HPP

#include <cmath>
#include "constants.hpp"
#include "Fwave.hpp"

namespace edge {
  namespace swe {
    namespace solvers {
      template< t_entityType   TL_T_EL,
                unsigned short TL_N_CRS >
      class FiniteVolume;
    }
  }
}

/**
 * Finite volume solver for the shallow water equations.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_N_CRS >
class edge::swe::solvers::FiniteVolume {
  private:
    //! number of dimensions
    static const unsigned short TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces per element
    static const unsigned short TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of shallow water quantities
    static const unsigned short TL_N_QTS = (TL_N_DIS == 1) ? 2 : 3;

    /**
     * Computes the CFL time step for the given element.
     *
     * @param i_h water height.
     * @param i_hu momentum.
     * @param i_volume volume of the element.
     * @param i_g gravity.
     * @param i_cfl cfl number.
     **/
    static double computeCflTimeStep( double i_h,
                                      double i_hu,
                                      double i_volume, 
                                      double i_g   = 9.80665,
                                      double i_cfl = 0.4 ) {
      // only elements with water are updated
      if( i_h > 0 ) {
        // compute particle velocity
        double l_u = i_hu / i_h;

        // compute maximum, absolute wave speed
        double l_s = std::abs( l_u ) + std::sqrt( i_g * i_h );

        // compute time step
        double l_dT  = ( i_volume / l_s );
               l_dT *= i_cfl;
        return l_dT;
      }
      else {
        return std::numeric_limits< double >::max();
      }
    }

  public:
    /**
     * Gets the time step statistics according to the CFL-criterion for the entire mesh across all concurrent runs.
     * The computation uses the minimum time step among all concurrent runs in an element.
     * 
     * @param i_nEls number of elements.
     * @param i_charsEl element characteristics.
     * @param i_dofs degrees of freedom for the shallow water equations.
     * @param o_minDt set to minimum occurring time step. This is the minimum among the elements, the minimum of the runs is computed before.
     * @param o_aveDt set to average occuring time step. This is the average among the elements, the average of the runs is computed before.
     * @param o_maxDt set to maximum occuring time step. This is the maximum among the elements, the maximum of the runs is computer before.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_CHARS_EL element characteristics, offereing member variable .volume.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_CHARS_EL >
    static void getTimeStepStatistics( TL_T_LID               i_nEls,
                                       TL_T_CHARS_EL const  * i_charsEl,
                                       TL_T_REAL     const (* i_dofs)[TL_N_QTS][1][TL_N_CRS],
                                       double               & o_minDt,
                                       double               & o_aveDt,
                                       double               & o_maxDt ) {
      // intialize statistics
      o_minDt = std::numeric_limits< double >::max();
      o_aveDt = 0;
      o_maxDt = 0;

      for( TL_T_LID l_el = 0; l_el < i_nEls; l_el++ ) {
        // compute minum CFL associated time step for all concurrent runs
        double l_cDt = std::numeric_limits< double >::max();
        for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
          l_cDt = std::min( l_cDt, computeCflTimeStep( i_dofs[l_el][0][0][l_ru],
                                                       i_dofs[l_el][1][0][l_ru],
                                                       i_charsEl[l_el].volume ) );
        }

        // add element to stats
        o_minDt  = std::min( o_minDt, l_cDt );
        o_aveDt += l_cDt;
        o_maxDt  = std::max( o_maxDt, l_cDt );
      }

      // average
      o_aveDt /= i_nEls;
    }

   /**
    * Computes the net-updates for the given faces using the f-wave solver.
    *
    * @param i_first first face.
    * @param i_size number of faces after first.
    * @param i_faEl faces adjacent to elements.
    * @param i_charsFa face characteristics.
    * @param i_dofs degrees of freedom (height, momentum).
    * @param i_bath bathymetry for the elements.
    * @param o_netUpdate will be set to the net-updates for the faces' adjacent elements.
    *
    * @paramt TL_T_LID integral type of local ids.
    * @paramt TL_T_REAL floating point precision.
    * @paramt TL_T_CHARS_FA struct of the face characterstics, offering member .outNormal.
    **/
   template< typename TL_T_LID,
             typename TL_T_REAL,
             typename TL_T_CHARS_FA >
   static void netUpdates( TL_T_LID               i_first,
                           TL_T_LID               i_size,
                           TL_T_LID      const (* i_faEl)[2],
                           TL_T_CHARS_FA const  * i_charsFa,
                           TL_T_REAL     const (* i_dofs)[TL_N_QTS][1][TL_N_CRS],
                           TL_T_REAL     const (* i_bath)[1][1],
                           TL_T_REAL           (* o_netUpdates)[4][1][TL_N_CRS] ) {
      // compute net-updates
      for( TL_T_LID l_fa = i_first; l_fa < i_first+i_size; l_fa++ ) {
        // determine neighbors
        TL_T_LID l_le, l_ri;

        if( i_charsFa[l_fa].outNormal[0] > 0 ) {
          l_le = i_faEl[l_fa][0];
          l_ri = i_faEl[l_fa][1];
        }
        else {
          l_le = i_faEl[l_fa][1];
          l_ri = i_faEl[l_fa][0];
        }

        // compute net-updates
        solvers::Fwave< TL_N_CRS >::computeNetUpdates( i_dofs[l_le][0][0], i_dofs[l_ri][0][0],
                                                       i_dofs[l_le][1][0], i_dofs[l_ri][1][0],
                                                       i_bath[l_le][0][0], i_bath[l_ri][0][0],
                                                       o_netUpdates[l_fa][0] );
      }
   }

    /**
     * Updates the elements with net-update contribution within a time step.
     *
     * @param i_first first face.
     * @param i_size number of faces after first.
     * @param i_dT time step.
     * @param i_elFa ids of faces adjacent to the elements.
     * @param i_elChars element characteristics.
     * @param i_netUpdates of face-local net-updates, private for conurrent runs.
     * @param io_dofs DOFs: shallow water quantities in the elements, private for concurrent runs.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point precision.
     * @paramt TL_T_CHARS_EL element characteristics, offering member variable .volume.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_CHARS_EL >
    static void update( TL_T_LID               i_first,
                        TL_T_LID               i_size,
                        TL_T_REAL              i_dT,
                        TL_T_LID      const (* i_elFa)[TL_N_FAS],
                        TL_T_CHARS_EL const  * i_elChars,
                        TL_T_REAL     const (* i_netUpdates)[4][1][TL_N_CRS],
                        TL_T_REAL           (* io_dofs)[TL_N_QTS][1][TL_N_CRS] ) {
      // update the elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
        // determine adjacent faces
        TL_T_LID l_left  = i_elFa[l_el][0];
        TL_T_LID l_right = i_elFa[l_el][1];

        // scale update (dt / dx)
        TL_T_REAL l_sca = i_dT * ( TL_T_REAL(1) / TL_T_REAL(i_elChars[l_el].volume) );

        for( int_cfr l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
          // right-going from the left
          io_dofs[l_el][0][0][l_ru] -= l_sca * i_netUpdates[l_left][2][0][l_ru];
          io_dofs[l_el][1][0][l_ru] -= l_sca * i_netUpdates[l_left][3][0][l_ru];

          // left-going from the right
          io_dofs[l_el][0][0][l_ru] -= l_sca * i_netUpdates[l_right][0][0][l_ru];
          io_dofs[l_el][1][0][l_ru] -= l_sca * i_netUpdates[l_right][1][0][l_ru];
        }
      }
    }
};

#endif

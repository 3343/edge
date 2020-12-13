/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
* Copyright (c) 2020, Friedrich Schiller University Jena
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
     * @param i_inDia in-circle diameter of the element.
     * @param i_g gravity.
     * @param i_cfl cfl number.
     **/
    static double computeCflTimeStep( double i_h,
                                      double i_hu,
                                      double i_inDia,
                                      double i_g   = 9.80665,
                                      double i_cfl = 0.9 ) {
      // only elements with water are updated
      if( i_h > 0 ) {
        // compute particle velocity
        double l_u = i_hu / i_h;

        // compute maximum, absolute wave speed
        double l_s = std::abs( l_u ) + std::sqrt( i_g * i_h );

        // compute time step
        double l_dT  = ( i_inDia / l_s );
               l_dT /= TL_N_DIS;
               l_dT *= i_cfl;
        return l_dT;
      }
      else {
        return std::numeric_limits< double >::max();
      }
    }

    /**
     * @brief Computes the left and right quantities by rotating the momentum and applying boundary conditions.
     *
     * @param i_spType sparse type of the face.
     * @param i_elsAd elements adjacent to the face.
     * @param i_n normal of the face.
     * @param i_dofs degrees of freedom.
     * @param i_bath bathymetry.
     * @param o_dofs will be set to the quantities (h and hu) on the left and rightside of the face. [*][][]: left/right [][*][]: quantity, [][][*]: fused run.
     * @param o_bath will be set to bathymetry on the left and right side of the face.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_SP integral type of the sparse type.
     * @paramt TL_T_REAL floating point precision.
     */
    template< typename TL_T_LID,
              typename TL_T_SP,
              typename TL_T_REAL >
    static void qtsLr( TL_T_SP            i_spType,
                       TL_T_LID           i_elsAd[2],
                       TL_T_REAL          i_n[2],
                       TL_T_REAL const (* i_dofs)[TL_N_QTS][1][TL_N_CRS],
                       TL_T_REAL const (* i_bath)[1][1],
                       TL_T_REAL          o_dofs[2][2][TL_N_CRS],
                       TL_T_REAL          o_bath[2] ) {
      for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
        TL_T_LID l_elAd = i_elsAd[l_sd];

        // default right element: no boundary condition
        if(       l_sd == 0 ||
            (    (i_spType & OUTFLOW)    != OUTFLOW
              && (i_spType & REFLECTING) != REFLECTING ) ) {
          // set bathymetry
          o_bath[l_sd] = i_bath[ l_elAd ][0][0];

          // set DOFs
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            // copy over water heights
            o_dofs[l_sd][0][l_cr] = i_dofs[l_elAd][0][0][l_cr];

            // derive face-normal momentum
            o_dofs[l_sd][1][l_cr] = 0;
            for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
              o_dofs[l_sd][1][l_cr] += TL_T_REAL(i_n[l_di]) * i_dofs[l_elAd][1+l_di][0][l_cr];
            }
          }
        }
      }

      if( (i_spType & OUTFLOW) == OUTFLOW ) {
        // set bathymetry
        o_bath[1] = o_bath[0];

        // set DOFs
        for( unsigned short l_qt = 0; l_qt < 2; l_qt++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            o_dofs[1][l_qt][l_cr] = o_dofs[0][l_qt][l_cr];
      }
      else if( (i_spType & REFLECTING) == REFLECTING ) {
        // set bathymetry
        o_bath[1] = o_bath[0];

        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          o_dofs[1][0][l_cr] =  o_dofs[0][0][l_cr];
          o_dofs[1][1][l_cr] = -o_dofs[0][1][l_cr];
        }
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
                                                       i_charsEl[l_el].inDia ) );
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
     * Computes the normal net-updates for the given faces using the f-wave solver.
     *
     * @param i_first first face.
     * @param i_size number of faces after first.
     * @param i_faEl faces adjacent to elements.
     * @param i_charsFa face characteristics.
     * @param i_dofs degrees of freedom (height, momentum).
     * @param i_bath bathymetry for the elements.
     * @param o_nusN will be set to the normal net-updates for the faces' adjacent elements.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point precision.
     * @paramt TL_T_CHARS_FA struct of the face characterstics, offering member .outNormal.
    **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_CHARS_FA >
    static void nusN( TL_T_LID               i_first,
                      TL_T_LID               i_size,
                      TL_T_LID      const (* i_faEl)[2],
                      TL_T_CHARS_FA const  * i_charsFa,
                      TL_T_REAL     const (* i_dofs)[TL_N_QTS][1][TL_N_CRS],
                      TL_T_REAL     const (* i_bath)[1][1],
                      TL_T_REAL           (* o_nusN)[2][2][TL_N_CRS] ) {
      // compute net-updates
      for( TL_T_LID l_fa = i_first; l_fa < i_first+i_size; l_fa++ ) {
        // bathymetry
        TL_T_REAL l_bath[2];

        // left and right DOFs
        TL_T_REAL l_dofs[2][2][TL_N_CRS];

        // adjacent elements
        const TL_T_LID *l_elsAd = i_faEl[l_fa];

        // normal
        TL_T_REAL l_n[TL_N_DIS];
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
          l_n[l_di] = i_charsFa[l_fa].outNormal[l_di];

        qtsLr( i_charsFa[l_fa].spType,
               l_elsAd,
               l_n,
               i_dofs,
               i_bath,
               l_dofs,
               l_bath );

        // init net-updates
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
          for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              o_nusN[l_fa][l_sd][l_qt][l_cr] = 0;
            }
          }
        }

        // take care of dry elements
        TL_T_REAL l_bathSol[2] = { l_bath[0], l_bath[1] };
        if( l_bath[0] < 0 && l_bath[1] < 0 ) {}
        else if( l_bath[1] < 0 ) { // reflecting boundary if left is dry
          l_bathSol[0] = l_bath[1];
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_dofs[0][0][l_cr] =  l_dofs[1][0][l_cr];
            l_dofs[0][1][l_cr] = -l_dofs[1][1][l_cr];
          }
        }
        else if( l_bath[0] < 0 ) { // reflecting boundary if right is dry
          l_bathSol[1] = l_bath[0];
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_dofs[1][0][l_cr] =  l_dofs[0][0][l_cr];
            l_dofs[1][1][l_cr] = -l_dofs[0][1][l_cr];
          }
        }
        else { // no net-update computation if both are dry
          continue;
        }

        // compute net-updates
        solvers::Fwave<
          TL_N_CRS
        >::nusN( l_dofs[0][0],    l_dofs[1][0],
                 l_dofs[0][1],    l_dofs[1][1],
                 l_bathSol[0],    l_bathSol[1],
                 o_nusN[l_fa][0], o_nusN[l_fa][1] );

        // reset net-updates for dry cells
        if( l_bath[0] >= 0 ) {
          for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              o_nusN[l_fa][0][l_qt][l_cr] = 0;
            }
          }
        }
        if( l_bath[1] >= 0 ) {
          for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              o_nusN[l_fa][1][l_qt][l_cr] = 0;
            }
          }
        }
      }
    }

    /**
     * Updates the elements with net-update contribution within a time step.
     *
     * @param i_first first face.
     * @param i_size number of faces after first.
     * @param i_dT time step.
     * @param i_faEl ids of elements adjacent to the faces.
     * @param i_elFa ids of faces adjacent to the elements.
     * @param i_charsFa face characteristics.
     * @param i_charsEl element characteristics.
     * @param i_nusN normal face-local net-updates.
     * @param io_dofs DOFs: shallow water quantities in the elements.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point precision.
     * @paramt TL_T_CHARS_FA face characteristics, offering member variabl .outNormal.
     * @paramt TL_T_CHARS_EL element characteristics, offering member variable .volume.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_EL >
    static void update( TL_T_LID               i_first,
                        TL_T_LID               i_size,
                        TL_T_REAL              i_dT,
                        TL_T_LID      const (* i_faEl)[2],
                        TL_T_LID      const (* i_elFa)[TL_N_FAS],
                        TL_T_CHARS_FA const  * i_charsFa,
                        TL_T_CHARS_EL const  * i_charsEl,
                        TL_T_REAL     const (* i_nusN)[2][2][TL_N_CRS],
                        TL_T_REAL           (* io_dofs)[TL_N_QTS][1][TL_N_CRS] ) {
      // update the elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
        // scale update (dt / dx)
        TL_T_REAL l_sca = i_dT * ( TL_T_REAL(1) / TL_T_REAL(i_charsEl[l_el].volume) );

        for( unsigned short l_ad = 0; l_ad < TL_N_FAS; l_ad++ ) {
          // adjacent face
          const TL_T_LID l_faAd = i_elFa[l_el][l_ad];

          // derive scaling of the net-update
          TL_T_REAL l_scaFa = l_sca * i_charsFa[l_faAd].area;

          // determine offset (left/right distinction) for net-updates
          unsigned short l_sd = ( i_faEl[l_faAd][0] != l_el );

          // normal
          TL_T_REAL l_n[TL_N_DIS];
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
          l_n[l_di] = i_charsFa[l_faAd].outNormal[l_di];

          // apply water height update
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            io_dofs[l_el][0][0][l_cr] -= l_scaFa * i_nusN[l_faAd][l_sd][0][l_cr];

          // apply momentum update
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
            TL_T_REAL l_scaN = l_scaFa * TL_T_REAL(l_n[l_di]);

            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
              io_dofs[l_el][l_di+1][0][l_cr] -= l_scaN * i_nusN[l_faAd][l_sd][1][l_cr];
          }
        }
      }
    }
};

#endif

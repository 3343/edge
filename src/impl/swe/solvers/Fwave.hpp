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
 * F-Wave solver for the one dimensional shallow water equations.
 **/

#ifndef EDGE_SWE_FWAVE_HPP
#define EDGE_SWE_FWAVE_HPP

#include <cmath>

namespace edge {
  namespace swe {
    namespace solvers {
      template< unsigned short TL_N_CRS >
      class Fwave;
    }
  }
}

/**
 * @brief One-dimensional f-Wave solver for the shallow water equations.
 *
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< unsigned short TL_N_CRS >
class edge::swe::solvers::Fwave {
  private:
    /**
     * @brief Derives approximate eigenvalues (Einfeldt speeds), as the extrema of characteristics speeds and Roe speeds.
     *
     * @param i_hL water height for the left element.
     * @param i_hR water height for the right element.
     * @param i_huL momentum for left element.
     * @param i_huR momentum for right element.
     * @param o_lam will be set to eigenvalues.
     * @param i_g gravity constant.
     *
     * @paramt TL_T_REAL floating point precision.
     */
    template< typename TL_T_REAL >
    static void evs( TL_T_REAL const i_hL[TL_N_CRS],
                     TL_T_REAL const i_hR[TL_N_CRS],
                     TL_T_REAL const i_uL[TL_N_CRS],
                     TL_T_REAL const i_uR[TL_N_CRS],
                     TL_T_REAL       o_lam[2][TL_N_CRS],
                     TL_T_REAL       i_g = TL_T_REAL(9.80665) ) {
      // compute characteristic speeds
      TL_T_REAL l_charL[TL_N_CRS], l_charR[TL_N_CRS];
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        l_charL[l_cr] = i_uL[l_cr] - std::sqrt( i_g * i_hL[l_cr] );
        l_charR[l_cr] = i_uR[l_cr] + std::sqrt( i_g * i_hR[l_cr] );
      }

      // compute roe averages
      TL_T_REAL l_hRoe[TL_N_CRS], l_uRoe[TL_N_CRS];
      TL_T_REAL l_roeL[TL_N_CRS], l_roeR[TL_N_CRS];

      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        TL_T_REAL l_hSqrtL = std::sqrt(i_hL[l_cr]);
        TL_T_REAL l_hSqrtR = std::sqrt(i_hR[l_cr]);

        l_hRoe[l_cr]  = TL_T_REAL(0.5) * ( i_hL[l_cr] + i_hL[l_cr] );
        l_uRoe[l_cr]  = l_hSqrtL * i_uL[l_cr] + l_hSqrtR * i_uR[l_cr];
        l_uRoe[l_cr] /= l_hSqrtL + l_hSqrtR;

        TL_T_REAL l_ghSqrtRoe = std::sqrt( i_g * l_hRoe[l_cr]);
        l_roeL[l_cr] = l_uRoe[l_cr] - l_ghSqrtRoe;
        l_roeR[l_cr] = l_uRoe[l_cr] + l_ghSqrtRoe;
      }

      // assign resulting values
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        o_lam[0][l_cr] = std::min( l_charL[l_cr], l_roeL[l_cr] );
        o_lam[1][l_cr] = std::max( l_charR[l_cr], l_roeR[l_cr] );
      }
    }

    /**
     * @brief Computes the wave strength in the wave decomposition.
     *
     * @param i_fJump jump in fluxes.
     * @param i_lam eigenvalues.
     * @param o_beta will be set to wave strengths.
     *
     * @paramt TL_T_REAL floating point precision.
     */
    template< typename TL_T_REAL >
    static void waveStrengths( TL_T_REAL const i_fJump[2][TL_N_CRS],
                               TL_T_REAL const i_lam[2][TL_N_CRS],
                               TL_T_REAL       o_beta[2][TL_N_CRS] ) {
      TL_T_REAL l_adMbc[TL_N_CRS];
      TL_T_REAL l_iR[2][2][TL_N_CRS];

      // compute inverse right eigenvector matrix
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        l_adMbc[l_cr] = TL_T_REAL(1.0) / ( i_lam[1][l_cr] - i_lam[0][l_cr] );

        l_iR[0][0][l_cr] =  l_adMbc[l_cr] * i_lam[1][l_cr];
        l_iR[0][1][l_cr] = -l_adMbc[l_cr];
        l_iR[1][0][l_cr] = -l_adMbc[l_cr] * i_lam[0][l_cr];
        l_iR[1][1][l_cr] =  l_adMbc[l_cr];
      }

      // wave strengths
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        o_beta[0][l_cr]  = l_iR[0][0][l_cr] * i_fJump[0][l_cr];
        o_beta[0][l_cr] += l_iR[0][1][l_cr] * i_fJump[1][l_cr];

        o_beta[1][l_cr]  = l_iR[1][0][l_cr] * i_fJump[0][l_cr];
        o_beta[1][l_cr] += l_iR[1][1][l_cr] * i_fJump[1][l_cr];
      }
    }

    /**
     * @brief Computes 1d net-updates based on the given eigenvector decomposition.
     *
     * @param i_lam eigenvalues.
     * @param i_beta wave strengths.
     * @param o_nusL will be set to left net-updates.
     * @param o_nusR will be set to right net-updates.
     *
     * @paramt TL_T_REAL floating point precision.
     */
    template< typename TL_T_REAL >
    static void nus( TL_T_REAL const i_lam[2][TL_N_CRS],
                     TL_T_REAL const i_beta[2][TL_N_CRS],
                     TL_T_REAL       o_nusL[2][TL_N_CRS],
                     TL_T_REAL       o_nusR[2][TL_N_CRS] ) {
      // init net-updates
      for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          o_nusL[l_qt][l_cr] = 0;
          o_nusR[l_qt][l_cr] = 0;
        }
      }

      // compute net-updates
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        if( i_lam[0][l_cr] < 0 ) {
          o_nusL[0][l_cr] += i_beta[0][l_cr];
          o_nusL[1][l_cr] += i_beta[0][l_cr] * i_lam[0][l_cr];
        }
        else {
          o_nusR[0][l_cr] += i_beta[0][l_cr];
          o_nusR[1][l_cr] += i_beta[0][l_cr] * i_lam[0][l_cr];
        }

        if( i_lam[1][l_cr] > 0 ) {
          o_nusR[0][l_cr] += i_beta[1][l_cr];
          o_nusR[1][l_cr] += i_beta[1][l_cr] * i_lam[1][l_cr];
        }
        else {
          o_nusL[0][l_cr] += i_beta[1][l_cr];
          o_nusL[1][l_cr] += i_beta[1][l_cr] * i_lam[1][l_cr];
        }
      }
    }

  public:
    /**
     * Computes the normal net-updates.
     *
     * @param i_hL water height for the left element.
     * @param i_hR water height for the right element.
     * @param i_huL momentum for left element.
     * @param i_huR momentum for right element.
     * @param i_bL bathymetry for the left element.
     * @param i_bR bathymetry for the right element.
     * @param o_nusL will be set to left net-updates.
     * @param o_nusR will be set to right net-updates.
     * @param i_g gravity constant.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void nusN( TL_T_REAL const i_hL[TL_N_CRS],
                      TL_T_REAL const i_hR[TL_N_CRS],
                      TL_T_REAL const i_huL[TL_N_CRS],
                      TL_T_REAL const i_huR[TL_N_CRS],
                      TL_T_REAL const i_bL,
                      TL_T_REAL const i_bR,
                      TL_T_REAL       o_nusL[2][TL_N_CRS],
                      TL_T_REAL       o_nusR[2][TL_N_CRS],
                      TL_T_REAL       i_g = TL_T_REAL(9.80665) ) {
      // particle veolocities
      TL_T_REAL l_uL[TL_N_CRS];
      TL_T_REAL l_uR[TL_N_CRS];

      // eigenvalues
      TL_T_REAL l_lam[2][TL_N_CRS];

      // jump in fluxes
      TL_T_REAL l_fJump[2][TL_N_CRS];

      // wave strengths
      TL_T_REAL l_beta[2][TL_N_CRS];

      // derive particle velocities
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        l_uL[l_cr] = i_huL[l_cr] / i_hL[l_cr];
        l_uR[l_cr] = i_huR[l_cr] / i_hR[l_cr];
      }

      // compute eigenvalues
      evs( i_hL,   i_hR,
           l_uL,   l_uR,
           l_lam );

      // compute jumps in fluxes
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        // jump in fluxes
        l_fJump[0][l_cr]  = i_huR[l_cr] - i_huL[l_cr];
        l_fJump[1][l_cr]  = i_huR[l_cr] * l_uR[l_cr]  + TL_T_REAL(0.5)  * i_g * i_hR[l_cr] * i_hR[l_cr];
        l_fJump[1][l_cr] -= i_huL[l_cr] * l_uL[l_cr]  + TL_T_REAL(0.5)  * i_g * i_hL[l_cr] * i_hL[l_cr];
        l_fJump[1][l_cr] += TL_T_REAL(0.5) * i_g * ( i_hR[l_cr] + i_hL[l_cr] ) * ( i_bR - i_bL );
      }

      // compute wave strengths
      waveStrengths( l_fJump,
                     l_lam,
                     l_beta );

      // compute net-updates
      nus( l_lam,
           l_beta,
           o_nusL,
           o_nusR );
    }
};

#endif

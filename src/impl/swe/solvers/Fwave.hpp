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
  public:
    /**
     * Updates the DoFs of an element with contribution of the face.
     *
     * @param i_hL water height for the left side element.
     * @param i_hR water height for the right side element.
     * @param i_huL momentum for left side element.
     * @param i_huR momentum for right side element.
     * @param i_bL b for the left side element.
     * @param i_bR b for the right side element.
     * @param o_netUpdates set to the net-update of the face.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void computeNetUpdates( const TL_T_REAL i_hL[TL_N_CRS],
                                   const TL_T_REAL i_hR[TL_N_CRS],
                                   const TL_T_REAL i_huL[TL_N_CRS],
                                   const TL_T_REAL i_huR[TL_N_CRS],
                                   const TL_T_REAL i_bL,
                                   const TL_T_REAL i_bR,
                                         TL_T_REAL o_netUpdates[4][TL_N_CRS] ) {
      // gravity constant
      TL_T_REAL l_g = 9.80665;
      // particle veolocities
      TL_T_REAL l_uL[TL_N_CRS];
      TL_T_REAL l_uR[TL_N_CRS];

      // eigenvalues
      TL_T_REAL l_lambdaL[TL_N_CRS];
      TL_T_REAL l_lambdaR[TL_N_CRS];

      // jump in fluxes
      TL_T_REAL l_fJump[2][TL_N_CRS];

      // scalar of the matrix inversion
      TL_T_REAL l_adMbc[TL_N_CRS];

      // inverse matrix
      TL_T_REAL l_iR[2][2][TL_N_CRS];

      // eigencoefficients
      TL_T_REAL l_beta[2][TL_N_CRS];

      // iterate over concurrent forward runs
      for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
        // particle velocity
        l_uL[l_ru] = i_huL[l_ru] / i_hL[l_ru];
        l_uR[l_ru] = i_huR[l_ru] / i_hR[l_ru];

        // u -/+ sqrt(g*h)
        l_lambdaL[l_ru] = l_uL[l_ru] - std::sqrt( l_g * i_hL[l_ru] );
        l_lambdaR[l_ru] = l_uR[l_ru] + std::sqrt( l_g * i_hR[l_ru] );

        // jump in fluxes
        l_fJump[0][l_ru]  = i_huR[l_ru] - i_huL[l_ru];
        l_fJump[1][l_ru]  = i_huR[l_ru] * l_uR[l_ru]  + TL_T_REAL(0.5)  * l_g * i_hR[l_ru] * i_hR[l_ru];
        l_fJump[1][l_ru] -= i_huL[l_ru] * l_uL[l_ru]  + TL_T_REAL(0.5)  * l_g * i_hL[l_ru] * i_hL[l_ru];
        l_fJump[1][l_ru] += TL_T_REAL(0.5) * l_g * ( i_hR[l_ru] + i_hL[l_ru] ) * ( i_bR - i_bL );

        // compute scalar for 2x2 matrix inverse
        l_adMbc[l_ru] = TL_T_REAL(1.0) / ( l_lambdaR[l_ru] - l_lambdaL[l_ru] );

        l_iR[0][0][l_ru] =  l_adMbc[l_ru] * l_lambdaR[l_ru];
        l_iR[0][1][l_ru] = -l_adMbc[l_ru];
        l_iR[1][0][l_ru] = -l_adMbc[l_ru] * l_lambdaL[l_ru];
        l_iR[1][1][l_ru] =  l_adMbc[l_ru];

        // compute eigen coefficients
        l_beta[0][l_ru]  = l_iR[0][0][l_ru] * l_fJump[0][l_ru];
        l_beta[0][l_ru] += l_iR[0][1][l_ru] * l_fJump[1][l_ru];

        l_beta[1][l_ru]  = l_iR[1][0][l_ru] * l_fJump[0][l_ru];
        l_beta[1][l_ru] += l_iR[1][1][l_ru] * l_fJump[1][l_ru];

        // update DOFs
        // TODO: This fails in supercritical states!
        o_netUpdates[0][l_ru] = l_beta[0][l_ru];
        o_netUpdates[1][l_ru] = l_beta[0][l_ru] * l_lambdaL[l_ru];

        o_netUpdates[2][l_ru] = l_beta[1][l_ru];
        o_netUpdates[3][l_ru] = l_beta[1][l_ru] * l_lambdaR[l_ru];
      }
    }
};

#endif

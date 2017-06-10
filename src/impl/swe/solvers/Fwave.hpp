/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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

#ifndef FWAVE_HPP
#define FWAVE_HPP

#include <algorithm>
#include "constants.hpp"

namespace edge {
  namespace swe {
    namespace solvers {
      class Fwave;
    }
  }
}

class edge::swe::solvers::Fwave {
  public:
    /**
     * Updates the DoFs of an element with contribution of the face.
     *
     * @param i_hL water height for the left/minus side element.
     * @param i_hR water height for the right/plus side element.
     * @param i_huL momentum for left/minus side element.
     * @param i_huR momentum for right/plus side element.
     * @param i_bL b for the left/minus side element.
     * @param i_bR b for the right/plus side element.
     * @param o_netUpdates set to the net-update of the face.
     **/
    static void computeNetUpdates( const t_elementModePrivate1   i_hL[ N_ELEMENT_MODES][N_CRUNS],
                                   const t_elementModePrivate1   i_hR[ N_ELEMENT_MODES][N_CRUNS],
                                   const t_elementModePrivate1   i_huL[N_ELEMENT_MODES][N_CRUNS],
                                   const t_elementModePrivate1   i_huR[N_ELEMENT_MODES][N_CRUNS],
                                   const t_elementModeShared1    i_bL,
                                   const t_elementModeShared1    i_bR,
                                         t_elementModePrivate1   o_netUpdates[4][N_FACE_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      (void) __builtin_assume_aligned(i_hL, ALIGNMENT.CRUNS);  (void) __builtin_assume_aligned(i_hR, ALIGNMENT.CRUNS);
      (void) __builtin_assume_aligned(i_huL, ALIGNMENT.CRUNS); (void) __builtin_assume_aligned(i_huR, ALIGNMENT.CRUNS);
      (void) __builtin_assume_aligned(o_netUpdates, ALIGNMENT.CRUNS);
#endif

      // particle veolocities
      t_elementModePrivate1 l_uL[N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));
      t_elementModePrivate1 l_uR[N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // eigenvalues
      t_elementModePrivate1 l_lambdaL[N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));
      t_elementModePrivate1 l_lambdaR[N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // jump in fluxes
      t_elementModePrivate1 l_fJump[2][N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // scalar of the matrix inversion
      t_elementModePrivate1 l_adMbc[N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // inverse matrix
      t_elementModePrivate1 l_iR[2][2][N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // eigencoefficients
      t_elementModePrivate1 l_beta[2][N_CRUNS] __attribute__((aligned(ALIGNMENT.BASE.STACK)));

      // iterate over concurrent forward runs
#pragma omp simd
      for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
        // particle velocity
        l_uL[l_run] = i_huL[0][l_run] / i_hL[0][l_run];
        l_uR[l_run] = i_huR[0][l_run] / i_hR[0][l_run];

        // u -/+ sqrt(g*h)
        l_lambdaL[l_run] = l_uL[l_run] - std::sqrt( 9.81 * i_hL[0][l_run] );
        l_lambdaR[l_run] = l_uR[l_run] + std::sqrt( 9.81 * i_hR[0][l_run] );

        // jump in fluxes
        l_fJump[0][l_run]  = i_huR[0][l_run] - i_huL[0][l_run];
        l_fJump[1][l_run]  = i_huR[0][l_run] * l_uR[l_run]  + (t_elementModePrivate1) 0.5  * (t_elementModePrivate1) 9.81 * i_hR[0][l_run] * i_hR[0][l_run];
        l_fJump[1][l_run] -= i_huL[0][l_run] * l_uL[l_run]  + (t_elementModePrivate1) 0.5  * (t_elementModePrivate1) 9.81 * i_hL[0][l_run] * i_hL[0][l_run];
        l_fJump[1][l_run] += (t_elementModePrivate1) 0.5    * (t_elementModePrivate1) 9.81 * ( i_hR[0][l_run] + i_hL[0][l_run] * ( i_bR - i_bL ) );

        // compute scalar for 2x2 matrix inverse
        l_adMbc[l_run] = (t_elementModePrivate1) 1.0 / ( l_lambdaR[l_run] - l_lambdaL[l_run] );

        l_iR[0][0][l_run] =  l_adMbc[l_run] * l_lambdaR[l_run];
        l_iR[0][1][l_run] = -l_adMbc[l_run];
        l_iR[1][0][l_run] = -l_adMbc[l_run] * l_lambdaL[l_run];
        l_iR[1][1][l_run] =  l_adMbc[l_run];

        // compute eigen coefficients
        l_beta[0][l_run]  = l_iR[0][0][l_run] * l_fJump[0][l_run];
        l_beta[0][l_run] += l_iR[0][1][l_run] * l_fJump[1][l_run];

        l_beta[1][l_run]  = l_iR[1][0][l_run] * l_fJump[0][l_run];
        l_beta[1][l_run] += l_iR[1][1][l_run] * l_fJump[1][l_run];

        // update DOFs
        // TODO: This fails in supercritical states!
        o_netUpdates[0][0][l_run] = l_beta[0][l_run];
        o_netUpdates[1][0][l_run] = l_beta[0][l_run] * l_lambdaL[l_run];

        o_netUpdates[2][0][l_run] = l_beta[1][l_run];
        o_netUpdates[3][0][l_run] = l_beta[1][l_run] * l_lambdaR[l_run];
      }
    }
};

#endif

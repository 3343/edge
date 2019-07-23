/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2019, Regents of the University of California
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
 * Viscoelastic attenuation.
 **/

#ifndef EDGE_SEISMIC_SETUPS_VISCO_ELASTICITY_HPP
#define EDGE_SEISMIC_SETUPS_VISCO_ELASTICITY_HPP

#include <cmath>
#include "io/logging.h"
#include "../common.hpp"

namespace edge {
  namespace seismic {
    namespace setups {
      class ViscoElasticity;
    }
  }
}

/**
 * Frequency-independent viscoelastic attenuation.
 *
 * Reference: Kaeser, Dumbser, Puente, Igel
 *            An arbitrary high-order Discontinuous Galerkin method for elastic waves on unstructured meshes - III. Viscoelastic attenuation
 *            2007
 * 
 **/
class edge::seismic::setups::ViscoElasticity {
  private:
    /**
     * Computes the anelastic coefficients for the given number of mechanisms and quality factor.
     * 
     * @param i_nMechs number of mechanisms.
     * @param i_qualityFactor (constant) quality factor.
     * @param i_freqs relaxation frequencies.
     * @param o_coeffs will be set to anelastic coefficients.
     **/
    static void anelasticCoeffPs( unsigned short i_nMechs,
                                  double         i_qualityFactor,
                                  double const * i_freqs,
                                  double       * o_coeffs );

    /**
     * Computes the elastic Lame paramters lambda and mu.
     *
     * Reference: Moczo, Bystricky, Kristek, Carcione, Bouchon
     *            Hybrid Modeling of P-SV Seismic Motion at Inhomogeneous Viscoelastic Topograhic Structures.
     *            1997
     *
     * @param i_nMechs number of relaxation mechanisms.
     * @param i_freqRef reference frequency.
     * @param i_freqs relaxation frequencies.
     * @param i_coeffsKappa anelastic kappa coefficients.
     * @param i_coeffsMu anelastic mu coefficients.
     * @param i_lamPhase lame parameter lambda w.r.t. the phase velocity.
     * @param i_muPhase lame parameter mu w.r.t. the phase velocity.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will be set to elastic mu.
     **/
    static void elasticLame( unsigned short         i_nMechs,
                             double                 i_freqRef,
                             double         const * i_freqs,
                             double         const * i_coeffsKappa,
                             double         const * i_coeffsMu,
                             double                 i_lamPhase,
                             double                 i_muPhase,
                             double               & o_lamElastic,
                             double               & o_muElastic );

    /**
     * Computes the anelastic coefficients in terms of the Lame parameters for the given number of mechanisms and quality factor.
     * 
     * @param i_nMechs number of mechanisms.
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param i_qp (constant) quality factor for p-waves.
     * @param i_qs (constant) quality factor for s-waves.
     * @param i_lamPhase lame parameter lambda of the phase velocities.
     * @param i_muPhase lame parameter mu of the phase velocities.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will besetto elastic mu.
     * @param o_coeffs will be set to anelastic coefficients.
     **/
    static void anelasticCoeffLame( unsigned short   i_nMechs,
                                    double           i_freqCen,
                                    double           i_freqRat,
                                    double           i_qp,
                                    double           i_qs,
                                    double           i_lamPhase,
                                    double           i_muPhase,
                                    double         & o_lamElastic,
                                    double         & o_muElastic,
                                    double         * o_coeffs );
  public:
    /**
     * @brief Computes the relaxation frequencies.
     * 
     * @param i_nMechs number of mechanisms (n).
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param o_freqs will be set to relaxation frequencies, 2n-1 in total.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void frequencies( unsigned short   i_nMechs,
                             double           i_freqCen,
                             double           i_freqRat,
                             TL_T_REAL      * o_freqs ) {
      o_freqs[0]  = 2 * M_PI * i_freqCen;
      o_freqs[0] /= std::sqrt( i_freqRat );

      unsigned short l_nFreqs = 2 * i_nMechs - 1;

      for( unsigned short l_fr = 1; l_fr < l_nFreqs; l_fr++ ) {
        o_freqs[l_fr]  = log( o_freqs[0] );
        o_freqs[l_fr] += ( l_fr * log( i_freqRat ) ) / ( 2 * (i_nMechs - 1) );
        o_freqs[l_fr]  = exp( o_freqs[l_fr] );
      }
    }

    /**
     * Computes the dense two-dimensional anelastic source matrix.
     *   The source matrix is split into two parts.
     *   1) The elastic update, for which only the non-zero part,
     *      acting on the elastic stress tensor, is returned.
     *   2) The relaxation frequencies (which are formally multiplied by indentity matrices).
     *
     * @param i_nMechs number of mechanisms.
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param i_qp (constant) quality factor for p-waves.
     * @param i_qs (constant) quality factor for s-waves.
     * @param i_lamPhase lame parameter lambda.
     * @param i_muPhase lame parameter mu.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will be set to elastic mu.
     * @param o_srcElasticStress will be set to source matrices, acting on the stresses (i_nMechs*3*3 entries).
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void src( unsigned short   i_nMechs,
                     double           i_freqCen,
                     double           i_freqRat,
                     double           i_qp,
                     double           i_qs,
                     double           i_lamPhase,
                     double           i_muPhase,
                     double         & o_lamElastic,
                     double         & o_muElastic,
                     TL_T_REAL     (* o_srcElasticStress)[3*3] ) {
      // compute anelastic coefficients
      double *l_coeffs = new double[ i_nMechs*2 ];

      anelasticCoeffLame( i_nMechs,
                          i_freqCen,
                          i_freqRat,
                          i_qp,
                          i_qs,
                          i_lamPhase,
                          i_muPhase,
                          o_lamElastic,
                          o_muElastic,
                          l_coeffs );

      // compute elastic part of the source matrix
      for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
        // set to zero
        for( unsigned short l_ro = 0; l_ro < 3; l_ro++ )
          for( unsigned short l_co = 0; l_co < 3; l_co++ )
            o_srcElasticStress[l_me][l_ro*3 + l_co] = 0;

        // compute nonzeros
        double l_nzs[3] = {0,0,0};
        l_nzs[0]  =   -o_lamElastic * l_coeffs[ l_me            ];
        l_nzs[0] -=  2*o_muElastic  * l_coeffs[ i_nMechs + l_me ];

        l_nzs[1]  = -2*o_muElastic  * l_coeffs[ i_nMechs + l_me ];

        l_nzs[2]  =   -o_lamElastic * l_coeffs[ l_me            ];

        // set the non-zero values
        o_srcElasticStress[l_me][0*3 + 0] = l_nzs[0];
        o_srcElasticStress[l_me][1*3 + 1] = l_nzs[0];

        o_srcElasticStress[l_me][2*3 + 2] = l_nzs[1];

        o_srcElasticStress[l_me][0*3 + 1] = l_nzs[2];
        o_srcElasticStress[l_me][1*3 + 0] = l_nzs[2];
      }

      delete[] l_coeffs;
    }

    /**
     * Computes the two-dimensional anelastic source matrix.
     * Only non-zeros are set (compressed sparse rows).
     *
     * @param i_nMechs number of mechanisms.
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param i_qp (constant) quality factor for p-waves.
     * @param i_qs (constant) quality factor for s-waves.
     * @param i_lamPhase lame parameter lambda.
     * @param i_muPhase lame parameter mu.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will be set to elastic mu.
     * @param o_srcElasticStress will be set to source matrices, acting on the stresses (i_nMechs*5 entries).
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void src( unsigned short   i_nMechs,
                     double           i_freqCen,
                     double           i_freqRat,
                     double           i_qp,
                     double           i_qs,
                     double           i_lamPhase,
                     double           i_muPhase,
                     double         & o_lamElastic,
                     double         & o_muElastic,
                     TL_T_REAL     (* o_srcElasticStress)[5] ) {
      // memory for dense source matrix
      TL_T_REAL (* l_src)[3*3] = (TL_T_REAL (*) [3*3]) new TL_T_REAL[i_nMechs * 3*3];

      // compute dense source matrix
      src( i_nMechs,
           i_freqCen,
           i_freqRat,
           i_qp,
           i_qs,
           i_lamPhase,
           i_muPhase,
           o_lamElastic,
           o_muElastic,
           l_src );

      // extract non-zeros
      for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
        o_srcElasticStress[l_me][ 0] = l_src[l_me][0*3 + 0];
        o_srcElasticStress[l_me][ 1] = l_src[l_me][0*3 + 1];
        o_srcElasticStress[l_me][ 2] = l_src[l_me][1*3 + 0];

        o_srcElasticStress[l_me][ 3] = l_src[l_me][1*3 + 1];
        o_srcElasticStress[l_me][ 4] = l_src[l_me][2*3 + 2];
      }

      // free memory
      delete[] l_src;
    }

    /**
     * Computes the anelastic part of the two-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only the last three columns are set.
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void star( const TL_T_REAL i_jacInv[2][2],
                            TL_T_REAL o_starA[2][3 * 2] ) {
      /*
       * The viscoelastic Jacobians are given as
       *
       *   0    0     0    -1     0 -- 0
       *   0    0     0     0     0 -- 1
       *
       *   0    0     0     0  -1/2 -- 2
       *   |    |     |     |     |
       *   0    1     2     3     4
       *
       *
       *
       *   0    0     0     0     0 -- 0
       *   0    0     0     0    -1 -- 1
       *
       *   0    0     0   -1/2    0 -- 2
       *   |    |     |     |     |
       *   0    1     2     3     4
       */
      // init to zero
      for( unsigned short l_di = 0; l_di < 2; l_di++ )
        for( unsigned short l_m0 = 0; l_m0 < 3; l_m0++ )
          for( unsigned short l_m1 = 0; l_m1 < 2; l_m1++ )
            o_starA[l_di][l_m0*2 + l_m1] = 0;

      // iterate over reference dimension
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        // set non-zero values corresponding to first Jacobian
        o_starA[l_di][0*2 + 3 - 3] += -1.0 * i_jacInv[0][l_di];
        o_starA[l_di][2*2 + 4 - 3] += -0.5 * i_jacInv[0][l_di];

        // set non-zero values corresponding to second Jacobian
        o_starA[l_di][1*2 + 4 - 3] += -1.0 * i_jacInv[1][l_di];
        o_starA[l_di][2*2 + 3 - 3] += -0.5 * i_jacInv[1][l_di];
      }
    }

    /**
     * Computes the anelastic part of the two-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only non-zeros are set (compressed sparse rows).
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void star( const TL_T_REAL i_jacInv[2][2],
                            TL_T_REAL o_starA[2][4] ) {
      // get dense star matrices
      TL_T_REAL l_starA[2][3 * 2];
      star( i_jacInv,
            l_starA );

      // extract potential (depends on jacobians) non-zeros
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        o_starA[l_di][0] = l_starA[l_di][0*2 + 3 - 3];
        o_starA[l_di][1] = l_starA[l_di][1*2 + 4 - 3];
        o_starA[l_di][2] = l_starA[l_di][2*2 + 3 - 3];
        o_starA[l_di][3] = l_starA[l_di][2*2 + 4 - 3];
      }
    }

    /**
     * Computes the two-dimensional anelastic flux solvers for a single face in face-aligned coordinates.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoR density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param o_fsMidL will be set to matrix for left element's contribution.
     * @param o_fsMidR will be set to matrix for right element's contribution.
     **/
    template< typename TL_T_REAL >
    static void fsMid( double    i_rhoL,
                      double    i_rhoR,
                      double    i_lamL,
                      double    i_lamR,
                      double    i_muL, 
                      double    i_muR,
                      TL_T_REAL o_fsMidL[3][5],
                      TL_T_REAL o_fsMidR[3][5] ) {
      // compute wave speeds
      double l_cpL = edge::seismic::common::getVelP( i_rhoL,
                                                    i_lamL,
                                                    i_muL );

      double l_cpR = edge::seismic::common::getVelP( i_rhoR,
                                                    i_lamR,
                                                    i_muR );

      double l_csL = edge::seismic::common::getVelS( i_rhoL,
                                                    i_muL );

      double l_csR = edge::seismic::common::getVelS( i_rhoR,
                                                    i_muR );

      // init matrices
      for( unsigned short l_ro = 0; l_ro < 3; l_ro++ ) {
        for( unsigned short l_co = 0; l_co < 5; l_co++ ) {
          o_fsMidL[l_ro][l_co] = 0;
          o_fsMidR[l_ro][l_co] = 0;
        }
      }

      // set non-zeros of left flux solver
      o_fsMidL[0][0]  = l_cpL;
      o_fsMidL[0][0] /= i_lamL + 2*i_muL + i_rhoR*l_cpL*l_cpR;

      o_fsMidL[0][3]  = -(i_lamL + 2*i_muL);
      o_fsMidL[0][3] /= i_lamL + 2*i_muL + i_rhoR * l_cpL * l_cpR;

      o_fsMidL[2][2]  = l_csL * l_csR;
      o_fsMidL[2][2] /= 2*i_muL*l_csR + 2*i_muR*l_csL;

      o_fsMidL[2][4]  = -i_muL * l_csR;
      o_fsMidL[2][4] /= 2*i_muL*l_csR + 2*i_muR*l_csL;

      // set non-zeros of right flux solver
      o_fsMidR[0][0]  = -l_cpL;
      o_fsMidR[0][0] /= i_lamL + 2*i_muL + i_rhoR*l_cpL*l_cpR;

      o_fsMidR[0][3]  = -l_cpL * ( i_lamR + 2*i_muR );
      o_fsMidR[0][3] /= l_cpR * ( i_lamL + 2*i_muL + i_rhoR * l_cpL * l_cpR );

      o_fsMidR[2][2]  = -l_csL * l_csR;
      o_fsMidR[2][2] /= 2*i_muL*l_csR + 2*i_muR*l_csL;

      o_fsMidR[2][4]  = -i_muR * l_csL;
      o_fsMidR[2][4] /= 2*i_muL * l_csR + 2*i_muR * l_csL;
    }

    /**
     * Computes the dense three-dimensional anelastic source matrix.
     *   The source matrix is split into two parts.
     *   1) The elastic update, for which only the non-zero part,
     *      acting on the elastic stress tensor, is returned.
     *   2) The relaxation frequencies (which are formally multiplied by indentity matrices).
     *
     * @param i_nMechs number of mechanisms.
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param i_qp (constant) quality factor for p-waves.
     * @param i_qs (constant) quality factor for s-waves.
     * @param i_lamPhase lame parameter lambda.
     * @param i_muPhase lame parameter mu.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will be set to elastic mu.
     * @param o_srcElasticStress will be set to source matrices, acting on the stresses (i_nMechs*6*6 entries).
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void src( unsigned short   i_nMechs,
                     double           i_freqCen,
                     double           i_freqRat,
                     double           i_qp,
                     double           i_qs,
                     double           i_lamPhase,
                     double           i_muPhase,
                     double         & o_lamElastic,
                     double         & o_muElastic,
                     TL_T_REAL     (* o_srcElasticStress)[6*6] ) {
      // compute anelastic coefficients
      double *l_coeffs = new double[ i_nMechs*2 ];

      anelasticCoeffLame( i_nMechs,
                          i_freqCen,
                          i_freqRat,
                          i_qp,
                          i_qs,
                          i_lamPhase,
                          i_muPhase,
                          o_lamElastic,
                          o_muElastic,
                          l_coeffs );

      // compute elastic part of the source matrix
      for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
        // set to zero
        for( unsigned short l_ro = 0; l_ro < 6; l_ro++ )
          for( unsigned short l_co = 0; l_co < 6; l_co++ )
            o_srcElasticStress[l_me][l_ro*6 + l_co] = 0;

        // compute nonzeros
        double l_nzs[3] = {0,0,0};
        l_nzs[0]  =   -o_lamElastic * l_coeffs[ l_me            ];
        l_nzs[0] -=  2*o_muElastic  * l_coeffs[ i_nMechs + l_me ];

        l_nzs[1]  = -2*o_muElastic  * l_coeffs[ i_nMechs + l_me ];

        l_nzs[2]  =   -o_lamElastic * l_coeffs[ l_me            ];

        // set the non-zero values
        o_srcElasticStress[l_me][0*6 + 0] = l_nzs[0];
        o_srcElasticStress[l_me][1*6 + 1] = l_nzs[0];
        o_srcElasticStress[l_me][2*6 + 2] = l_nzs[0];

        o_srcElasticStress[l_me][3*6 + 3] = l_nzs[1];
        o_srcElasticStress[l_me][4*6 + 4] = l_nzs[1];
        o_srcElasticStress[l_me][5*6 + 5] = l_nzs[1];

        o_srcElasticStress[l_me][0*6 + 1] = l_nzs[2];
        o_srcElasticStress[l_me][0*6 + 2] = l_nzs[2];
        o_srcElasticStress[l_me][1*6 + 0] = l_nzs[2];
        o_srcElasticStress[l_me][1*6 + 2] = l_nzs[2];
        o_srcElasticStress[l_me][2*6 + 0] = l_nzs[2];
        o_srcElasticStress[l_me][2*6 + 1] = l_nzs[2];
      }

      delete[] l_coeffs;
    }

    /**
     * Computes the three-dimensional anelastic source matrix.
     * Only non-zeros are set (compressed sparse rows).
     * 
     * @param i_nMechs number of mechanisms.
     * @param i_freqCen central frequency.
     * @param i_freqRat frequency ratio.
     * @param i_qp (constant) quality factor for p-waves.
     * @param i_qs (constant) quality factor for s-waves.
     * @param i_lamPhase lame parameter lambda.
     * @param i_muPhase lame parameter mu.
     * @param o_lamElastic will be set to elastic lambda.
     * @param o_muElastic will be set to elastic mu.
     * @param o_srcElasticStress will be set to source matrices, acting on the stresses (i_nMechs*6*6 entries).
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void src( unsigned short   i_nMechs,
                     double           i_freqCen,
                     double           i_freqRat,
                     double           i_qp,
                     double           i_qs,
                     double           i_lamPhase,
                     double           i_muPhase,
                     double         & o_lamElastic,
                     double         & o_muElastic,
                     TL_T_REAL     (* o_srcElasticStress)[12] ) {
      // memory for dense source matrix
      TL_T_REAL (* l_src)[6*6] = (TL_T_REAL (*) [6*6]) new TL_T_REAL[i_nMechs * 6*6];

      // compute dense source matrix
      src( i_nMechs,
           i_freqCen,
           i_freqRat,
           i_qp,
           i_qs,
           i_lamPhase,
           i_muPhase,
           o_lamElastic,
           o_muElastic,
           l_src );

      // extract non-zeros
      for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
        o_srcElasticStress[l_me][ 0] = l_src[l_me][0*6 + 0];
        o_srcElasticStress[l_me][ 1] = l_src[l_me][0*6 + 1];
        o_srcElasticStress[l_me][ 2] = l_src[l_me][0*6 + 2];

        o_srcElasticStress[l_me][ 3] = l_src[l_me][1*6 + 0];
        o_srcElasticStress[l_me][ 4] = l_src[l_me][1*6 + 1];
        o_srcElasticStress[l_me][ 5] = l_src[l_me][1*6 + 2];


        o_srcElasticStress[l_me][ 6] = l_src[l_me][2*6 + 0];
        o_srcElasticStress[l_me][ 7] = l_src[l_me][2*6 + 1];
        o_srcElasticStress[l_me][ 8] = l_src[l_me][2*6 + 2];

        o_srcElasticStress[l_me][ 9] = l_src[l_me][3*6 + 3];
        o_srcElasticStress[l_me][10] = l_src[l_me][4*6 + 4];
        o_srcElasticStress[l_me][11] = l_src[l_me][5*6 + 5];
      }

      // free memory
      delete[] l_src;
    }

    /**
     * Computes the anelastic part of the three-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only the last three columns are set.
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void star( const TL_T_REAL i_jacInv[3][3],
                            TL_T_REAL o_starA[3][6 * 3] ) {
      /*
       * The viscoelastic Jacobians are given as (see Kaeser et. al, 2007):
       *
       *   0    0    0     0    0    0    -1     0    0 -- 0
       *   0    0    0     0    0    0     0     0    0 -- 1
       *   0    0    0     0    0    0     0     0    0 -- 2
       *
       *   0    0    0     0    0    0     0  -1/2    0 -- 3
       *   0    0    0     0    0    0     0     0    0 -- 4
       *   0    0    0     0    0    0     0     0 -1/2 -- 5
       *   |    |    |     |    |    |     |     |    |
       *   0    1    2     3    4    5     6     7    8
       *
       *
       *
       *   0    0    0     0    0    0     0     0    0 -- 0
       *   0    0    0     0    0    0     0    -1    0 -- 1
       *   0    0    0     0    0    0     0     0    0 -- 2
       *
       *   0    0    0     0    0    0   -1/2    0    0 -- 3
       *   0    0    0     0    0    0     0     0 -1/2 -- 4
       *   0    0    0     0    0    0     0     0    0 -- 5
       *   |    |    |     |    |    |     |     |    |
       *   0    1    2     3    4    5     6     7    8
       *
       *
       *
       *   0    0    0     0    0    0     0     0    0 -- 0
       *   0    0    0     0    0    0     0     0    0 -- 1
       *   0    0    0     0    0    0     0     0   -1 -- 2
       *
       *   0    0    0     0    0    0     0     0    0 -- 3
       *   0    0    0     0    0    0     0  -1/2    0 -- 4
       *   0    0    0     0    0    0  -1/2     0    0 -- 5
       *   |    |    |     |    |    |     |     |    |
       *   0    1    2     3    4    5     6     7    8
       */
      // init to zero
      for( unsigned short l_di = 0; l_di < 3; l_di++ )
        for( unsigned short l_m0 = 0; l_m0 < 6; l_m0++ )
          for( unsigned short l_m1 = 0; l_m1 < 3; l_m1++ )
            o_starA[l_di][l_m0*3 + l_m1] = 0;

      // iterate over reference dimension
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        // set non-zero values corresponding to first Jacobian
        o_starA[l_di][0*3 + 6 - 6] += -1.0 * i_jacInv[0][l_di];
        o_starA[l_di][3*3 + 7 - 6] += -0.5 * i_jacInv[0][l_di];
        o_starA[l_di][5*3 + 8 - 6] += -0.5 * i_jacInv[0][l_di];

        // set non-zero values corresponding to second Jacobian
        o_starA[l_di][1*3 + 7 - 6] += -1.0 * i_jacInv[1][l_di];
        o_starA[l_di][3*3 + 6 - 6] += -0.5 * i_jacInv[1][l_di];
        o_starA[l_di][4*3 + 8 - 6] += -0.5 * i_jacInv[1][l_di];

        // set non-zero values corresponding to third Jacobian
        o_starA[l_di][2*3 + 8 - 6] += -1.0 * i_jacInv[2][l_di];
        o_starA[l_di][4*3 + 7 - 6] += -0.5 * i_jacInv[2][l_di];
        o_starA[l_di][5*3 + 6 - 6] += -0.5 * i_jacInv[2][l_di];
    }
  }

  /**
   *  Computes the anelastic part of the three-dimensional star matrices, excluding the scaling with the relaxation frequencies.
   *  Only non-zeros are set (compressed sparse rows).
   *        
   * @param i_jacInv inverse jacobian matrix of the element.
   * @param o_starA will be set to assembled star matrix.
   *
   * @paramt TL_T_REAL floating point precision.
   **/
  template< typename TL_T_REAL >
  static void star( const TL_T_REAL i_jacInv[3][3],
                          TL_T_REAL o_starA[3][9] ) {
    // get dense star matrices
    TL_T_REAL l_starA[3][6 * 3];
    star( i_jacInv,
          l_starA );

    // extract potential (depends on jacobians) non-zeros
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      o_starA[l_di][0] = l_starA[l_di][0*3 + 6 - 6];
      o_starA[l_di][1] = l_starA[l_di][1*3 + 7 - 6];
      o_starA[l_di][2] = l_starA[l_di][2*3 + 8 - 6];
      o_starA[l_di][3] = l_starA[l_di][3*3 + 6 - 6];
      o_starA[l_di][4] = l_starA[l_di][3*3 + 7 - 6];
      o_starA[l_di][5] = l_starA[l_di][4*3 + 7 - 6];
      o_starA[l_di][6] = l_starA[l_di][4*3 + 8 - 6];
      o_starA[l_di][7] = l_starA[l_di][5*3 + 6 - 6];
      o_starA[l_di][8] = l_starA[l_di][5*3 + 8 - 6];
    }
  }

  /**
   * Computes the three-dimensional anelastic flux solvers for a single face in face-aligned coordinates.
   *
   * @param i_rhoL density of the left element.
   * @param i_rhoR density of the right element.
   * @param i_lamL Lame parameter lambda of the left element.
   * @param i_lamR Lame parameter lambda of the right element.
   * @param i_muL Lame parameter mu of the left element.
   * @param i_muR Lame parameter mu of the right element.
   * @param o_fsMidL will be set to matrix for left element's contribution.
   * @param o_fsMidR will be set to matrix for right element's contribution.
   *
   * @paramt TL_T_REAL real type.
   **/
  template< typename TL_T_REAL >
  static void fsMid( double    i_rhoL,
                     double    i_rhoR,
                     double    i_lamL,
                     double    i_lamR,
                     double    i_muL, 
                     double    i_muR,
                     TL_T_REAL o_fsMidL[6][9],
                     TL_T_REAL o_fsMidR[6][9] ) {
    // compute wave speeds
    double l_cpL = edge::seismic::common::getVelP( i_rhoL,
                                                   i_lamL,
                                                   i_muL );

    double l_cpR = edge::seismic::common::getVelP( i_rhoR,
                                                   i_lamR,
                                                   i_muR );

    double l_csL = edge::seismic::common::getVelS( i_rhoL,
                                                   i_muL );

    double l_csR = edge::seismic::common::getVelS( i_rhoR,
                                                   i_muR );

    // init matrices
    for( unsigned short l_ro = 0; l_ro < 6; l_ro++ ) {
      for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
        o_fsMidL[l_ro][l_co] = 0;
        o_fsMidR[l_ro][l_co] = 0;
      }
    }

    // set non-zeros of left flux solver
    o_fsMidL[0][0]  = l_cpL * l_cpR;
    o_fsMidL[0][0] /= i_lamL * l_cpR + i_lamR * l_cpL + 2 * i_muL * l_cpR + 2 * i_muR * l_cpL;

    o_fsMidL[0][6]  = -l_cpR * ( i_lamL + 2 * i_muL);
    o_fsMidL[0][6] /= i_lamL * l_cpR + i_lamR * l_cpL + 2 * i_muL * l_cpR + 2 * i_muR * l_cpL;

    o_fsMidL[3][3]  = l_csL * l_csR;
    o_fsMidL[3][3] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidL[3][7]  = -i_muL * l_csR;
    o_fsMidL[3][7] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidL[5][5]  = l_csL * l_csR;
    o_fsMidL[5][5] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidL[5][8]  = -i_muL * l_csR;
    o_fsMidL[5][8] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    // set non-zeros of right flux solver
    o_fsMidR[0][0]  = -l_cpL * l_cpR;
    o_fsMidR[0][0] /= i_lamL * l_cpR + i_lamR * l_cpL + 2 * i_muL * l_cpR + 2 * i_muR * l_cpL;

    o_fsMidR[0][6]  = -l_cpL * ( i_lamR + 2 * i_muR);
    o_fsMidR[0][6] /= i_lamL * l_cpR + i_lamR * l_cpL + 2 * i_muL * l_cpR + 2 * i_muR * l_cpL;

    o_fsMidR[3][3]  = -l_csL * l_csR;
    o_fsMidR[3][3] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidR[3][7]  = -i_muR * l_csL;
    o_fsMidR[3][7] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidR[5][5]  = -l_csL * l_csR;
    o_fsMidR[5][5] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;

    o_fsMidR[5][8]  = -i_muR * l_csL;
    o_fsMidR[5][8] /= 2 * i_muL * l_csR + 2 * i_muR * l_csL;
  }
};


#endif
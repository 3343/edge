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
#ifndef EDGE_SEISMIC_SETUPS_VISCO_ELASTICITY_H
#define EDGE_SEISMIC_SETUPS_VISCO_ELASTICITY_H

#include <cmath>
#include <limits>
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
    static void fsMid( double i_rhoL,
                       double i_rhoR,
                       double i_lamL,
                       double i_lamR,
                       double i_muL,
                       double i_muR,
                       double o_fsMidL[3][5],
                       double o_fsMidR[3][5] );

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
     **/
    static void fsMid( double i_rhoL,
                       double i_rhoR,
                       double i_lamL,
                       double i_lamR,
                       double i_muL,
                       double i_muR,
                       double o_fsMidL[6][9],
                       double o_fsMidR[6][9] );

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
     **/
    static void star( double const i_jacInv[2][2],
                      double       o_starA[2][3 * 2] );

    /**
     * Computes the anelastic part of the two-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only non-zeros are set (compressed sparse rows).
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     **/
    static void star( double const i_jacInv[2][2],
                      double       o_starA[2][4] );

    /**
     * Computes the anelastic part of the three-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only the last three columns are set.
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     **/
    static void star( double const i_jacInv[3][3],
                      double       o_starA[3][6 * 3] );

    /**
     * Computes the anelastic part of the three-dimensional star matrices, excluding the scaling with the relaxation frequencies.
     * Only non-zeros are set (compressed sparse rows).
     *
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starA will be set to assembled star matrix.
     **/
    static void star( double const i_jacInv[3][3],
                      double       o_starA[3][9] );

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
     * Sets up the two-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_tA transformation of stresses from face-aligned coordinates to Cartesian.
     * @param i_tm1 transformation of elastic quantities from Cartesian coordinates to face-aligned.
     * @param o_fsAl will be set to anelastic flux solver, applied to the left element's DOFs.
     * @param o_fsAr will be set to anelastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     **/
    static void fs( double       i_rhoL,
                    double       i_rhoR,
                    double       i_lamL,
                    double       i_lamR,
                    double       i_muL,
                    double       i_muR,
                    double const i_tA[3][3],
                    double const i_tm1[5][5],
                    double       o_fsAl[3*5],
                    double       o_fsAr[3*5],
                    bool         i_freeSurface );

    /**
     * Sets up the three-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_tA transformation of stresses from face-aligned coordinates to Cartesian.
     * @param i_tm1 transformation of elastic quantities from Cartesian coordinates to face-aligned.
     * @param o_fsAl will be set to anelastic flux solver, applied to the left element's DOFs.
     * @param o_fsAr will be set to anelastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     **/
    static void fs( double       i_rhoL,
                    double       i_rhoR,
                    double       i_lamL,
                    double       i_lamR,
                    double       i_muL,
                    double       i_muR,
                    double const i_tA[6][6],
                    double const i_tm1[9][9],
                    double       o_fsAl[6*9],
                    double       o_fsAr[6*9],
                    bool         i_freeSurface );
};

#endif
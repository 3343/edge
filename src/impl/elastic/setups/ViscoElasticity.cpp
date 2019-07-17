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
#include "ViscoElasticity.h"

// ignore shadowing issues of Eigen library
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wshadow"
#include "submodules/eigen/Eigen/Dense"
#pragma GCC diagnostic pop

void edge::seismic::setups::ViscoElasticity::elasticLame( unsigned short         i_nMechs,
                                                          double                 i_freqRef,
                                                          double         const * i_freqs,
                                                          double         const * i_coeffsKappa,
                                                          double         const * i_coeffsMu,
                                                          double                 i_lamPhase,
                                                          double                 i_muPhase,
                                                          double               & o_lamElastic,
                                                          double               & o_muElastic ) {
  // compute parameters theta1, theta2 and r
  double l_theta1[2] = {1, 1};
  double l_theta2[2] = {0, 0};

  for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
    double l_freqRat = i_freqRef / i_freqs[l_me*2];
    double l_sca = 1.0 / (1 + l_freqRat*l_freqRat);

    l_theta1[0] -= i_coeffsKappa[l_me] * l_sca;
    l_theta1[1] -= i_coeffsMu[l_me]    * l_sca;

    l_theta2[0] += i_coeffsKappa[l_me] * l_freqRat * l_sca;
    l_theta2[1] += i_coeffsMu[l_me]    * l_freqRat * l_sca;
  }

  double l_r[2] = {std::numeric_limits< double >::max(), std::numeric_limits< double >::max()};
  for( unsigned short l_la = 0; l_la < 2; l_la++ ) {
    l_r[l_la] = l_theta1[l_la] * l_theta1[l_la] + l_theta2[l_la] * l_theta2[l_la];
    l_r[l_la] = std::sqrt( l_r[l_la] );
  }

  // compute elastic Lame parameters
  o_muElastic = i_muPhase * ( l_r[1] + l_theta1[1] ) / ( 2 * l_r[1] * l_r[1] );

  double l_kappa = (i_lamPhase + 2 * i_muPhase);
  l_kappa *= ( l_r[0] + l_theta1[0] ) / ( 2 * l_r[0] * l_r[0] );

  o_lamElastic = l_kappa - 2 * o_muElastic;
}

void edge::seismic::setups::ViscoElasticity::anelasticCoeffPs( unsigned short i_nMechs,
                                                               double         i_qualityFactor,
                                                               double const * i_freqs,
                                                               double       * o_coeffs ) {
  // number of frequencies
  unsigned short l_nFreqs = 2 * i_nMechs - 1;

  // inverse quality factor
  double l_qInv =  1.0 / i_qualityFactor;

  // allocate memory for matrix entries
  double *l_matRaw = new double[l_nFreqs * i_nMechs];

  // assemble matrix
  for( unsigned short l_fr = 0; l_fr < l_nFreqs; l_fr++ ) {
    for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
      double l_mf = i_freqs[ l_me * 2 ];

      l_matRaw[ l_fr * i_nMechs + l_me ]  = l_mf * i_freqs[l_fr] + l_mf*l_mf * l_qInv;
      l_matRaw[ l_fr * i_nMechs + l_me ] /= i_freqs[l_fr]*i_freqs[l_fr] + l_mf*l_mf;        }
  }

  // "assemble" right-hand-side (we assume constant quality factors)
  double * l_rhsRaw = new double[ l_nFreqs ];

  for( unsigned short l_fr = 0; l_fr < l_nFreqs; l_fr++ ) {
    l_rhsRaw[l_fr] = l_qInv;
  }

  // convert to Eigen datastructures
  Eigen::Map <
    Eigen::Matrix <
    double,
    Eigen::Dynamic,
    Eigen::Dynamic,
    Eigen::RowMajor > > l_mat( l_matRaw, l_nFreqs, i_nMechs );

  Eigen::Map <
    Eigen::Matrix <
    double,
    Eigen::Dynamic,
    1,
    Eigen::ColMajor > > l_rhs( l_rhsRaw, l_nFreqs );

  // least squares fit
  Eigen::Map <
    Eigen::Matrix <
    double,
    Eigen::Dynamic,
    1,
    Eigen::ColMajor > > l_fit( o_coeffs, i_nMechs );

  // solve problem
  l_fit = l_mat.bdcSvd( Eigen::ComputeThinU | Eigen::ComputeThinV ).solve( l_rhs );

  delete[] l_matRaw;
  delete[] l_rhsRaw;
}

void edge::seismic::setups::ViscoElasticity::anelasticCoeffLame( unsigned short   i_nMechs,
                                                                 double           i_freqCen,
                                                                 double           i_freqRat,
                                                                 double           i_qp,
                                                                 double           i_qs,
                                                                 double           i_lamPhase,
                                                                 double           i_muPhase,
                                                                 double         & o_lamElastic,
                                                                 double         & o_muElastic,
                                                                 double         * o_coeffs ) {
  // get relaxation frequencies
  double *l_freqs = new double[2*i_nMechs - 1];

  frequencies( i_nMechs,
               i_freqCen,
               i_freqRat,
               l_freqs );

  // compute p/s version
  for( unsigned short l_wt = 0; l_wt < 2; l_wt++ )
    anelasticCoeffPs( i_nMechs,
                      (l_wt == 0) ? i_qp : i_qs,
                      l_freqs,
                      o_coeffs+(l_wt*i_nMechs) );

  // compute elastic lame parameters
  double l_freqRef = 2 * M_PI * i_freqCen;

  elasticLame( i_nMechs,
               l_freqRef,
               l_freqs,
               o_coeffs+0,
               o_coeffs+i_nMechs,
               i_lamPhase,
               i_muPhase,
               o_lamElastic,
               o_muElastic );

  // adjust lambda coefficients
  for( unsigned short l_me = 0; l_me < i_nMechs; l_me++ ) {
    double l_sca = 2.0 * o_muElastic / o_lamElastic;
    o_coeffs[l_me]  = ( 1.0 + l_sca ) * o_coeffs[l_me];
    o_coeffs[l_me] -= l_sca * o_coeffs[i_nMechs+l_me];
  }

  delete[] l_freqs;
}
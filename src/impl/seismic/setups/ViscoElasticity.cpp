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
#include "linalg/Matrix.h"

#include "ViscoElasticity.h"

// ignore shadowing issues of Eigen library
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
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

void edge::seismic::setups::ViscoElasticity::star( double const i_jacInv[2][2],
                                                   double       o_starA[2][3 * 2] ) {
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

void edge::seismic::setups::ViscoElasticity::star( double const i_jacInv[2][2],
                                                   double       o_starA[2][4] ) {
  // get dense star matrices
  double l_starA[2][3 * 2];
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

void edge::seismic::setups::ViscoElasticity::star( const double i_jacInv[3][3],
                                                         double o_starA[3][6 * 3] ) {
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

void edge::seismic::setups::ViscoElasticity::star( double const i_jacInv[3][3],
                                                   double       o_starA[3][9] ) {
  // get dense star matrices
  double l_starA[3][6 * 3];
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

void edge::seismic::setups::ViscoElasticity::fsMid( double i_rhoL,
                                                    double i_rhoR,
                                                    double i_lamL,
                                                    double i_lamR,
                                                    double i_muL,
                                                    double i_muR,
                                                    double o_fsMidL[3][5],
                                                    double o_fsMidR[3][5] ) {
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

void edge::seismic::setups::ViscoElasticity::fsMid( double i_rhoL,
                                                    double i_rhoR,
                                                    double i_lamL,
                                                    double i_lamR,
                                                    double i_muL,
                                                    double i_muR,
                                                    double o_fsMidL[6][9],
                                                    double o_fsMidR[6][9] ) {
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

void edge::seismic::setups::ViscoElasticity::fs( double       i_rhoL,
                                                 double       i_rhoR,
                                                 double       i_lamL,
                                                 double       i_lamR,
                                                 double       i_muL,
                                                 double       i_muR,
                                                 double const i_tA[3][3],
                                                 double const i_tm1[5][5],
                                                 double       o_fsAl[3*5],
                                                 double       o_fsAr[3*5],
                                                 bool         i_freeSurface ) {
  // intermediate matrices
  double l_fsMidA[2][3][5];
  double l_tmp[2][3][5];

  // compute mid part of the solver
  seismic::setups::ViscoElasticity::fsMid( i_rhoL,
                                           i_rhoR,
                                           i_lamL,
                                           i_lamR,
                                           i_muL,
                                           i_muR,
                                           l_fsMidA[0],
                                           l_fsMidA[1] );

  // apply back-rotation to face-aligned coordinate system and mid part
  for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
    linalg::Matrix::matMulB0( 3, 5, 3,
                              i_tA[0],
                              l_fsMidA[l_sd][0],
                              l_tmp[l_sd][0] );
  }

  // rotate left solver to face-aligned coordinates (element exists per definition of our boundary conditions)
  linalg::Matrix::matMulB0( 3, 5, 5,
                            l_tmp[0][0],
                            i_tm1[0],
                            o_fsAl );

  /**
   * Free surface boundary conditions mirror (rotated) x-components of the stress tensor.
   * Reference: Eq. (51) in
   *            Martin Kaeser and Michael Dumbser
   *            An arbitrary high-order discontinuous Galerkin method for elastic
   *            waves on unstructured meshes – I. The two-dimensional isotropic
   *            case with external source terms
   **/
  if( i_freeSurface == true ) {
    // multiply tmpR with diag( -1, 1, -1, 1, 1 ) from the right.
    //                           0      2
    for( unsigned short l_ro = 0; l_ro < 3; l_ro++ ) {
      l_tmp[1][l_ro][0] *= -1;
      l_tmp[1][l_ro][2] *= -1;
    }
  }

  // rotate right solver to face-aligned coordinates
  linalg::Matrix::matMulB0( 3, 5, 5,
                            l_tmp[1][0],
                            i_tm1[0],
                            o_fsAr );
}

void edge::seismic::setups::ViscoElasticity::fs( double       i_rhoL,
                                                 double       i_rhoR,
                                                 double       i_lamL,
                                                 double       i_lamR,
                                                 double       i_muL,
                                                 double       i_muR,
                                                 double const i_tA[6][6],
                                                 double const i_tm1[9][9],
                                                 double       o_fsAl[6*9],
                                                 double       o_fsAr[6*9],
                                                 bool         i_freeSurface ) {
  // intermediate matrices
  double l_fsMidA[2][6][9];
  double l_tmp[2][6][9];

  // compute mid part of the solver
  fsMid( i_rhoL,
         i_rhoR,
         i_lamL,
         i_lamR,
         i_muL,
         i_muR,
         l_fsMidA[0],
         l_fsMidA[1] );

  // apply back-rotation to face-aligned coordinate system and mid part
  for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
    linalg::Matrix::matMulB0( 6, 9, 6,
                              i_tA[0],
                              l_fsMidA[l_sd][0],
                              l_tmp[l_sd][0] );
  }

  // rotate left solver to face-aligned coordinates (element exists per definition of our boundary conditions)
  linalg::Matrix::matMulB0( 6, 9, 9,
                            l_tmp[0][0],
                            i_tm1[0],
                            o_fsAl );

  /*
  * Free surface boundary conditions mirror (rotated) x-components of the stress tensor.
  * Reference: Eq. (41) in
  *            Dumbser, Michael, and Martin Kaeser.
  *            "An arbitrary high-order discontinuous Galerkin method for elastic waves on
  *            unstructured meshes—II. The three-dimensional isotropic case."
  *            Geophysical Journal International 167.1 (2006): 319-336.
  */
  if( i_freeSurface == true ) {
    // multiply right side with diag (-1, 1, 1, -1, 1, -1, 1, 1, 1 ) from the right.
    //                                 0         3      5
    for( unsigned short l_ro = 0; l_ro < 6; l_ro++ ) {
      l_tmp[1][l_ro][0] *= -1;
      l_tmp[1][l_ro][3] *= -1;
      l_tmp[1][l_ro][5] *= -1;
    }
  }

  // rotate right solver to face-aligned coordinates
  linalg::Matrix::matMulB0( 6, 9, 9,
                            l_tmp[1][0],
                            i_tm1[0],
                            o_fsAr );
}
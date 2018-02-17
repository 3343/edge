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
 * Unit tests for friction laws
 **/

#include <catch.hpp>
#include "constants.hpp"
#define private public
#include "FrictionLaws.hpp"
#undef private

#include <iostream>

TEST_CASE( "Perturbation of middle states through linear slip weaking in 2D.", "[FrictionLaws][pertLSW2D]" ) {
  /*
   * right equals minus
   *
   * ql = 1.0E2, -2.0E3, -4.0E5, 10.0, -15.0
   * qr = -0.1E2, -0.9E3,  5.0E5,  4.0,   5.0
   * (see internal boundary for derivation)
   * ms = -48059955, -16042393.1267, 92512100, 6.99999656679, -4.95133141038
   */

  // set up left and right vars plus middle state
  double l_qL[5] = {  1.0E2, -2.0E3, -4.0E5, 10.0, -15.0 };
  double l_qR[5] = { -0.1E2, -0.9E3,  5.0E5,  4.0,   5.0 };
  double l_ms[5][1] = { {-48059955}, {-16042393.1267}, {92512100}, {6.99999656679}, {-4.95133141038} };

  // set up mat vars
  double l_lam = 32080749540.0;
  double l_mu  = 32019625230.0;
  double l_rho = 2670;

  double l_cs  = std::sqrt( l_mu / l_rho );
  double l_cp  = l_lam + 2.0 * l_mu;
         l_cp /= l_rho;
         l_cp  = std::sqrt( l_cp );

  // set up friction law
  double l_dt, l_csDmuM, l_csDmuP;
  double l_mus[1], l_mud[1], l_dcInv[1], l_sn0[1], l_ss0[1], l_muf[1], l_dd[1];
  double l_msM[5][1], l_msP[5][1];

  l_dt = 0.004;

  l_csDmuM = l_cs / l_mu;
  l_csDmuP = l_cs / l_mu;

  l_mus[0]   = 0.677;
  l_mud[0]   = 0.55;
  l_dcInv[0] = 2.5;
  l_sn0[0]   = -120E6;
  l_ss0[0]   = 81.6E6;
  l_muf[0]   = 0.59;
  l_dd[0]    = 0.2;

  // compare the middle states with formulas in
  // Reference:
  //   Dynamic rupture modeling on unstructured meshes using a discontinuous Galerkin method
  //   Puente, Ampuero, Kaeser
  //   Journal of Geophysical Research, Vol. 114
  // Remark: plus equals left

  double l_msRef[2][5];
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    l_msRef[0][l_qt] = l_msRef[1][l_qt] = l_ms[l_qt][0];
  }

  // get strength
  double l_strength = l_muf[0]*(l_ms[0][0] + l_sn0[0]);
  // adjust normal stress (compression is positive in the paper
  l_strength *= -1;

  // set perturbed shear stress
  double l_ssp = l_strength - l_ss0[0];
         l_ssp = std::min( l_ms[2][0], l_ssp );
  l_msRef[0][2] = l_msRef[1][2] = l_ssp;

  // set perturbed velocities
  double l_vp = l_qL[4] + l_cs / l_mu * ( l_ssp - l_qL[2] );
  double l_vm = l_qR[4] - l_cs / l_mu * ( l_ssp - l_qR[2] );

  l_msRef[0][4] = l_vm;
  l_msRef[1][4] = l_vp;

  // get the slip rate
  double l_slipRateRef  = ( 2.0 * l_cs ) / l_mu;
         l_slipRateRef *= ( l_ssp - l_ms[2][0] );
         l_slipRateRef  = std::abs( l_slipRateRef );

  // get the slip
  double l_slipRef = l_dd[0] + l_dt * l_slipRateRef;

  // get the new friction coefficient
  double l_muRef  = (l_mus[0] - l_mud[0]) * l_dcInv[0];
         l_muRef *= l_slipRef;
         l_muRef  = l_mus[0] - l_muRef;
         l_muRef  = std::max( l_mud[0], l_muRef );

  double l_sr[2], l_ss[2];
  bool l_per[1];

  // get solution through our solver
  edge::elastic::solvers::FrictionLaws< 2, 1 >::linSlipWeak( l_dt,
                                                             l_csDmuM, l_csDmuP,
                                                             l_mus, l_mud, l_dcInv,
                                                             l_sn0, l_ss0,
                                                             l_ms,
                                                             l_dd, l_muf,
                                                             l_sr, l_ss,
                                                             l_msM, l_msP, l_per );

  // check that the fault failed
  REQUIRE( l_per[0] );

  // check the perturbed middle states
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    REQUIRE( l_msM[l_qt][0] == Approx( l_msRef[0][l_qt] ) );
    REQUIRE( l_msP[l_qt][0] == Approx( l_msRef[1][l_qt] ) );
  }

  // check slip
  REQUIRE( l_slipRef == Approx( l_dd[0] ) );

  // check friction coefficinet
  REQUIRE( l_muf[0] == Approx(l_muRef) );
}

TEST_CASE( "Perturbation of middle states through linear slip weaking in 3D.", "[FrictionLaws][pertLSW3D]" ) {
  /*
   * left equals plus
   * right equals minus
   *
   * ql = 1.0E6, -2.0E6, -4.0E5, 5.0E5, -2.1E6, 4.9E6, 1.0, -0.2,  1.0
   * qr = -0.1E7, -0.9E6,  5.0E6, 2.0E6,  3.4E6, 9.3E5, 1.2,  0.3, -0.3
   * (see internal boundary for derivation)
   * ms = 1138813.59552691 -1953669.95223839 -353669.952238394 3283577.90819717 -2100000 -2378532.97028903 1.00866501844737 0.101050690844916 0.212809035238327
   */

  double l_qP[9] = {  1.0E6, -2.0E6, -4.0E5, 5.0E5, -2.1E6, 4.9E6, 1.0, -0.2,  1.0 };
  double l_qM[9] = { -0.1E7, -0.9E6,  5.0E6, 2.0E6,  3.4E6, 9.3E5, 1.2,  0.3, -0.3 };
  double l_ms[9][1] = { {1138813.59552691}, {-1953669.95223839}, {-353669.952238394}, // normal stresses
                        {3283577.90819717}, {         -2100000}, {-2378532.97028903}, // shear stresses
                        {1.00866501844737}, {0.101050690844916}, {0.212809035238327}  // particle velocities
                      };

  // set up material parameters
  double l_lamP = 32080749540.0;
  double l_muP  = 32019625230.0;
  double l_rhoP = 2670;

  double l_lamM = l_lamP / 2;
  double l_muM = l_muP / 2;
  double l_rhoM = 2600;

  // set up velocities
  double l_csP  = std::sqrt( l_muP / l_rhoP );

  double l_csM  = std::sqrt( l_muM / l_rhoM );

  double l_cpM  = l_lamM + 2.0 * l_muM;
         l_cpM /= l_rhoM;
         l_cpM  = std::sqrt( l_cpM );

  double l_cpP  = l_lamP + 2.0 * l_muP;
         l_cpP /= l_rhoP;
         l_cpP  = std::sqrt( l_cpP );

  // set up friction law
  double l_dt, l_csDmuM, l_csDmuP;
  double l_mus[1], l_mud[1], l_dcInv[1], l_sn0[1], l_ss0[2][1], l_muf[1], l_dd[2][1];
  double l_msM[9][1], l_msP[9][1];

  l_csDmuM = l_csM / l_muM;
  l_csDmuP = l_csP / l_muP;

  l_dt = 0.004;

  l_mus[0]    = 0.677;
  l_mud[0]    = 0.55;
  l_dcInv[0]  = 2.5;
  l_sn0[0]    = -120E6;
  l_ss0[0][0] = 40E6;
  l_ss0[1][0] = std::sqrt( 81.6E6 * 81.6E6 - 40.0E6 * 40.0E6 );
  l_muf[0]    = 0.59;
  l_dd[0][0]  = 0.2;
  l_dd[1][0]  = -0.1;

  double l_ss0A[1];
         l_ss0A[0] = std::sqrt( l_ss0[0][0] * l_ss0[0][0] + l_ss0[1][0] * l_ss0[1][0] );

  // compare the perturbed middle states with explicit formulas in
  // Reference: Pelties, Gabriel, Ampuero
  //            Verification of an ADER-DG method for complex dynamic rupture problems.
  //            Geosci. Model Dev., 7, 2014

  double l_msRef[2][9];
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
   l_msRef[0][l_qt] = l_ms[l_qt][0];
   l_msRef[1][l_qt] = l_ms[l_qt][0];
  }

  double l_strength  = l_muf[0] * ( std::abs( l_ms[0][0] + l_sn0[0] ) );
  double l_ssA       = ( l_ms[3][0] + l_ss0[0][0] ) * ( l_ms[3][0] + l_ss0[0][0] );
         l_ssA      += ( l_ms[5][0] + l_ss0[1][0] ) * ( l_ms[5][0] + l_ss0[1][0] );
         l_ssA       = std::sqrt( l_ssA );

  // set perturbed shear stresses
  for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
    l_msRef[l_sd][3]  = l_msRef[l_sd][3] + l_ss0[0][0];
    l_msRef[l_sd][3] *= l_strength / l_ssA;
    l_msRef[l_sd][3] -= l_ss0[0][0];

    l_msRef[l_sd][5]  = l_msRef[l_sd][5] + l_ss0[1][0];
    l_msRef[l_sd][5] *= l_strength / l_ssA;
    l_msRef[l_sd][5] -= l_ss0[1][0];
  }

  // set perturbed fault-parallel particle velocities
  l_msRef[0][7] = l_qP[7] + ( l_msRef[0][3] - l_qP[3] ) / ( l_csP * l_rhoP );
  l_msRef[0][8] = l_qP[8] + ( l_msRef[0][5] - l_qP[5] ) / ( l_csP * l_rhoP );

  l_msRef[1][7] = l_qM[7] - ( l_msRef[0][3] - l_qM[3] ) / ( l_csM * l_rhoM );
  l_msRef[1][8] = l_qM[8] - ( l_msRef[0][5] - l_qM[5] ) / ( l_csM * l_rhoM );

  // set slip rates
  double l_srRef[2];
  l_srRef[0]  = 1.0 / ( l_csP * l_rhoP ) + 1.0 / ( l_csM * l_rhoM );
  l_srRef[0] *= ( l_msRef[0][3] - l_ms[3][0] );

  l_srRef[1]  = 1.0 / ( l_csP * l_rhoP ) + 1.0 / ( l_csM * l_rhoM );
  l_srRef[1] *= ( l_msRef[0][5] - l_ms[5][0] );

  // call our friction solvers
  double l_sr[2][1], l_tr[2][1];
  bool l_per[1];

  edge::elastic::solvers::FrictionLaws<
    3,
    1
   >::linSlipWeak( l_dt,
                   l_csDmuM,
                   l_csDmuP,
                   l_mus,
                   l_mud,
                   l_dcInv,
                   l_sn0,
                   l_ss0,
                   l_ss0A,
                   l_ms,
                   l_dd,
                   l_muf,
                   l_sr,
                   l_tr,
                   l_msM,
                   l_msP,
                   l_per );

  // check that the fault failed
  REQUIRE( l_per[0] );

  // compare the perturbed middle states
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    REQUIRE( l_msP[l_qt][0] == Approx( l_msRef[0][l_qt] ) );
    REQUIRE( l_msM[l_qt][0] == Approx( l_msRef[1][l_qt] ) );
  }

  // compare the slip rates
  REQUIRE( l_sr[0][0] == Approx( l_srRef[0] ) );
  REQUIRE( l_sr[1][0] == Approx( l_srRef[1] ) );
}

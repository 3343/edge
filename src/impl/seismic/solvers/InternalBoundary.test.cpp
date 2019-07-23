/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Unit tests for internal boundary conditions at faces through quadrature rules.
 **/
#include <catch.hpp>
#include "TimePred.hpp"
#define private public
#include "common.hpp"
#undef private
#include "../common.hpp"
#define private public
#include "InternalBoundary.hpp"
#undef private
#include "FrictionLaws.hpp"

TEST_CASE( "Internal Boundary: Middle state jump 3D.", "[intBnd][jump3D]" ) {
 /**
  * Input:
  *  3588.3       1210
  * -18768.3     -12637.3
  * -3733.8      -5187.24
  * -147941      -158914
  * -438.478     -3754.98
  * -260.437      1846.82
  * -0.000583019  0.000496199
  * -0.0145165    0.0165153
  * -2.26128e-05 -0.000285888
  *
  * rho:    2670
  * lambda: 32043759360
  * mu:     32038120320
  *
  * Derivation in Mathematica:
  *   strengths =  Rm1.(( {
  *          {3588.3},
  *          {-18768.3},
  *          {-3733.8},
  *          {-147941.},
  *          {-438.478},
  *          {-260.437},
  *          {-0.000583019},
  *          {-0.0145165},
  *          {-2.26128*10^-5}
  *         } ) - ( {
  *          {1210},
  *          {-12637.3},
  *          {-5187.24},
  *          {-158914},
  *          {-3754.98},
  *           {1846.82},
  *          {0.000496199},
  *          {0.0165153},
  *          {-0.000285888}
  *         } ));
  *   strengths = 
  *     strengths /. 
  *      Thread[{laml, lamr, mul, mur, rhol, rhor} -> {32043759360., 
  *         32043759360., 32038120320., 32038120320., 2670., 2670}];
  *   eigenvectors = 
  *     R /. Thread[{laml, lamr, mul, mur, rhol, rhor} -> {32043759360., 
  *         32043759360., 32038120320., 32038120320., 2670., 2670}];
  *   strengths[[1, 1]]*eigenvectors[[All, 1]] + 
  *    strengths[[2, 1]]*eigenvectors[[All, 2]] + 
  *    strengths[[3, 1]]*eigenvectors[[All, 3]]
  *
  * Output:
  * {-7455.39, -2485.42, -2485.42, -138018., 0., 163.872, -0.00046538, \
-0.0149227, 0.000017718}
  **/

  double l_lam = 32043759360;
  double l_mu  = 32038120320;
  double l_rho = 2670;


  double l_msJump1[9][9];

  edge::seismic::solvers::InternalBoundarySolvers<TET4>::solvMsJumpL( l_lam, l_lam,
                                                                      l_mu,  l_mu,
                                                                      l_rho, l_rho,
                                                                      l_msJump1 );

  // set left and right
  double l_qL[9] = {3588.3, -18768.3, -3733.8, -147941, -438.478, -260.437, -0.000583019, -0.0145165, -2.26128e-05};
  double l_qR[9] = {1210, -12637.3, -5187.24, -158914, -3754.98, 1846.82, 0.000496199, 0.0165153, -0.000285888  };

  // compute jump
  double l_jump[9];
  for( unsigned short l_qt =0; l_qt < 9; l_qt++ ) {
    l_jump[l_qt] = l_qL[l_qt] - l_qR[l_qt];
  }

  // derive the update of the middle state solver
  double l_update[9];
  edge::linalg::Matrix::matMulB0( 9, 1, 9,
                                  l_msJump1[0], l_jump, l_update );

  // compare the results
  REQUIRE( l_update[0] == Approx( -7455.39    ) );
  REQUIRE( l_update[1] == Approx( -2485.42    ) );
  REQUIRE( l_update[2] == Approx( -2485.42    ) );
  REQUIRE( l_update[3] == Approx( -138018.    ) );
  REQUIRE( l_update[4] == Approx(  0.         ) );
  REQUIRE( l_update[5] == Approx(  163.872    ) );
  REQUIRE( l_update[6] == Approx( -0.00046538 ) );
  REQUIRE( l_update[7] == Approx( -0.0149227  ) );
  REQUIRE( l_update[8] == Approx( 0.000017718 ) );
}

TEST_CASE( "Internal Boundary: 2D solvers with middle state perturbation.", "[intBnd][solvers2D]" ) {
  // middle states
  double l_msJ[5][5];
  // flux computation
  double l_msF[5][5];

  double l_lam = 32080749540.0;
  double l_mu  = 32019625230.0;
  double l_rho = 2670;

  // get the solvers
  edge::seismic::solvers::InternalBoundarySolvers<
    TRIA3 >::solvMsJumpL( l_lam, l_lam,
                          l_mu,  l_mu,
                          l_rho, l_rho,
                          l_msJ );

  edge::seismic::solvers::InternalBoundarySolvers<
    TRIA3 >::solvMsFlux( 1.0, 0.0, // no trafo
                         l_lam, l_mu, l_rho,
                         l_msF );

  // quantities on the left and right side
  double l_qL[5] = {  1.0E2, -2.0E3, -4.0E5, 10.0, -15.0 };
  double l_qR[5] = { -0.1E2, -0.9E3,  5.0E5,  4.0,   5.0 };

  // jump in quantities
  double l_jump[5];
  for( unsigned l_qt = 0; l_qt < 5; l_qt++ ) {
    l_jump[l_qt] = l_qR[l_qt] - l_qL[l_qt];
  }

  // get the middle state
  double l_ms[5] = { 1.0E2, -2.0E3, -4.0E5, 10.0, -15.0 };
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_ms[l_q1] += l_msJ[l_q1][l_q2] * l_jump[l_q2];
    }
  }

  // compare the middle states with explicit formulas in
  // Reference:
  //   Dynamic rupture modeling on unstructured meshes using a discontinuous Galerkin method
  //   Puente, Ampuero, Kaeser
  //   Journal of Geophysical Research, Vol. 114
  // Remark: plus equals left
  double l_cs  = std::sqrt( l_mu / l_rho );
  double l_cp  = l_lam + 2.0 * l_mu;
         l_cp /= l_rho;
         l_cp  = std::sqrt( l_cp );
  double l_msRef[5];

  l_msRef[0] =  ( l_qR[0] + l_qL[0] );
  l_msRef[0] += ( (l_lam + 2.0 * l_mu) / l_cp ) * ( l_qR[3] - l_qL[3] );
  l_msRef[0] /= 2;

  l_msRef[1]  = (l_lam / l_cp) * ( l_qR[3] - l_qL[3] );
  l_msRef[1] += ( l_lam / ( l_lam + 2.0 * l_mu ) ) * ( l_qR[0] + l_qR[1] );
  l_msRef[1] += 2.0 * l_qL[1];
  l_msRef[1] /= 2;

  l_msRef[2]  = (l_qR[2] + l_qL[2] );
  l_msRef[2] += ( l_mu / l_cs ) * ( l_qR[4] - l_qL[4] );
  l_msRef[2] /= 2;

  l_msRef[3]  = l_qR[3] + l_qL[3];
  l_msRef[3] += ( l_cp / (l_lam + 2.0 * l_mu) ) * ( l_qR[0] - l_qL[0] );
  l_msRef[3] /= 2;

  l_msRef[4]  = l_qR[4] + l_qL[4];
  l_msRef[4] += ( l_cs / l_mu ) * ( l_qR[2] - l_qL[2] );
  l_msRef[4] /= 2;

  // check the middle states
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    REQUIRE( l_ms[l_qt] == Approx( l_msRef[l_qt] ) );
  }

  // calling the friction solver with + equals left, gives the following perturbed middle states
  // minus:         plus:
  // -4.806e+07    -4.806e+07
  // -1.60424e+07  -1.60424e+07
  // 1.75554e+07    1.75554e+07
  // 7              7
  // 3.15542        -13.0581
  //

  double l_msLp[5] = { -4.806e+07, -1.60424e+07, 1.75554e+07, 7, -13.0581 };
  double l_msRp[5] = { -4.806e+07, -1.60424e+07, 1.75554e+07, 7,  3.15542 };

  // the flux is obtained by multiplying the perturbed middle states with the Jacobian
  double l_fluxL[5] = {0,0,0,0,0};
  double l_fluxR[5] = {0,0,0,0,0};

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_fluxL[l_q1] += l_msF[l_q1][l_q2] * l_msLp[l_q2];
      l_fluxR[l_q1] += l_msF[l_q1][l_q2] * l_msRp[l_q2];
    }
  }

  // set up reference jacobian
  double l_jacRef[5][5] = { {0,0,0,0,0},
                            {0,0,0,0,0},
                            {0,0,0,0,0},
                            {0,0,0,0,0},
                            {0,0,0,0,0} };
  l_jacRef[0][3] = -(l_lam + 2.0 * l_mu);
  l_jacRef[1][3] = -l_lam;
  l_jacRef[2][4] = -l_mu;
  l_jacRef[3][0] = -1.0 / l_rho;
  l_jacRef[4][2] = -1.0 / l_rho;

  // check the jacobian with the flux solver
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_jacRef[l_q1][l_q2] == Approx( l_msF[l_q1][l_q2] ) );
    }
  }

  // multiply with jacobian and get ref flux (yep, that's totally redundant)
  double l_fluxLref[5] = {0,0,0,0,0};
  double l_fluxRref[5] = {0,0,0,0,0};

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_fluxLref[l_q1] += l_jacRef[l_q1][l_q2] * l_msLp[l_q2];
      l_fluxRref[l_q1] += l_jacRef[l_q1][l_q2] * l_msRp[l_q2];
    }
  }

  // check the fluxes
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    REQUIRE( l_fluxL[l_qt] == l_fluxLref[l_qt] );
    REQUIRE( l_fluxR[l_qt] == l_fluxRref[l_qt] );
  }

  /*
   * Do all of the above by calling the integrated friction solver.
   */
  double l_dt = 0.004;
  double l_tm1[5][5] = { {1,0,0,0,0},
                         {0,1,0,0,0},
                         {0,0,1,0,0},
                         {0,0,0,1,0},
                         {0,0,0,0,1} };
  double l_dofsL[5][1][1] = { {{ 1.0E2}}, {{-2.0E3}}, {{-4.0E5}}, {{10.0}}, {{-15.0}} };
  double l_dofsR[5][1][1] = { {{-0.1E2}}, {{-0.9E3}}, {{ 5.0E5}}, {{ 4.0}}, {{  5.0}} };

  edge::seismic::solvers::t_LinSlipWeakGlobal<double,1>      l_gl;
  edge::seismic::solvers::t_LinSlipWeakFace<double>          l_fa;
  edge::seismic::solvers::t_LinSlipWeakSubFace<double, 2, 1> l_sf[1];

  // assign the values of our example
  l_gl.mus[0]   = 0.677;
  l_gl.mud[0]   = 0.55;
  l_gl.dcInv[0] = 2.5;

  l_fa.lEqM = false;
  l_fa.csDmuM = l_cs / l_mu;
  l_fa.csDmuP = l_cs / l_mu;

  l_sf[0].sn0[0]    = -120E6;
  l_sf[0].ss0[0][0] = 81.6E6;
  l_sf[0].muf[0]    = 0.59;
  l_sf[0].dd[0][0]  = 0.2;

  // generate a struct out of it
  struct {
    edge::seismic::solvers::t_LinSlipWeakGlobal<double,1>       *gl;
    edge::seismic::solvers::t_LinSlipWeakFace<double>           *fa;
    edge::seismic::solvers::t_LinSlipWeakSubFace<double, 2, 1> (*sf)[1];
  } l_faData;
  l_faData.gl = &l_gl;
  l_faData.fa = &l_fa;
  l_faData.sf = &l_sf;

  // call the solver
  double l_netUpsL[5][1][1];
  double l_netUpsR[5][1][1];
  bool l_per[1][1];

  edge::seismic::solvers::InternalBoundary<
    TRIA3,
    5,
    1,
    1 >::netUpdates<
           double,
           edge::seismic::solvers::FrictionLaws< 2, 1 >
         > (  l_tm1,
              l_msJ,
              l_msF,
              l_msF,
              l_dofsL,
              l_dofsR,
              l_netUpsL,
              l_netUpsR,
              l_per,
              1.0,
             &l_faData );

  // check for perturbation
  REQUIRE( l_per[0][0] == true );

  // check the updates
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    REQUIRE( l_netUpsL[l_qt][0][0] == Approx(-l_fluxLref[l_qt]) );
    REQUIRE( l_netUpsR[l_qt][0][0] == Approx( l_fluxRref[l_qt]) );
  }
}

TEST_CASE( "Internal Boundary: 3D solvers with middle state perturbation.", "[intBnd][solvers3D]" ) {
  // middle states
  double l_msJ[9][9];
  // flux computation
  double l_msFL[9][9];
  double l_msFR[9][9];

  double l_lamL = 32080749540.0;
  double l_muL  = 32019625230.0;
  double l_rhoL = 2670;

  double l_lamR = l_lamL / 2;
  double l_muR = l_muL / 2;
  double l_rhoR = 2600;

  // get the solvers
  edge::seismic::solvers::InternalBoundarySolvers<
    TET4
  >::solvMsJumpL( l_lamL, l_lamR,
                  l_muL,  l_muR,
                  l_rhoL, l_rhoR,
                  l_msJ );

  edge::seismic::solvers::InternalBoundarySolvers<
    TET4
   >::solvMsFlux( 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0,
                  l_lamL, l_muL, l_rhoL,
                  l_msFL );
  edge::seismic::solvers::InternalBoundarySolvers<
    TET4
   >::solvMsFlux( 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0,
                  l_lamR, l_muR, l_rhoR,
                  l_msFR );

  // quantities on the left and right side
  double l_qL[9] = {  1.0E6, -2.0E6, -4.0E5, 5.0E5, -2.1E6, 4.9E6, 1.0, -0.2,  1.0 };
  double l_qR[9] = { -0.1E7, -0.9E6,  5.0E6, 2.0E6,  3.4E6, 9.3E5, 1.2,  0.3, -0.3 };

  // jump in quantities
  double l_jump[9];
  for( unsigned l_qt = 0; l_qt < 9; l_qt++ ) {
    l_jump[l_qt] = l_qR[l_qt] - l_qL[l_qt];
  }

  // get the middle state
  double l_ms[9] = {  1.0E6, -2.0E6, -4.0E5, 5.0E5, -2.1E6, 4.9E6, 1.0, -0.2,  1.0 };
  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      l_ms[l_q1] += l_msJ[l_q1][l_q2] * l_jump[l_q2];
    }
  }

  // compare the middle states with explicit formulas in
  // Reference: Pelties, Gabriel, Ampuero
  //            Verification of an ADER-DG method for complex dynamic rupture problems.
  //            Geosci. Model Dev., 7, 2014
  // Remark: plus equals left
  double l_csL  = std::sqrt( l_muL / l_rhoL );

  double l_csR  = std::sqrt( l_muR / l_rhoR );

  double l_cpL  = l_lamL + 2.0 * l_muL;
         l_cpL /= l_rhoL;
         l_cpL  = std::sqrt( l_cpL );

  double l_cpR  = l_lamR + 2.0 * l_muR;
         l_cpR /= l_rhoR;
         l_cpR  = std::sqrt( l_cpR );

  double l_msRef[9];

  l_msRef[0]  = l_qR[0] - l_qL[0] + l_cpR * l_rhoR * ( l_qR[6] - l_qL[6] );
  l_msRef[0] *= l_cpL * l_rhoL;
  l_msRef[0] /= l_cpL * l_rhoL + l_cpR * l_rhoR;
  l_msRef[0] += l_qL[0];

  l_msRef[3]  = l_qR[3] - l_qL[3] + l_csR * l_rhoR * ( l_qR[7] - l_qL[7] );
  l_msRef[3] *= l_csL * l_rhoL;
  l_msRef[3] /= l_csL * l_rhoL + l_csR * l_rhoR;
  l_msRef[3] += l_qL[3];

  l_msRef[5]  = l_qR[5] - l_qL[5] + l_csR * l_rhoR * ( l_qR[8] - l_qL[8] );
  l_msRef[5] *= l_csL * l_rhoL;
  l_msRef[5] /= l_csL * l_rhoL + l_csR * l_rhoR;
  l_msRef[5] += l_qL[5];

  l_msRef[6]  = l_msRef[0] - l_qL[0];
  l_msRef[6] /= l_cpL * l_rhoL;
  l_msRef[6] += l_qL[6];

  l_msRef[7]  = l_msRef[3] - l_qL[3];
  l_msRef[7] /= l_csL * l_rhoL;
  l_msRef[7] += l_qL[7];

  l_msRef[8]  = l_msRef[5] - l_qL[5];
  l_msRef[8] /= l_csL * l_rhoL;
  l_msRef[8] += l_qL[8];

  REQUIRE( l_ms[0] == Approx( l_msRef[0] ) );
  REQUIRE( l_ms[3] == Approx( l_msRef[3] ) );
  REQUIRE( l_ms[5] == Approx( l_msRef[5] ) );
  REQUIRE( l_ms[6] == Approx( l_msRef[6] ) );
  REQUIRE( l_ms[7] == Approx( l_msRef[7] ) );
  REQUIRE( l_ms[8] == Approx( l_msRef[8] ) );

  //
  // After perturbation of the middle states through the 3D linear slip weakening
  // friciton law with the paramters below, gives the following middle states
  // left:              right
  // 1138813.59552691   1138813.59552691
  // -1953669.95223839 -1953669.95223839
  // -353669.952238394 -353669.952238394
  // -2635018.36571765 -2635018.36571765
  // -2100000          -2100000
  // -11778726.9409157 -11778726.9409157
  // 1.00866501844737   1.00866501844737
  // -0.539059827293307 1.01840882886168
  // -0.803844704037184 1.66980053099153


  double l_msLp[9] = {  1138813.59552691, -1953669.95223839,  -353669.952238394,    // normal stresses
                       -2635018.36571765, -2100000,           -11778726.9409157,    // shear stresses
                        1.00866501844737, -0.539059827293307, -0.803844704037184 }; // particle velocities

  double l_msRp[9] = {  1138813.59552691, -1953669.95223839,  -353669.952238394,    // normal stresses
                       -2635018.36571765, -2100000,           -11778726.9409157,    // shear stresses
                        1.00866501844737,  1.01840882886168,   1.66980053099153 };  // particle velocities


  // the flux is obtained by multiplying the perturbed middle states with the Jacobian
  double l_fluxL[9] = {0,0,0, 0,0,0, 0,0,0};
  double l_fluxR[9] = {0,0,0, 0,0,0, 0,0,0};

  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      l_fluxL[l_q1] += l_msFL[l_q1][l_q2] * l_msLp[l_q2];
      l_fluxR[l_q1] += l_msFR[l_q1][l_q2] * l_msRp[l_q2];
    }
  }

  // set up reference jacobians
  double l_jacRefL[9][9] = { {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0} };
  l_jacRefL[0][6] = -(l_lamL + 2.0 * l_muL);
  l_jacRefL[1][6] = -l_lamL;
  l_jacRefL[2][6] = -l_lamL;

  l_jacRefL[3][7] = -l_muL;
  l_jacRefL[5][8] = -l_muL;

  l_jacRefL[6][0] = -1.0 / l_rhoL;
  l_jacRefL[7][3] = -1.0 / l_rhoL;
  l_jacRefL[8][5] = -1.0 / l_rhoL;

  double l_jacRefR[9][9] = { {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0},
                             {0,0,0, 0,0,0, 0,0,0} };
  l_jacRefR[0][6] = -(l_lamR + 2.0 * l_muR);
  l_jacRefR[1][6] = -l_lamR;
  l_jacRefR[2][6] = -l_lamR;

  l_jacRefR[3][7] = -l_muR;
  l_jacRefR[5][8] = -l_muR;

  l_jacRefR[6][0] = -1.0 / l_rhoR;
  l_jacRefR[7][3] = -1.0 / l_rhoR;
  l_jacRefR[8][5] = -1.0 / l_rhoR;


  // check the jacobian with the flux solver
  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      REQUIRE( l_jacRefL[l_q1][l_q2] == Approx( l_msFL[l_q1][l_q2] ) );
      REQUIRE( l_jacRefR[l_q1][l_q2] == Approx( l_msFR[l_q1][l_q2] ) );
    }
  }


  // multiply with jacobian and get ref flux (yep, that's totally redundant, once again..)
  double l_fluxLref[9] = {0,0,0, 0,0,0, 0,0,0};
  double l_fluxRref[9] = {0,0,0, 0,0,0, 0,0,0};

  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      l_fluxLref[l_q1] += l_jacRefL[l_q1][l_q2] * l_msLp[l_q2];
      l_fluxRref[l_q1] += l_jacRefR[l_q1][l_q2] * l_msRp[l_q2];
    }
  }

  // check the fluxes
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    REQUIRE( l_fluxL[l_qt] == l_fluxLref[l_qt] );
    REQUIRE( l_fluxR[l_qt] == l_fluxRref[l_qt] );
  }

  /*
   * Do all of the above by calling the integrated friction solver.
   */
  double l_dt = 0.004;
  double l_tm1[9][9] = { {1,0,0,0,0,0,0,0,0},
                         {0,1,0,0,0,0,0,0,0},
                         {0,0,1,0,0,0,0,0,0},
                         {0,0,0,1,0,0,0,0,0},
                         {0,0,0,0,1,0,0,0,0},
                         {0,0,0,0,0,1,0,0,0},
                         {0,0,0,0,0,0,1,0,0},
                         {0,0,0,0,0,0,0,1,0},
                         {0,0,0,0,0,0,0,0,1} };
  double l_dofsL[9][1][1] = { {{ 1.0E6}}, {{-2.0E6}}, {{-4.0E5}},
                              {{ 5.0E5}}, {{-2.1E6}}, {{ 4.9E6}},
                              {{   1.0}}, {{-0.2  }}, {{   1.0}} };
  double l_dofsR[9][1][1] = { {{-0.1E7}}, {{-0.9E6}}, {{ 5.0E6}},
                              {{ 2.0E6}}, {{ 3.4E6}}, {{ 9.3E5}},
                              {{   1.2}}, {{   0.3}}, {{  -0.3}} };

  edge::seismic::solvers::t_LinSlipWeakGlobal<double,1>      l_gl;
  edge::seismic::solvers::t_LinSlipWeakFace<double>          l_fa;
  edge::seismic::solvers::t_LinSlipWeakSubFace<double, 3, 1> l_sf[1];


  // assign the values of our example
  l_gl.mus[0]    = 0.677;
  l_gl.mud[0]    = 0.55;
  l_gl.dcInv[0]  = 2.5;

  l_fa.lEqM = false;
  l_fa.csDmuP = l_csL / l_muL;
  l_fa.csDmuM = l_csR / l_muR;

  l_sf[0].sn0[0]    = -120E6;
  l_sf[0].ss0[0][0] = 40E6;
  l_sf[0].ss0[1][0] = std::sqrt( 81.6E6 * 81.6E6 - 40.0E6 * 40.0E6 );
  l_sf[0].ss0A[0] = std::sqrt(   l_sf[0].ss0[0][0] * l_sf[0].ss0[0][0]
                               + l_sf[0].ss0[1][0] * l_sf[0].ss0[1][0] );
  l_sf[0].muf[0]    = 0.59;
  l_sf[0].dd[0][0]  = 0.2;
  l_sf[0].dd[1][0]  = -0.1;

  // generate a struct out of that
  struct {
    edge::seismic::solvers::t_LinSlipWeakGlobal<double,1>       *gl;
    edge::seismic::solvers::t_LinSlipWeakFace<double>           *fa;
    edge::seismic::solvers::t_LinSlipWeakSubFace<double, 3, 1> (*sf)[1];
  } l_faData;
  l_faData.gl = &l_gl;
  l_faData.fa = &l_fa;
  l_faData.sf = &l_sf;

  // call the solver
  double l_netUpsL[9][1][1];
  double l_netUpsR[9][1][1];
  bool l_per[1][1];

  edge::seismic::solvers::InternalBoundary<
    TET4,
    9,
    1,
    1 >::netUpdates<
           double,
           edge::seismic::solvers::FrictionLaws< 3, 1 >
         > (  l_tm1,
              l_msJ,
              l_msFL,
              l_msFR,
              l_dofsL,
              l_dofsR,
              l_netUpsL,
              l_netUpsR,
              l_per,
              1.0,
             &l_faData );

  // check for perturbation
  REQUIRE( l_per[0][0] == true );

  // check the updates
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    REQUIRE( l_netUpsL[l_qt][0][0] == Approx(-l_fluxLref[l_qt]) );
    REQUIRE( l_netUpsR[l_qt][0][0] == Approx( l_fluxRref[l_qt]) );
  }
}

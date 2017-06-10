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

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64


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

  edge::elastic::solvers::InternalBoundarySolvers<TET4>::solvMsJumpL( l_lam, l_lam,
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

TEST_CASE( "Eval of internal boundary conditions through quadrature rules, quad4r", "[intBnd][quad4rEval]" ) {
  // face- and vertex-ids
  unsigned short l_faIdL, l_faIdR, l_vIdR;

  // normals and tangents
  double l_n[2];

  // lame paramters
  double l_lam[2];
  double l_mu[2];
  double l_rho[2];

  /*
   * Quadrature-free reference solution of first unit test:
   *   element type: quad4r
   *   order: 2 -> 4 element modes
   *   quantities: 5
   *   concurrent forward runs: 2
   *   lame parameters left element:  2.2, 1.3, 1.1
   *   lame parameters right element: 2.3, 1.4, 1.0
   */

  // set vertex id to only possible value
  l_vIdR  = 0;

  // set up lame parameters and density
  l_lam[0] = 2.2; l_mu[0] = 1.3; l_rho[0] = 1.1;
  l_lam[1] = 2.3; l_mu[1] = 1.4; l_rho[1] = 1.0;

  // set up normal
  l_n[0] = 0.5;
  l_n[1] = -0.8660254038;

  // set up shared basis
  edge::dg::Basis l_basis1( QUAD4R, 2 );

  // set up quad-free flux matrices
  double l_fluxMatrices1[CE_N_FLUX_MATRICES(QUAD4R)][4][4];
  l_basis1.getFluxMm1Dense( 4, (double *) l_fluxMatrices1, true );

  // set up quadrature-free flux solvers
  double l_fluxSolversL1[2][5][5];
  double l_fluxSolversR1[2][5][5];
  edge::elastic::solvers::common::setupSolver2d( l_rho[0], l_rho[1],
                                                 l_lam[0], l_lam[1],
                                                 l_mu[ 0], l_mu[ 1],
                                                 l_n[0], l_n[1], 0,
                                                 l_fluxSolversL1[0], l_fluxSolversL1[1] );

  edge::elastic::solvers::common::setupSolver2d(  l_rho[1], l_rho[0],
                                                  l_lam[1], l_lam[0],
                                                  l_mu[ 1],  l_mu[0],
                                                 -l_n[0], -l_n[1], 0,
                                                  l_fluxSolversR1[0], l_fluxSolversR1[1] );

  // set up DOFs
  double l_dofsL1[5][4][2];
  double l_dofsR1[5][4][2];

  for( unsigned short l_q = 0; l_q < 5; l_q++ ) {
    for( unsigned short l_md = 0; l_md < 4; l_md++ ) {
      for( unsigned short l_ru = 0; l_ru < 2; l_ru++ ) {
        l_dofsL1[l_q][l_md][l_ru] = (-0.5 - l_q + l_md) * (l_ru+1);
        l_dofsR1[l_q][l_md][l_ru] = ( 0.5 + l_q - l_md) * (l_ru+1);
      }
    }
  }

  // trafo from physical to face-coordinates
  double l_trafoInv[5][5];
  edge::elastic::common::setupTrafoInv2d( l_n[0], l_n[1], l_trafoInv );

  // middle state solver
  double l_msJump1[5][5];
  edge::elastic::solvers::InternalBoundarySolvers<QUAD4R>::solvMsJumpL( l_lam[0], l_lam[1],
                                                                        l_mu[0],  l_mu[1],
                                                                        l_rho[0], l_rho[1],
                                                                        l_msJump1 );

  // flux solvers
  double l_msFluxL1[5][5];
  double l_msFluxR1[5][5];
  edge::elastic::solvers::InternalBoundarySolvers<QUAD4R>::solvMsFlux( l_n[0], l_n[1],
                                                                       l_lam[0], l_mu[0], l_rho[0],
                                                                       l_msFluxL1 );
  edge::elastic::solvers::InternalBoundarySolvers<QUAD4R>::solvMsFlux( l_n[0], l_n[1],
                                                                       l_lam[1], l_mu[1], l_rho[1],
                                                                       l_msFluxR1 );

  // change sign of flux solvers, which don't have a build-in minus for internal boundaries
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_msFluxL1[l_q1][l_q2] *= -1;
      l_msFluxR1[l_q1][l_q2] *= -1;
    }
  }

  // get inverse mass matrix
  double l_massInvDense1[4][4];
  l_basis1.getMassInvDense( 4, l_massInvDense1[0], true );
  double l_massInvDiag1[4];
  for( int_md l_md = 0; l_md < 4; l_md++ ) l_massInvDiag1[l_md] = l_massInvDense1[l_md][l_md];

  // get face-quad points, weights and evaluated basis
  double l_ptsFaces1[8][2][2];
  double l_weightsFaces1[2];
  double l_basisFaces1[8][2][4];
  edge::dg::QuadratureEval< QUAD4R, 2 >::faces( l_ptsFaces1,
                                                l_weightsFaces1,
                                                l_basisFaces1 );

  /*
   * Check for quadrature in space
   */
  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 4; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR = 0; l_faIdR < 4; l_faIdR++ ) {
      // compute reference solution
      double l_refTmp1[5][4][2];
      double l_refUpdateL1[5][4][2];
      double l_refUpdateR1[5][4][2];


      // left element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_dofsL1[0][0],
                                             l_fluxMatrices1[l_faIdL][0],
                                             l_refTmp1[0][0] );

      // left element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversL1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );

      // left element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_dofsR1[0][0],
                                             l_fluxMatrices1[4+l_faIdL*4+l_faIdR][0],
                                             l_refTmp1[0][0] );

      // left element, neighboring contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversL1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );


      // right element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_dofsR1[0][0],
                                             l_fluxMatrices1[l_faIdR][0],
                                             l_refTmp1[0][0] );

      // right element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversR1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      // right element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_dofsL1[0][0],
                                             l_fluxMatrices1[4+l_faIdR*4+l_faIdL][0],
                                             l_refTmp1[0][0] );

      // right element, neighboring contribution multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversR1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      /*
       * Solution through numerical quadrature
       */
      // compute solution through numerical quadrature
      double l_quadUpdateL1[5][4][2];
      double l_quadUpdateR1[5][4][2];

      edge::elastic::solvers::InternalBoundary<QUAD4R, 5, 2, 2>::evalSpace( l_faIdL, l_faIdR, l_vIdR,
                                                                            l_massInvDiag1, l_weightsFaces1, l_basisFaces1,
                                                                            l_trafoInv,
                                                                            l_msJump1,
                                                                            l_msFluxL1, l_msFluxR1,
                                                                            l_dofsL1, l_dofsR1,
                                                                            l_quadUpdateL1, l_quadUpdateR1 );

      /*
       * Compare the two solutions
       */
      for( int_qt l_qt = 0; l_qt < 5; l_qt++ ) {
        for( int_md l_md = 0; l_md < 4; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < 2; l_ru++ ) {
            REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
            REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
          }
        }
      }

    } // remote face id
  } // local face id

  /*
   * Now check for quadrature in space and time.
   */
  // time integrated dofs
  double l_tIntL1[5][4][2];
  double l_tIntR1[5][4][2];
  // time derivatives
  double l_tDerL1[2][5][4][2];
  double l_tDerR1[2][5][4][2];
  // time step
  double l_dT1 = 0.4572;
  // scratch memory
  double l_scatchMem1[4][5][4][2];

  // set up star matrix
  double l_starL1[2][5][5];
  double l_starR1[2][5][5];
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
    for( unsigned short l_qt1 = 0; l_qt1 < 5; l_qt1++ ) {
      for( unsigned short l_qt2 = 0; l_qt2 < 5; l_qt2++ ) {
        l_starL1[l_di][l_qt1][l_qt2] = l_di - l_qt1 * l_qt2;
        l_starR1[l_di][l_qt1][l_qt2] = l_di * l_qt1 - l_qt2 - 0.5;
      }
    }
  }

  // set up quad-free transposed stiffness matrices (mutiplied with inverse mass matrix)
  double l_stiffT1[2][4][4];
  l_basis1.getStiffMm1Dense( 4, (double *) l_stiffT1, true,  true );

  // perform time predictions
  edge::elastic::solvers::TimePred< QUAD4R, 5, 2, 2 >::ckVanilla( l_dT1,
                                                                  l_stiffT1,
                                                                  l_starL1,
                                                                  l_dofsL1,
                                                                  l_scatchMem1[0],
                                                                  l_tDerL1,
                                                                  l_tIntL1 );

  edge::elastic::solvers::TimePred< QUAD4R, 5, 2, 2 >::ckVanilla( l_dT1,
                                                                  l_stiffT1,
                                                                  l_starR1,
                                                                  l_dofsR1,
                                                                  l_scatchMem1[0],
                                                                  l_tDerR1,
                                                                  l_tIntR1 );

  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 4; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR =0; l_faIdR < 4; l_faIdR++ ) {
      // compute reference solution
      double l_refTmp1[5][4][2];
      double l_refUpdateL1[5][4][2];
      double l_refUpdateR1[5][4][2];


      // left element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_tIntL1[0][0],
                                             l_fluxMatrices1[l_faIdL][0],
                                             l_refTmp1[0][0] );

      // left element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversL1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );

      // left element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_tIntR1[0][0],
                                             l_fluxMatrices1[4+l_faIdL*4+l_faIdR][0],
                                             l_refTmp1[0][0] );

      // left element, neighboring contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversL1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );


      // right element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_tIntR1[0][0],
                                             l_fluxMatrices1[l_faIdR][0],
                                             l_refTmp1[0][0] );

      // right element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversR1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      // right element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 2,
                                             5, 4, 4,
                                             l_tIntL1[0][0],
                                             l_fluxMatrices1[4+l_faIdR*4+l_faIdL][0],
                                             l_refTmp1[0][0] );

      // right element, neighboring contribution multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 2,
                                             5, 4, 5,
                                             l_fluxSolversR1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      /*
       * Solution through numerical quadrature
       */
      // get line points and weights
      double l_ptsLine1[2];
      double l_weightsLine1[2];
      edge::dg::QuadratureEval< QUAD4R, 2 >::line( l_ptsLine1,
                                                   l_weightsLine1 );

      // compute solution through numerical quadrature
      double l_quadUpdateL1[5][4][2];
      double l_quadUpdateR1[5][4][2];

      edge::elastic::solvers::InternalBoundary<QUAD4R, 5, 2, 2>::evalSpaceTime( l_faIdL, l_faIdR, l_vIdR,
                                                                                l_massInvDiag1,
                                                                                l_dT1, l_ptsLine1, l_weightsLine1,
                                                                                l_weightsFaces1, l_basisFaces1,
                                                                                l_trafoInv,
                                                                                l_msJump1,
                                                                                l_msFluxL1, l_msFluxR1,
                                                                                l_tDerL1, l_tDerR1,
                                                                                l_scatchMem1,
                                                                                l_quadUpdateL1, l_quadUpdateR1 );

      /*
       * Compare the two solutions
       */
      for( int_qt l_qt = 0; l_qt < 5; l_qt++ ) {
        for( int_md l_md = 0; l_md < 4; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < 2; l_ru++ ) {
            REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
            REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
          }
        }
      }

    } // remote face id
  } // local face id
}


TEST_CASE( "Eval of internal boundary conditions through quadrature rules, tria3", "[intBnd][tria3Eval]" ) {
  // face- and vertex-ids
  unsigned short l_faIdL, l_faIdR, l_vIdR;

  // normals and tangents
  double l_n[2];

  // lame paramters
  double l_lam[2];
  double l_mu[2];
  double l_rho[2];

  /*
   * Quadrature-free reference solution of first unit test:
   *   element type: tria3
   *   order: 5 -> 15 element modes
   *   quantities: 5
   *   concurrent forward runs: 8
   *   lame parameters left element:  2.9, 1.5, 1.4
   *   lame parameters right element: 3.3, 1.8, 1.2
   */

  // set vertex id to only possible value
  l_vIdR  = 0;

  // set up lame parameters and density
  l_lam[0] = 2.9; l_mu[0] = 1.5; l_rho[0] = 1.4;
  l_lam[1] = 3.3; l_mu[1] = 1.8; l_rho[1] = 1.2;

  // set up normal
  l_n[0] = -0.1;
  l_n[1] = -0.99498743710662;

  // set up shared basis
  edge::dg::Basis l_basis1( TRIA3, 5 );

  // set up quad-free flux matrices
  double l_fluxMatrices1[CE_N_FLUX_MATRICES(TRIA3)][15][15];
  l_basis1.getFluxMm1Dense( 15, (double *) l_fluxMatrices1, true );

  // set up quadrature-free flux solvers
  double l_fluxSolversL1[2][5][5];
  double l_fluxSolversR1[2][5][5];
  edge::elastic::solvers::common::setupSolver2d( l_rho[0], l_rho[1],
                                                 l_lam[0], l_lam[1],
                                                 l_mu[ 0], l_mu[ 1],
                                                 l_n[0], l_n[1], 0,
                                                 l_fluxSolversL1[0], l_fluxSolversL1[1] );

  edge::elastic::solvers::common::setupSolver2d(  l_rho[1], l_rho[0],
                                                  l_lam[1], l_lam[0],
                                                  l_mu[ 1],  l_mu[0],
                                                 -l_n[0], -l_n[1], 0,
                                                  l_fluxSolversR1[0], l_fluxSolversR1[1] );

  // trafo from physical to face-coordinates
  double l_trafoInv[5][5];
  edge::elastic::common::setupTrafoInv2d( l_n[0], l_n[1], l_trafoInv );

  // middle state solver
  double l_msJump1[5][5];
  edge::elastic::solvers::InternalBoundarySolvers<TRIA3>::solvMsJumpL( l_lam[0], l_lam[1],
                                                                       l_mu[0],  l_mu[1],
                                                                       l_rho[0], l_rho[1],
                                                                       l_msJump1 );

  // flux solvers
  double l_msFluxL1[5][5];
  double l_msFluxR1[5][5];
  edge::elastic::solvers::InternalBoundarySolvers<TRIA3>::solvMsFlux( l_n[0], l_n[1],
                                                                      l_lam[0], l_mu[0], l_rho[0],
                                                                      l_msFluxL1 );
  edge::elastic::solvers::InternalBoundarySolvers<TRIA3>::solvMsFlux( l_n[0], l_n[1],
                                                                      l_lam[1], l_mu[1], l_rho[1],
                                                                      l_msFluxR1 );

  // change sign of flux solvers, which don't have a build-in minus for intenral boundaries
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_msFluxL1[l_q1][l_q2] *= -1;
      l_msFluxR1[l_q1][l_q2] *= -1;
    }
  }

  // get inverse mass matrix
  double l_massInvDense1[15][15];
  l_basis1.getMassInvDense( 15, l_massInvDense1[0], true );
  double l_massInvDiag1[15];
  for( int_md l_md = 0; l_md < 15; l_md++ ) l_massInvDiag1[l_md] = l_massInvDense1[l_md][l_md];

  // get face-quad points, weights and evaluated basis
  double l_ptsFaces1[6][5][2];
  double l_weightsFaces1[5];
  double l_basisFaces1[6][5][15];
  edge::dg::QuadratureEval< TRIA3, 5 >::faces( l_ptsFaces1,
                                               l_weightsFaces1,
                                               l_basisFaces1 );

  // set up DOFs
  double l_dofsL1[5][15][8];
  double l_dofsR1[5][15][8];

  for( unsigned short l_q = 0; l_q < 5; l_q++ ) {
    for( unsigned short l_md = 0; l_md < 15; l_md++ ) {
      for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
        l_dofsL1[l_q][l_md][l_ru] = (-0.5 - l_q + l_md) * (l_ru+1);
        l_dofsR1[l_q][l_md][l_ru] = ( 0.5 + l_q - l_md) * (l_ru+1);
      }
    }
  }

  /*
   * Check for quadrature in space
   */
  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 3; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR = 0; l_faIdR < 3; l_faIdR++ ) {
      // compute reference solution
      double l_refTmp1[5][15][8];
      double l_refUpdateL1[5][15][8];
      double l_refUpdateR1[5][15][8];


      // left element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_dofsL1[0][0],
                                             l_fluxMatrices1[l_faIdL][0],
                                             l_refTmp1[0][0] );

      // left element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversL1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );

      // left element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_dofsR1[0][0],
                                             l_fluxMatrices1[3+l_faIdL*3+l_faIdR][0],
                                             l_refTmp1[0][0] );

      // left element, neighboring contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversL1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );


      // right element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_dofsR1[0][0],
                                             l_fluxMatrices1[l_faIdR][0],
                                             l_refTmp1[0][0] );

      // right element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversR1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      // right element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_dofsL1[0][0],
                                             l_fluxMatrices1[3+l_faIdR*3+l_faIdL][0],
                                             l_refTmp1[0][0] );

      // right element, neighboring contribution multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversR1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      /*
       * Solution through numerical quadrature
       */
      // compute solution through numerical quadrature
      double l_quadUpdateL1[5][15][8];
      double l_quadUpdateR1[5][15][8];

      edge::elastic::solvers::InternalBoundary<TRIA3, 5, 5, 8>::evalSpace( l_faIdL, l_faIdR, l_vIdR,
                                                                           l_massInvDiag1, l_weightsFaces1, l_basisFaces1,
                                                                           l_trafoInv,
                                                                           l_msJump1,
                                                                           l_msFluxL1, l_msFluxR1,
                                                                           l_dofsL1, l_dofsR1,
                                                                           l_quadUpdateL1, l_quadUpdateR1 );

      /*
       * Compare the two solutions
       */
      for( int_qt l_qt = 0; l_qt < 5; l_qt++ ) {
        for( int_md l_md = 0; l_md < 15; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
            REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
            REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
          }
        }
      }

    } // remote face id
  } // local face id

  /*
   * Now check for quadrature in space and time.
   */
  // time integrated dofs
  double l_tIntL1[5][15][8];
  double l_tIntR1[5][15][8];
  // time derivatives
  double l_tDerL1[5][5][15][8];
  double l_tDerR1[5][5][15][8];
  // time step
  double l_dT1 = 0.73;
  // scratch memory
  double l_scatchMem1[4][5][15][8];

  // set up star matrix
  double l_starL1[2][5][5];
  double l_starR1[2][5][5];
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
    for( unsigned short l_qt1 = 0; l_qt1 < 5; l_qt1++ ) {
      for( unsigned short l_qt2 = 0; l_qt2 < 5; l_qt2++ ) {
        l_starL1[l_di][l_qt1][l_qt2] = l_di - l_qt1 * l_qt2;
        l_starR1[l_di][l_qt1][l_qt2] = l_di * l_qt1 - l_qt2 - 0.5;
      }
    }
  }

  // set up quad-free transposed stiffness matrices (mutiplied with inverse mass matrix)
  double l_stiffT1[2][15][15];
  l_basis1.getStiffMm1Dense( 15, (double *) l_stiffT1, true,  true );

  // perform time predictions
  edge::elastic::solvers::TimePred< TRIA3, 5, 5, 8 >::ckVanilla( l_dT1,
                                                                 l_stiffT1,
                                                                 l_starL1,
                                                                 l_dofsL1,
                                                                 l_scatchMem1[0],
                                                                 l_tDerL1,
                                                                 l_tIntL1 );

  edge::elastic::solvers::TimePred< TRIA3, 5, 5, 8 >::ckVanilla( l_dT1,
                                                                 l_stiffT1,
                                                                 l_starR1,
                                                                 l_dofsR1,
                                                                 l_scatchMem1[0],
                                                                 l_tDerR1,
                                                                 l_tIntR1 );

  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 3; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR =0; l_faIdR < 3; l_faIdR++ ) {
      // compute reference solution
      double l_refTmp1[5][15][8];
      double l_refUpdateL1[5][15][8];
      double l_refUpdateR1[5][15][8];


      // left element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_tIntL1[0][0],
                                             l_fluxMatrices1[l_faIdL][0],
                                             l_refTmp1[0][0] );

      // left element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversL1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );

      // left element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_tIntR1[0][0],
                                             l_fluxMatrices1[3+l_faIdL*3+l_faIdR][0],
                                             l_refTmp1[0][0] );

      // left element, neighboring contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversL1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateL1[0][0] );


      // right element, local contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_tIntR1[0][0],
                                             l_fluxMatrices1[l_faIdR][0],
                                             l_refTmp1[0][0] );

      // right element, local contribution: multiply with flux solver
      edge::linalg::Matrix::matMulB0FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversR1[0][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      // right element, neighboring contribution: multiply with flux matrix
      edge::linalg::Matrix::matMulB0FusedAC( 8,
                                             5, 15, 15,
                                             l_tIntL1[0][0],
                                             l_fluxMatrices1[3+l_faIdR*3+l_faIdL][0],
                                             l_refTmp1[0][0] );

      // right element, neighboring contribution multiply with flux solver
      edge::linalg::Matrix::matMulB1FusedBC( 8,
                                             5, 15, 5,
                                             l_fluxSolversR1[1][0],
                                             l_refTmp1[0][0],
                                             l_refUpdateR1[0][0] );

      /*
       * Solution through numerical quadrature
       */
      // get line points and weights
      double l_ptsLine1[5];
      double l_weightsLine1[5];
      edge::dg::QuadratureEval< TRIA3, 5 >::line( l_ptsLine1,
                                                  l_weightsLine1 );

      // compute solution through numerical quadrature
      double l_quadUpdateL1[5][15][8];
      double l_quadUpdateR1[5][15][8];

      edge::elastic::solvers::InternalBoundary<TRIA3, 5, 5, 8>::evalSpaceTime( l_faIdL, l_faIdR, l_vIdR,
                                                                               l_massInvDiag1,
                                                                               l_dT1, l_ptsLine1, l_weightsLine1,
                                                                               l_weightsFaces1, l_basisFaces1,
                                                                               l_trafoInv,
                                                                               l_msJump1,
                                                                               l_msFluxL1, l_msFluxR1,
                                                                               l_tDerL1, l_tDerR1,
                                                                               l_scatchMem1,
                                                                               l_quadUpdateL1, l_quadUpdateR1 );

      /*
       * Compare the two solutions
       */
      for( int_qt l_qt = 0; l_qt < 5; l_qt++ ) {
        for( int_md l_md = 0; l_md < 15; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
            REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
            REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
          }
        }
      }

    } // remote face id
  } // local face id
}

TEST_CASE( "Eval of internal boundary conditions through quadrature rules, hex8r", "[intBnd][hex8rEval]" ) {
  // face- and vertex-ids
  unsigned short l_faIdL, l_faIdR;

  // normals and tangents
  double l_n[3];
  double l_s[3];
  double l_t[3];

  // lame paramters
  double l_lam[2];
  double l_mu[2];
  double l_rho[2];

  /*
   * Quadrature-free reference solution of first unit test:
   *   element type: hex8r
   *   order: 3 -> 27 element modes
   *   quantities: 9
   *   concurrent forward runs: 8
   *   lame parameters left element:  2.9, 1.5, 2.5
   *   lame parameters right element: 3.4, 1.2, 1.3
   */

  // set up lame parameters and density
  l_lam[0] = 2.9; l_mu[0] = 1.5; l_rho[0] = 2.5;
  l_lam[1] = 3.4; l_mu[1] = 1.2; l_rho[1] = 1.3;

  // set up basis (normal and tangents)
  l_n[0] = 0;
  l_n[1] = -1;
  l_n[2] = 0;

  l_s[0] = 0;
  l_s[1] = 0;
  l_s[2] = -1;

  l_t[0] = -1;
  l_t[1] = 0;
  l_t[2] = 0;

  // set up shared dg basis
  edge::dg::Basis l_basis1( HEX8R, 3 );

  // set up quad-free flux matrices
  double l_fluxMatrices1[CE_N_FLUX_MATRICES(HEX8R)][27][27];
  l_basis1.getFluxMm1Dense( 27, (double *) l_fluxMatrices1, true );
  // set up quadrature-free flux solvers
  double l_fluxSolversL1[2][9][9];
  double l_fluxSolversR1[2][9][9];
  edge::elastic::solvers::common::setupSolver3d( l_rho[0], l_rho[1],
                                                 l_lam[0], l_lam[1],
                                                 l_mu[ 0], l_mu[ 1],
                                                 l_n[0], l_n[1], l_n[2],
                                                 l_s[0], l_s[1], l_s[2],
                                                 l_t[0], l_t[1], l_t[2],
                                                 l_fluxSolversL1[0], l_fluxSolversL1[1] );

  edge::elastic::solvers::common::setupSolver3d( l_rho[1], l_rho[0],
                                                 l_lam[1], l_lam[0],
                                                 l_mu[ 1],  l_mu[0],
                                                 -l_n[0], -l_n[1], -l_n[2],
                                                 l_s[0], l_s[1], l_s[2],
                                                 l_t[0], l_t[1], l_t[2],
                                                 l_fluxSolversR1[0], l_fluxSolversR1[1] );

  // trafo from physical to face-coordinates
  double l_trafoInv[9][9];
  edge::elastic::common::setupTrafoInv3d( l_n[0], l_n[1], l_n[2],
                                          l_s[0], l_s[1], l_s[2],
                                          l_t[0], l_t[1], l_t[2],
                                          l_trafoInv );

  // middle state solver
  double l_msJump1[9][9];
  edge::elastic::solvers::InternalBoundarySolvers<HEX8R>::solvMsJumpL( l_lam[0], l_lam[1],
                                                                       l_mu[0],  l_mu[1],
                                                                       l_rho[0], l_rho[1],
                                                                       l_msJump1 );

  // flux solvers
  double l_msFluxL1[9][9];
  double l_msFluxR1[9][9];
  edge::elastic::solvers::InternalBoundarySolvers<HEX8R>::solvMsFlux( l_n[0], l_n[1], l_n[2],
                                                                      l_s[0], l_s[1], l_s[2],
                                                                      l_t[0], l_t[1], l_t[2],
                                                                      l_lam[0], l_mu[0], l_rho[0],
                                                                      l_msFluxL1 );
  edge::elastic::solvers::InternalBoundarySolvers<HEX8R>::solvMsFlux( l_n[0], l_n[1], l_n[2],
                                                                      l_s[0], l_s[1], l_s[2],
                                                                      l_t[0], l_t[1], l_t[2],
                                                                      l_lam[1], l_mu[1], l_rho[1],
                                                                      l_msFluxR1 );

  // change sign of flux solvers, which don't have a build-in minus for intenral boundaries
  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      l_msFluxL1[l_q1][l_q2] *= -1;
      l_msFluxR1[l_q1][l_q2] *= -1;
    }
  }

  // get inverse mass matrix
  double l_massInvDense1[27][27];
  l_basis1.getMassInvDense( 27, l_massInvDense1[0], true );
  double l_massInvDiag1[27];
  for( int_md l_md = 0; l_md < 27; l_md++ ) l_massInvDiag1[l_md] = l_massInvDense1[l_md][l_md];

  // get face-quad points, weights and evaluated basis
  double l_ptsFaces1[12][9][3];
  double l_weightsFaces1[12];
  double l_basisFaces1[12][9][27];
  edge::dg::QuadratureEval< HEX8R, 3 >::faces( l_ptsFaces1,
                                               l_weightsFaces1,
                                               l_basisFaces1 );

  // set up DOFs
  double l_dofsL1[9][27][8];
  double l_dofsR1[9][27][8];

  for( unsigned short l_q = 0; l_q < 9; l_q++ ) {
    for( unsigned short l_md = 0; l_md < 27; l_md++ ) {
      for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
        l_dofsL1[l_q][l_md][l_ru] = (-0.5 - l_q + l_md) - (l_ru+1.3);
        l_dofsR1[l_q][l_md][l_ru] = ( 0.5 + l_q - l_md) + (l_ru+1.2);
      }
    }
  }


  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 6; l_faIdL++ ) {
    // derive remote face id
    if(      l_faIdL == 0 ) l_faIdR = 5;
    else if( l_faIdL == 1 ) l_faIdR = 3;
    else if( l_faIdL == 2 ) l_faIdR = 4;
    else if( l_faIdL == 3 ) l_faIdR = 1;
    else if( l_faIdL == 4 ) l_faIdR = 2;
    else if( l_faIdL == 5 ) l_faIdR = 0;

    // compute reference solution
    double l_refTmp1[9][27][8];
    double l_refUpdateL1[9][27][8];
    double l_refUpdateR1[9][27][8];


    // left element, local contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_dofsL1[0][0],
                                           l_fluxMatrices1[l_faIdL][0],
                                           l_refTmp1[0][0] );

    // left element, local contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB0FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversL1[0][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateL1[0][0] );

    // derive neighboring contrib flux matrix id of left element
    unsigned short l_fIdL = 6 + l_faIdL;

    // left element, neighboring contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_dofsR1[0][0],
                                           l_fluxMatrices1[l_fIdL][0],
                                           l_refTmp1[0][0] );

    // left element, neighboring contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB1FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversL1[1][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateL1[0][0] );


    // right element, local contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_dofsR1[0][0],
                                           l_fluxMatrices1[l_faIdR][0],
                                           l_refTmp1[0][0] );

    // right element, local contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB0FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversR1[0][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateR1[0][0] );

    // derive neighboring contrib flux matrix id of right element
    unsigned short l_fIdR = 6 + l_faIdR;

    // right element, neighboring contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_dofsL1[0][0],
                                           l_fluxMatrices1[l_fIdR][0],
                                           l_refTmp1[0][0] );

    // right element, neighboring contribution multiply with flux solver
    edge::linalg::Matrix::matMulB1FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversR1[1][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateR1[0][0] );

    /*
     * Solution through numerical quadrature
     */
    // compute solution through numerical quadrature
    double l_quadUpdateL1[9][27][8];
    double l_quadUpdateR1[9][27][8];

    edge::elastic::solvers::InternalBoundary<HEX8R, 9, 3, 8>::evalSpace( l_faIdL, l_faIdR, 0,
                                                                         l_massInvDiag1, l_weightsFaces1, l_basisFaces1,
                                                                         l_trafoInv,
                                                                         l_msJump1,
                                                                         l_msFluxL1, l_msFluxR1,
                                                                         l_dofsL1, l_dofsR1,
                                                                         l_quadUpdateL1, l_quadUpdateR1 );

    /*
     * Compare the two solutions
     */
    for( int_qt l_qt = 0; l_qt < 9; l_qt++ ) {
      for( int_md l_md = 0; l_md < 27; l_md++ ) {
        for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
          REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
          REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
        }
      }
    }
  } // local face id

  /*
   * Now check for quadrature in space and time.
   */
  // time integrated dofs
  double l_tIntL1[9][27][8];
  double l_tIntR1[9][27][8];
  // time derivatives
  double l_tDerL1[3][9][27][8];
  double l_tDerR1[3][9][27][8];
  // time step
  double l_dT1 = 0.512;
  // scratch memory
  double l_scatchMem1[4][9][27][8];

  // set up star matrix
  double l_starL1[3][9][9];
  double l_starR1[3][9][9];
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    for( unsigned short l_qt1 = 0; l_qt1 < 9; l_qt1++ ) {
      for( unsigned short l_qt2 = 0; l_qt2 < 9; l_qt2++ ) {
        l_starL1[l_di][l_qt1][l_qt2] = l_di + l_qt1 - l_qt2;
        l_starR1[l_di][l_qt1][l_qt2] = l_di - l_qt1 + l_qt2 - 0.5;
      }
    }
  }

  // set up quad-free transposed stiffness matrices (mutiplied with inverse mass matrix)
  double l_stiffT1[3][27][27];
  l_basis1.getStiffMm1Dense( 27, (double *) l_stiffT1, true,  true );

  // perform time predictions
  edge::elastic::solvers::TimePred< HEX8R, 9, 3, 8 >::ckVanilla( l_dT1,
                                                                 l_stiffT1,
                                                                 l_starL1,
                                                                 l_dofsL1,
                                                                 l_scatchMem1[0],
                                                                 l_tDerL1,
                                                                 l_tIntL1 );

  edge::elastic::solvers::TimePred< HEX8R, 9, 3, 8 >::ckVanilla( l_dT1,
                                                                 l_stiffT1,
                                                                 l_starR1,
                                                                 l_dofsR1,
                                                                 l_scatchMem1[0],
                                                                 l_tDerR1,
                                                                 l_tIntR1 );

  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 6; l_faIdL++ ) {
    // derive remote face id
    if(      l_faIdL == 0 ) l_faIdR = 5;
    else if( l_faIdL == 1 ) l_faIdR = 3;
    else if( l_faIdL == 2 ) l_faIdR = 4;
    else if( l_faIdL == 3 ) l_faIdR = 1;
    else if( l_faIdL == 4 ) l_faIdR = 2;
    else if( l_faIdL == 5 ) l_faIdR = 0;

    // compute reference solution
    double l_refTmp1[9][27][8];
    double l_refUpdateL1[9][27][8];
    double l_refUpdateR1[9][27][8];


    // left element, local contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_tIntL1[0][0],
                                           l_fluxMatrices1[l_faIdL][0],
                                           l_refTmp1[0][0] );

    // left element, local contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB0FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversL1[0][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateL1[0][0] );

    // derive neighboring contrib flux matrix id of left element
    unsigned short l_fIdL = 6 + l_faIdL;

    // left element, neighboring contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_tIntR1[0][0],
                                           l_fluxMatrices1[l_fIdL][0],
                                           l_refTmp1[0][0] );

    // left element, neighboring contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB1FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversL1[1][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateL1[0][0] );


    // right element, local contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_tIntR1[0][0],
                                           l_fluxMatrices1[l_faIdR][0],
                                           l_refTmp1[0][0] );

    // right element, local contribution: multiply with flux solver
    edge::linalg::Matrix::matMulB0FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversR1[0][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateR1[0][0] );

    // derive neighboring contrib flux matrix id of left element
    unsigned short l_fIdR = 6 + l_faIdR;

    // right element, neighboring contribution: multiply with flux matrix
    edge::linalg::Matrix::matMulB0FusedAC( 8,
                                           9, 27, 27,
                                           l_tIntL1[0][0],
                                           l_fluxMatrices1[l_fIdR][0],
                                           l_refTmp1[0][0] );

    // right element, neighboring contribution multiply with flux solver
    edge::linalg::Matrix::matMulB1FusedBC( 8,
                                           9, 27, 9,
                                           l_fluxSolversR1[1][0],
                                           l_refTmp1[0][0],
                                           l_refUpdateR1[0][0] );

    /*
     * Solution through numerical quadrature
     */
    // get line points and weights
    double l_ptsLine1[3];
    double l_weightsLine1[3];
    edge::dg::QuadratureEval< HEX8R, 3 >::line( l_ptsLine1,
                                                l_weightsLine1 );

    // compute solution through numerical quadrature
    double l_quadUpdateL1[9][27][8];
    double l_quadUpdateR1[9][27][8];

    edge::elastic::solvers::InternalBoundary<HEX8R, 9, 3, 8>::evalSpaceTime( l_faIdL, l_faIdR, 0,
                                                                             l_massInvDiag1,
                                                                             l_dT1, l_ptsLine1, l_weightsLine1,
                                                                             l_weightsFaces1, l_basisFaces1,
                                                                             l_trafoInv,
                                                                             l_msJump1,
                                                                             l_msFluxL1, l_msFluxR1,
                                                                             l_tDerL1, l_tDerR1,
                                                                             l_scatchMem1,
                                                                             l_quadUpdateL1, l_quadUpdateR1 );

    /*
     * Compare the two solutions
     */
    for( int_qt l_qt = 0; l_qt < 9; l_qt++ ) {
      for( int_md l_md = 0; l_md < 27; l_md++ ) {
        for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
          REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
          REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
        }
      }
    }
  } // local face id
}

TEST_CASE( "Eval of internal boundary conditions through quadrature rules, tet4", "[intBnd][tet4Eval]" ) {
  // face- and vertex-ids
  unsigned short l_faIdL, l_faIdR, l_veIdR;
  // left veId (determine at runtime)
  unsigned short l_veIdL = std::numeric_limits<unsigned short>::max();

  // normals and tangents
  double l_n[3];
  double l_s[3];
  double l_t[3];

  // lame paramters
  double l_lam[2];
  double l_mu[2];
  double l_rho[2];

  /*
   * Quadrature-free reference solution of first unit test:
   *   element type: tet4
   *   order: 4 -> 20 element modes
   *   quantities: 9
   *   concurrent forward runs: 8
   *   lame parameters left element:  32E6, 29E6, 2670
   *   lame parameters right element: 33E6, 24E6, 2530
   */

  // set up lame parameters and density
  l_lam[0] = 32E6; l_mu[0] = 29E6; l_rho[0] = 2670;
  l_lam[1] = 33E6; l_mu[1] = 34E6; l_rho[1] = 2530;

  // set up orthormal basis (normal and tangents)
  // WolframAlpha:
  //   In: gram schmidt {{1.2,1,1.9},{2.4,2.1,0},{5.6,1.9,3.1}}
  //   Out: Orthonormal basis:
  //        (0.487869, 0.406558, 0.77246) | (0.573171, 0.518229, -0.634755) | (0.658376, -0.752429, -0.0198008)
  l_n[0] = 0.487869;
  l_n[1] = 0.406558;
  l_n[2] = 0.77246;

  l_s[0] = 0.573171;
  l_s[1] = 0.518229;
  l_s[2] = -0.634755;

  l_t[0] = 0.658376;
  l_t[1] = -0.752429;
  l_t[2] = -0.0198008;

  // set up shared basis
  edge::dg::Basis l_basis1( TET4, 4 );

  // set up quad-free flux matrices
  double l_fluxMatrices1[CE_N_FLUX_MATRICES(TET4)][20][20];
  l_basis1.getFluxMm1Dense( 20, (double *) l_fluxMatrices1, true );

  // set up quadrature-free flux solvers
  double l_fluxSolversL1[2][9][9];
  double l_fluxSolversR1[2][9][9];
  edge::elastic::solvers::common::setupSolver3d( l_rho[0], l_rho[1],
                                                 l_lam[0], l_lam[1],
                                                 l_mu[ 0], l_mu[ 1],
                                                 l_n[0], l_n[1], l_n[2],
                                                 l_s[0], l_s[1], l_s[2],
                                                 l_t[0], l_t[1], l_t[2],
                                                 l_fluxSolversL1[0], l_fluxSolversL1[1] );

  edge::elastic::solvers::common::setupSolver3d( l_rho[1], l_rho[0],
                                                 l_lam[1], l_lam[0],
                                                 l_mu[ 1],  l_mu[0],
                                                 -l_n[0], -l_n[1], -l_n[2],
                                                 l_s[0], l_s[1], l_s[2],
                                                 l_t[0], l_t[1], l_t[2],
                                                 l_fluxSolversR1[0], l_fluxSolversR1[1] );

  // trafo from physical to face-coordinates
  double l_trafoInv[9][9];
  edge::elastic::common::setupTrafoInv3d( l_n[0], l_n[1], l_n[2],
                                          l_s[0], l_s[1], l_s[2],
                                          l_t[0], l_t[1], l_t[2],
                                          l_trafoInv );

  // middle state solver
  double l_msJump1[9][9];
  edge::elastic::solvers::InternalBoundarySolvers<TET4>::solvMsJumpL( l_lam[0], l_lam[1],
                                                                      l_mu[0],  l_mu[1],
                                                                      l_rho[0], l_rho[1],
                                                                      l_msJump1 );

  // flux solvers
  double l_msFluxL1[9][9];
  double l_msFluxR1[9][9];
  edge::elastic::solvers::InternalBoundarySolvers<TET4>::solvMsFlux( l_n[0], l_n[1], l_n[2],
                                                                     l_s[0], l_s[1], l_s[2],
                                                                     l_t[0], l_t[1], l_t[2],
                                                                     l_lam[0], l_mu[0], l_rho[0],
                                                                     l_msFluxL1 );
  edge::elastic::solvers::InternalBoundarySolvers<TET4>::solvMsFlux( l_n[0], l_n[1], l_n[2],
                                                                     l_s[0], l_s[1], l_s[2],
                                                                     l_t[0], l_t[1], l_t[2],
                                                                     l_lam[1], l_mu[1], l_rho[1],
                                                                     l_msFluxR1 );

  // change sign of flux solvers, which don't have a build-in minus for intenral boundaries
  for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
      l_msFluxL1[l_q1][l_q2] *= -1;
      l_msFluxR1[l_q1][l_q2] *= -1;
    }
  }

  // get inverse mass matrix
  double l_massInvDense1[20][20];
  l_basis1.getMassInvDense( 20, l_massInvDense1[0], true );
  double l_massInvDiag1[20];
  for( int_md l_md = 0; l_md < 20; l_md++ ) l_massInvDiag1[l_md] = l_massInvDense1[l_md][l_md];

  // get face-quad points, weights and evaluated basis
  double l_ptsFaces1[16][16][3];
  double l_weightsFaces1[16];
  double l_basisFaces1[16][16][20];
  edge::dg::QuadratureEval< TET4, 4 >::faces( l_ptsFaces1,
                                              l_weightsFaces1,
                                              l_basisFaces1 );

  // set up DOFs
  double l_dofsL1[9][20][8];
  double l_dofsR1[9][20][8];

  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
        l_dofsL1[l_qt][l_md][l_ru] = (-0.5 - l_qt + 0.2*l_md) + (l_ru+1.4);
        l_dofsR1[l_qt][l_md][l_ru] = ( 0.5 + l_qt - 0.3*l_md) - (l_ru-1.9);
      }
    }
  }

  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 4; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR = 0; l_faIdR < 4; l_faIdR++ ) {
      // iterate over face's vertex combinations
      for( l_veIdR = 0; l_veIdR < 3; l_veIdR++ ) {
        // compute reference solution
        double l_refTmp1[9][20][8];
        double l_refUpdateL1[9][20][8];
        double l_refUpdateR1[9][20][8];


        // left element, local contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_dofsL1[0][0],
                                               l_fluxMatrices1[l_faIdL][0],
                                               l_refTmp1[0][0] );

        // left element, local contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB0FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversL1[0][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateL1[0][0] );

        // derive neighboring contrib flux matrix id of left element
        // jump over local flux matrices
        unsigned short l_fIdL = 4;
        // jump over local face
        l_fIdL += l_faIdL * 4 * 3;

        // jump over neighboring face
        l_fIdL += l_faIdR * 3;

        // jump over vertices
        l_fIdL += l_veIdR;

        // left element, neighboring contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_dofsR1[0][0],
                                               l_fluxMatrices1[l_fIdL][0],
                                               l_refTmp1[0][0] );

        // left element, neighboring contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB1FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversL1[1][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateL1[0][0] );


        // right element, local contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_dofsR1[0][0],
                                               l_fluxMatrices1[l_faIdR][0],
                                               l_refTmp1[0][0] );

        // right element, local contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB0FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversR1[0][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateR1[0][0] );

        // derive neighboring contrib flux matrix id of right element
        // jump over local flux matrices
        unsigned short l_fIdR = 4;
        // jump over local face
        l_fIdR += l_faIdR * 4 * 3;

        // jump over neighboring face
        l_fIdR += l_faIdL * 3;

        // left vertex id is identical to right (counter-clockwise -> clockwise)
        l_veIdL = l_veIdR;

        // jump over vertices
        l_fIdR += l_veIdL;

        // right element, neighboring contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_dofsL1[0][0],
                                               l_fluxMatrices1[l_fIdR][0],
                                               l_refTmp1[0][0] );

        // right element, neighboring contribution multiply with flux solver
        edge::linalg::Matrix::matMulB1FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversR1[1][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateR1[0][0] );

        /*
         * Solution through numerical quadrature
         */
        // compute solution through numerical quadrature
        double l_quadUpdateL1[9][20][8];
        double l_quadUpdateR1[9][20][8];

        edge::elastic::solvers::InternalBoundary<TET4, 9, 4, 8>::evalSpace( l_faIdL, l_faIdR, l_veIdR,
                                                                            l_massInvDiag1, l_weightsFaces1, l_basisFaces1,
                                                                            l_trafoInv,
                                                                            l_msJump1,
                                                                            l_msFluxL1, l_msFluxR1,
                                                                            l_dofsL1, l_dofsR1,
                                                                            l_quadUpdateL1, l_quadUpdateR1 );

        /*
         * Compare the two solutions
         */
        for( int_qt l_qt = 0; l_qt < 9; l_qt++ ) {
          for( int_md l_md = 0; l_md < 20; l_md++ ) {
            for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
              REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
              REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
            }
          }
        }

      } // vertex id
    } // remote face id
  } // local face id

  /*
   * Now check for quadrature in space and time.
   */
  // time integrated dofs
  double l_tIntL1[9][20][8];
  double l_tIntR1[9][20][8];
  // time derivatives
  double l_tDerL1[4][9][20][8];
  double l_tDerR1[4][9][20][8];
  // time step
  double l_dT1 = 0.1;
  // scratch memory
  double l_scatchMem1[4][9][20][8];

  // set up star matrix
  double l_starL1[3][9][9];
  double l_starR1[3][9][9];
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    for( unsigned short l_qt1 = 0; l_qt1 < 9; l_qt1++ ) {
      for( unsigned short l_qt2 = 0; l_qt2 < 9; l_qt2++ ) {
        l_starL1[l_di][l_qt1][l_qt2] = l_di + l_qt1 - l_qt2 + 0.25;
        l_starR1[l_di][l_qt1][l_qt2] = l_di - l_qt1 + l_qt2 - 0.5;
      }
    }
  }

  // set up quad-free transposed stiffness matrices (mutiplied with inverse mass matrix)
  double l_stiffT1[3][20][20];
  l_basis1.getStiffMm1Dense( 20, (double *) l_stiffT1, true,  true );

  // perform time predictions
  edge::elastic::solvers::TimePred< TET4, 9, 4, 8 >::ckVanilla( l_dT1,
                                                                l_stiffT1,
                                                                l_starL1,
                                                                l_dofsL1,
                                                                l_scatchMem1[0],
                                                                l_tDerL1,
                                                                l_tIntL1 );

  edge::elastic::solvers::TimePred< TET4, 9, 4, 8 >::ckVanilla( l_dT1,
                                                                l_stiffT1,
                                                                l_starR1,
                                                                l_dofsR1,
                                                                l_scatchMem1[0],
                                                                l_tDerR1,
                                                                l_tIntR1 );

  // iterate over local face ids
  for( l_faIdL = 0; l_faIdL < 4; l_faIdL++ ) {
    // iterate over remote face ids
    for( l_faIdR = 0; l_faIdR < 4; l_faIdR++ ) {
      // iterate over face's vertex combinations
      for( l_veIdR = 0; l_veIdR < 3; l_veIdR++ ) {
        // compute reference solution
        double l_refTmp1[9][20][8];
        double l_refUpdateL1[9][20][8];
        double l_refUpdateR1[9][20][8];


        // left element, local contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_tIntL1[0][0],
                                               l_fluxMatrices1[l_faIdL][0],
                                               l_refTmp1[0][0] );

        // left element, local contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB0FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversL1[0][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateL1[0][0] );

        // derive neighboring contrib flux matrix id of left element
        // jump over local flux matrices
        unsigned short l_fIdL = 4;
        // jump over local face
        l_fIdL += l_faIdL * 4 * 3;

        // jump over neighboring face
        l_fIdL += l_faIdR * 3;

        // jump over vertices
        l_fIdL += l_veIdR;

        // left element, neighboring contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_tIntR1[0][0],
                                               l_fluxMatrices1[l_fIdL][0],
                                               l_refTmp1[0][0] );

        // left element, neighboring contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB1FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversL1[1][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateL1[0][0] );


        // right element, local contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_tIntR1[0][0],
                                               l_fluxMatrices1[l_faIdR][0],
                                               l_refTmp1[0][0] );

        // right element, local contribution: multiply with flux solver
        edge::linalg::Matrix::matMulB0FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversR1[0][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateR1[0][0] );


        // derive neighboring contrib flux matrix id of right element
        // jump over local flux matrices
        unsigned short l_fIdR = 4;
        // jump over local face
        l_fIdR += l_faIdR * 4 * 3;

        // jump over neighboring face
        l_fIdR += l_faIdL * 3;

        // left vertex id is identical to right (counter-clockwise -> clockwise)
        l_veIdL = l_veIdR;

        // jump over vertices
        l_fIdR += l_veIdL;

        // right element, neighboring contribution: multiply with flux matrix
        edge::linalg::Matrix::matMulB0FusedAC( 8,
                                               9, 20, 20,
                                               l_tIntL1[0][0],
                                               l_fluxMatrices1[l_fIdR][0],
                                               l_refTmp1[0][0] );

        // right element, neighboring contribution multiply with flux solver
        edge::linalg::Matrix::matMulB1FusedBC( 8,
                                               9, 20, 9,
                                               l_fluxSolversR1[1][0],
                                               l_refTmp1[0][0],
                                               l_refUpdateR1[0][0] );

        /*
         * Solution through numerical quadrature
         */
        // get line points and weights
        double l_ptsLine1[4];
        double l_weightsLine1[4];
        edge::dg::QuadratureEval< TET4, 4 >::line( l_ptsLine1,
                                                   l_weightsLine1 );

        // compute solution through numerical quadrature
        double l_quadUpdateL1[9][20][8];
        double l_quadUpdateR1[9][20][8];

        edge::elastic::solvers::InternalBoundary<TET4, 9, 4, 8 >::evalSpaceTime( l_faIdL, l_faIdR, l_veIdR,
                                                                                 l_massInvDiag1,
                                                                                 l_dT1, l_ptsLine1, l_weightsLine1,
                                                                                 l_weightsFaces1, l_basisFaces1,
                                                                                 l_trafoInv,
                                                                                 l_msJump1,
                                                                                 l_msFluxL1, l_msFluxR1,
                                                                                 l_tDerL1, l_tDerR1,
                                                                                 l_scatchMem1,
                                                                                 l_quadUpdateL1, l_quadUpdateR1 );

        /*
         * Compare the two solutions
         */
        for( int_qt l_qt = 0; l_qt < 9; l_qt++ ) {
          for( int_md l_md = 0; l_md < 20; l_md++ ) {
            for( int_cfr l_ru = 0; l_ru < 8; l_ru++ ) {
              REQUIRE( l_refUpdateL1[l_qt][l_md][l_ru] == Approx( l_quadUpdateL1[l_qt][l_md][l_ru] ) );
              REQUIRE( l_refUpdateR1[l_qt][l_md][l_ru] == Approx( l_quadUpdateR1[l_qt][l_md][l_ru] ) );
            }
          }
        }

      } // vertex id
    } // remote face id
  } // local face id

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
  edge::elastic::solvers::InternalBoundarySolvers<
    TRIA3 >::solvMsJumpL( l_lam, l_lam,
                          l_mu,  l_mu,
                          l_rho, l_rho,
                          l_msJ );

  edge::elastic::solvers::InternalBoundarySolvers<
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
  double l_massI[1] = {1.0};
  double l_dt = 0.004;
  double l_weightsFaces[1] = {1.0};
  double l_basisFaces[6][1][1] = { {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}} };
  double l_tm1[5][5] = { {1,0,0,0,0},
                         {0,1,0,0,0},
                         {0,0,1,0,0},
                         {0,0,0,1,0},
                         {0,0,0,0,1} };
  double l_dofsL[5][1][1] = { {{ 1.0E2}}, {{-2.0E3}}, {{-4.0E5}}, {{10.0}}, {{-15.0}} };
  double l_dofsR[5][1][1] = { {{-0.1E2}}, {{-0.9E3}}, {{ 5.0E5}}, {{ 4.0}}, {{  5.0}} };

  edge::elastic::solvers::t_LinSlipWeakGlobal<double,1>            l_gl;
  edge::elastic::solvers::t_LinSlipWeakFace<double>                l_fa;
  edge::elastic::solvers::t_LinSlipWeakFaceQuadPoint<double, 2, 1> l_qp[1];

  // assign the values of our example
  l_gl.mus[0]   = 0.677;
  l_gl.mud[0]   = 0.55;
  l_gl.dcInv[0] = 2.5;

  l_fa.lEqM = false;
  l_fa.csDmuM = l_cs / l_mu;
  l_fa.csDmuP = l_cs / l_mu;

  l_qp[0].sn0[0]    = -120E6;
  l_qp[0].ss0[0][0] = 81.6E6;
  l_qp[0].muf[0]    = 0.59;
  l_qp[0].dd[0][0]  = 0.2;

  // generate a struct out out
  struct {
    edge::elastic::solvers::t_LinSlipWeakGlobal<double,1>             *gl;
    edge::elastic::solvers::t_LinSlipWeakFace<double>                 *fa;
    edge::elastic::solvers::t_LinSlipWeakFaceQuadPoint<double, 2, 1> (*qp)[1];
  } l_faData;
  l_faData.gl = &l_gl;
  l_faData.fa = &l_fa;
  l_faData.qp = &l_qp;

  // call the solver
  double l_upL[5][1][1];
  double l_upR[5][1][1];

  edge::elastic::solvers::InternalBoundary<
    TRIA3,
    5,
    1,
    1 >::evalSpace<
      double,
       edge::elastic::solvers::FrictionLaws< 2, 1 >
    >(  0, 0, 0, // irrelevant for FV
        l_massI,
        l_weightsFaces,
        l_basisFaces,
        l_tm1,
        l_msJ,
        l_msF,
        l_msF,
        l_dofsL,
        l_dofsR,
        l_upL,
        l_upR,
        1.0,
       &l_faData );

  // check the updates
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    REQUIRE( l_upL[l_qt][0][0] == Approx(-l_fluxLref[l_qt]) );
    REQUIRE( l_upR[l_qt][0][0] == Approx( l_fluxRref[l_qt]) );
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
  edge::elastic::solvers::InternalBoundarySolvers<
    TET4
  >::solvMsJumpL( l_lamL, l_lamR,
                  l_muL,  l_muR,
                  l_rhoL, l_rhoR,
                  l_msJ );

  edge::elastic::solvers::InternalBoundarySolvers<
    TET4
   >::solvMsFlux( 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0,
                  l_lamL, l_muL, l_rhoL,
                  l_msFL );
  edge::elastic::solvers::InternalBoundarySolvers<
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
  double l_massI[1] = {1.0};
  double l_dt = 0.004;
  double l_weightsFaces[1] = {1.0};
  double l_basisFaces[16][1][1] = { {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}},
                                    {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}},
                                    {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}},
                                    {{1.0}}, {{1.0}}, {{1.0}}, {{1.0}} };
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

  edge::elastic::solvers::t_LinSlipWeakGlobal<double,1>            l_gl;
  edge::elastic::solvers::t_LinSlipWeakFace<double>                l_fa;
  edge::elastic::solvers::t_LinSlipWeakFaceQuadPoint<double, 3, 1> l_qp[1];


  // assign the values of our example
  l_gl.mus[0]    = 0.677;
  l_gl.mud[0]    = 0.55;
  l_gl.dcInv[0]  = 2.5;

  l_fa.lEqM = false;
  l_fa.csDmuP = l_csL / l_muL;
  l_fa.csDmuM = l_csR / l_muR;

  l_qp[0].sn0[0]    = -120E6;
  l_qp[0].ss0[0][0] = 40E6;
  l_qp[0].ss0[1][0] = std::sqrt( 81.6E6 * 81.6E6 - 40.0E6 * 40.0E6 );
  l_qp[0].ss0A[0] = std::sqrt(   l_qp[0].ss0[0][0] * l_qp[0].ss0[0][0]
                               + l_qp[0].ss0[1][0] * l_qp[0].ss0[1][0] );
  l_qp[0].muf[0]    = 0.59;
  l_qp[0].dd[0][0]  = 0.2;
  l_qp[0].dd[1][0]  = -0.1;

  // generate a struct out of that
  struct {
    edge::elastic::solvers::t_LinSlipWeakGlobal<double,1>             *gl;
    edge::elastic::solvers::t_LinSlipWeakFace<double>                 *fa;
    edge::elastic::solvers::t_LinSlipWeakFaceQuadPoint<double, 3, 1> (*qp)[1];
  } l_faData;
  l_faData.gl = &l_gl;
  l_faData.fa = &l_fa;
  l_faData.qp = &l_qp;

  // call the solver
  double l_upL[9][1][1];
  double l_upR[9][1][1];

  edge::elastic::solvers::InternalBoundary<
    TET4,
    9,
    1,
    1 >::evalSpace<
      double,
      edge::elastic::solvers::FrictionLaws< 3, 1 >
    >(  0, 0, 0, // irrelevant for FV
        l_massI,
        l_weightsFaces,
        l_basisFaces,
        l_tm1,
        l_msJ,
        l_msFL,
        l_msFR,
        l_dofsL,
        l_dofsR,
        l_upL,
        l_upR,
        1.0,
       &l_faData );

  // check the updates
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    REQUIRE( l_upL[l_qt][0][0] == Approx(-l_fluxLref[l_qt]) );
    REQUIRE( l_upR[l_qt][0][0] == Approx( l_fluxRref[l_qt]) );
  }
}

#endif

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
 * Unit tests for evals through quadrature rules.
 **/

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64

#include <catch.hpp>
#define private public
#include "QuadratureEval.hpp"
#undef private

TEST_CASE( "Derivation of quad points and weights", "[quadratureEval][pts]" ) {
  double l_sum = 0;

  /*
   * quad4r, order2
   */
  real_mesh l_quad4rPts1[8][2][2];
  real_mesh l_quad4rWeights1[4];
  real_mesh l_quad4rBasisEval1[8][2][4];
  edge::dg::QuadratureEval<QUAD4R, 2>::faces( l_quad4rPts1, l_quad4rWeights1, l_quad4rBasisEval1 );
  l_sum = 0;
  for( unsigned short l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TET4, 2 ); l_qp++ ) {
    l_sum += l_quad4rWeights1[l_qp];
  }
  REQUIRE( l_sum == Approx(1.0) );

  // face 0 & 2
  REQUIRE( l_quad4rPts1[0][0][1] == Approx(0) );
  REQUIRE( l_quad4rPts1[0][1][1] == Approx(0) );

  REQUIRE( l_quad4rPts1[2][0][1] == Approx(1) );
  REQUIRE( l_quad4rPts1[2][1][1] == Approx(1) );

  REQUIRE( l_quad4rPts1[0][0][0] == Approx( l_quad4rPts1[2][1][0] ) );
  REQUIRE( l_quad4rPts1[0][1][0] == Approx( l_quad4rPts1[2][0][0] ) );

  // face 1 & 3
  REQUIRE( l_quad4rPts1[1][0][0] == Approx(1) );
  REQUIRE( l_quad4rPts1[1][1][0] == Approx(1) );

  REQUIRE( l_quad4rPts1[3][0][0] == Approx(0) );
  REQUIRE( l_quad4rPts1[3][1][0] == Approx(0) );

  REQUIRE( l_quad4rPts1[1][0][1] == Approx( l_quad4rPts1[3][1][1] ) );
  REQUIRE( l_quad4rPts1[1][1][1] == Approx( l_quad4rPts1[3][0][1] ) );

  /*
   * hex8r, order 2
   */
  real_mesh l_hex8rPts1[12][4][3];
  real_mesh l_hex8rWeights1[4];
  real_mesh l_hex8rBasisEval1[12][4][8];
  edge::dg::QuadratureEval<HEX8R, 2>::faces( l_hex8rPts1, l_hex8rWeights1, l_hex8rBasisEval1 );
  l_sum = 0;
  for( unsigned short l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TET4, 2 ); l_qp++ ) {
    l_sum += l_hex8rWeights1[l_qp];
  }
  REQUIRE( l_sum == Approx(1.0) );

  // face 0 & 5
  for( unsigned short l_qp = 0; l_qp < 4; l_qp++ ) {
    REQUIRE( l_hex8rPts1[ 0][l_qp][0] == Approx(l_hex8rPts1[11][l_qp][0]) );
    REQUIRE( l_hex8rPts1[ 0][l_qp][1] == Approx(l_hex8rPts1[11][l_qp][1]) );
    REQUIRE( l_hex8rPts1[ 0][l_qp][2] == Approx(0)                        );
    REQUIRE( l_hex8rPts1[11][l_qp][2] == Approx(1)                        );
  }

  // face 1 & 3
  for( unsigned short l_qp = 0; l_qp < 4; l_qp++ ) {
    REQUIRE( l_hex8rPts1[ 1][l_qp][0] == Approx(l_hex8rPts1[ 9][l_qp][0]) );
    REQUIRE( l_hex8rPts1[ 1][l_qp][2] == Approx(l_hex8rPts1[ 9][l_qp][2]) );
    REQUIRE( l_hex8rPts1[ 1][l_qp][1] == Approx(0)                        );
    REQUIRE( l_hex8rPts1[ 9][l_qp][1] == Approx(1)                        );
  }

  // face 2 & 4
  for( unsigned short l_qp = 0; l_qp < 4; l_qp++ ) {
    REQUIRE( l_hex8rPts1[ 2][l_qp][1] == Approx(l_hex8rPts1[10][l_qp][1]) );
    REQUIRE( l_hex8rPts1[ 2][l_qp][2] == Approx(l_hex8rPts1[10][l_qp][2]) );
    REQUIRE( l_hex8rPts1[ 2][l_qp][0] == Approx(1)                        );
    REQUIRE( l_hex8rPts1[10][l_qp][0] == Approx(0)                        );
  }

  /*
   * tet4, order 2
   */
  real_mesh l_tet4Pts1[16][4][3];
  real_mesh l_tet4Weights1[4];
  real_mesh l_tet4BasisEval1[16][4][4];

  edge::dg::QuadratureEval<TET4, 2>::faces( l_tet4Pts1, l_tet4Weights1, l_tet4BasisEval1 );
  // check that the sum of the quad points gives the surface of the "reference face"
  l_sum = 0;
  for( unsigned short l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TET4, 2 ); l_qp++ ) {
    l_sum += l_tet4Weights1[l_qp];
  }
  REQUIRE( l_sum == Approx(0.5) );

  /*
   * tet4, order 6
   */
  real_mesh l_tet4Pts2[16][36][3];
  real_mesh l_tet4Weights2[36];
  real_mesh l_tet4BasisEval2[16][36][56];

  edge::dg::QuadratureEval<TET4, 6>::faces( l_tet4Pts2, l_tet4Weights2, l_tet4BasisEval2 );
  // check that the sum of the quad points gives the surface of the "reference face"
  l_sum = 0;
  for( unsigned short l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TET4, 6 ); l_qp++ ) {
    l_sum += l_tet4Weights2[l_qp];
  }
  REQUIRE( l_sum == Approx(0.5) );
}

#endif

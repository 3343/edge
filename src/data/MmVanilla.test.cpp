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
 * Unit tests for the vanilla matrix kernels.
 **/
#include <catch.hpp>
#include "MmVanilla.hpp"

TEST_CASE( "Vanilla: GEMM.", "[mmVanilla][gemm]" ) {
  // vanilla kernels
  edge::data::MmVanilla< float > l_van1;

  // add GEMM kernel
  l_van1.add( 3, 3, 3,
              1, 1, 1,
              0, 0, 1 );

  l_van1.add( 3, 3, 3,
              1, 1, 1,
              0, 1, 1 );

  float l_mat1[3][3] = { { 1, 2, 3 },
                         { 4, 5, 6 },
                         { 7, 8, 9 } };

  float l_res1[3][3];

  // run first kernel
  l_van1.m_kernels[0]( l_mat1[0], l_mat1[0], l_res1[0] );

  // check result
  REQUIRE( l_res1[0][0] == Approx( 30 ) );
  REQUIRE( l_res1[0][1] == Approx( 36 ) );
  REQUIRE( l_res1[0][2] == Approx( 42 ) );

  REQUIRE( l_res1[1][0] == Approx( 66 ) );
  REQUIRE( l_res1[1][1] == Approx( 81 ) );
  REQUIRE( l_res1[1][2] == Approx( 96 ) );

  REQUIRE( l_res1[2][0] == Approx( 102 ) );
  REQUIRE( l_res1[2][1] == Approx( 126 ) );
  REQUIRE( l_res1[2][2] == Approx( 150 ) );


  // run second kernel (adds results)
  l_van1.m_kernels[1]( l_mat1[0], l_mat1[0], l_res1[0] );

  // check result
  REQUIRE( l_res1[0][0] == Approx( 2*30 ) );
  REQUIRE( l_res1[0][1] == Approx( 2*36 ) );
  REQUIRE( l_res1[0][2] == Approx( 2*42 ) );

  REQUIRE( l_res1[1][0] == Approx( 2*66 ) );
  REQUIRE( l_res1[1][1] == Approx( 2*81 ) );
  REQUIRE( l_res1[1][2] == Approx( 2*96 ) );

  REQUIRE( l_res1[2][0] == Approx( 2*102 ) );
  REQUIRE( l_res1[2][1] == Approx( 2*126 ) );
  REQUIRE( l_res1[2][2] == Approx( 2*150 ) );
}

/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Tests the f-wave solver.
 **/
#include <catch.hpp>
#define private public
#include "Fwave.hpp"
#undef private

TEST_CASE( "Test the derivation of eigenvalues.", "[Fwave][eigenValue]" ) {
   /*
    * Test case:
    *  h: 10 | 9
    *  u: -3 | 3
    *
    * char speeds: -3 - sqrt(10) = -6.16227766016837933 | 3 + 3 = 6
    * roe height: 9.5
    * roe velocity: (sqrt(10) * -3 + 3 * 3) / ( sqrt(10) + sqrt(9) ) = -0.0790021169691720
    * roe speeds: -0.079002116969172024 - sqrt(9.5) =- 3.1612091184536602 | -0.079002116969172024 + sqrt(9.5) = 3.00320488451531620
    */
  double l_hL[4]   = { 10,  10,  10,  10 };
  double l_hR[4]   = {  9,   9,   9,   9 };
  double l_uL[4]   = { -3,  -3,  -3,  -3 };
  double l_uR[4]   = {  3,   3,   3,   3 };
  double l_lam[2][4] = { { -1,  -1,  -1,  -1 },
                         { -1,  -1,  -1,  -1 } };

  edge::swe::solvers::Fwave< 4 >::evs( l_hL,   l_hR,
                                       l_uL,   l_uR,
                                       l_lam,
                                       1.0 );

  for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
    for( unsigned short l_cr = 0; l_cr < 4; l_cr++ ) {
      REQUIRE( l_lam[l_sd][l_cr] == Approx( -6.16227766016837933 ) );
      REQUIRE( l_lam[l_sd][l_cr] == Approx(  6 )                   );
    }
  }
}

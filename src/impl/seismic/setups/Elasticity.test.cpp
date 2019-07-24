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
 * Unit tests for the elastic setup.
 **/
#include <catch.hpp>
#define private public
#include "Elasticity.h"
#undef private

TEST_CASE( "Dense derivation of the elastic part of the star matrix in 2D", "[Elasticity][star2dDense]" ) {
#include "Elasticity.test.inc"

  // set up matrix structures
  double l_starE[2][5*5];
  double l_jacInv[2][2] = { {0.4, 0.3}, {0.5, 0.1} };

  // compute star matrices
  edge::seismic::setups::Elasticity::star( 2600,
                                           20.8E9,
                                           9.9E9,
                                           l_jacInv,
                                           l_starE );

  // check results
  for( unsigned short l_di = 0; l_di < 2; l_di++ )
    for( unsigned short l_m0 = 0; l_m0 < 5; l_m0++ )
      for( unsigned short l_m1 = 0; l_m1 < 5; l_m1++ )
        REQUIRE( l_starE[l_di][l_m0*5 + l_m1] == Approx(l_starE2d[l_di][l_m0][l_m1]) );
}

TEST_CASE( "Sparse derivation of the elastic part of the star matrix in 2D", "[Elasticity][star2dSparse]" ) {
#include "Elasticity.test.inc"

  // set up matrix structures
  double l_starE[2][10];
  double l_jacInv[2][2] = { {0.4, 0.3}, {0.5, 0.1} };

  // compute star matrices
  edge::seismic::setups::Elasticity::star( 2600,
                                           20.8E9,
                                           9.9E9,
                                           l_jacInv,
                                           l_starE );

  // check results
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
    unsigned short l_nz = 0;
    for( unsigned short l_m0 = 0; l_m0 < 5; l_m0++ ) {
      for( unsigned short l_m1 = 0; l_m1 < 5; l_m1++ ) {
        if( l_starE2d[l_di][l_m0][l_m1] != 0 ) {
          REQUIRE( l_starE[l_di][l_nz] == Approx(l_starE2d[l_di][l_m0][l_m1]) );
          l_nz++;
        }
      }
    }
    REQUIRE( l_nz == 10 );
  }
}

TEST_CASE( "Derivation of the elastic part of the star matrix in 3D", "[Elasticity][star3d]" ) {
#include "Elasticity.test.inc"

  // set up matrix structures
  double l_starE[3][9*9];
  double l_jacInv[3][3] = { {0.4, 0.3, 0.5}, {0.5, 0.1, 0.3}, {0.2, 0.4, 0.2} };

  // compute star matrices
  edge::seismic::setups::Elasticity::star( 2600,
                                           20.8E9,
                                           9.9E9,
                                           l_jacInv,
                                           l_starE );

  // check results
  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    for( unsigned short l_m0 = 0; l_m0 < 9; l_m0++ )
      for( unsigned short l_m1 = 0; l_m1 < 9; l_m1++ )
        REQUIRE( l_starE[l_di][l_m0*9 + l_m1] == Approx(l_starE3d[l_di][l_m0][l_m1]) );
}

TEST_CASE( "Sparse derivation of the elastic part of the star matrix in 3D", "[Elasticity][star3dSparse]" ) {
#include "Elasticity.test.inc"

  // set up matrix structures
  double l_starE[3][24];
  double l_jacInv[3][3] = { {0.4, 0.3, 0.5}, {0.5, 0.1, 0.3}, {0.2, 0.4, 0.2} };

  // compute star matrices
  edge::seismic::setups::Elasticity::star( 2600,
                                           20.8E9,
                                           9.9E9,
                                           l_jacInv,
                                           l_starE );

  // check results
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    unsigned short l_nz = 0;
    for( unsigned short l_m0 = 0; l_m0 < 9; l_m0++ ) {
      for( unsigned short l_m1 = 0; l_m1 < 9; l_m1++ ) {
        if( l_starE3d[l_di][l_m0][l_m1] != 0 ) {
          REQUIRE( l_starE[l_di][l_nz] == Approx(l_starE3d[l_di][l_m0][l_m1]) );
          l_nz++;
        }
      }
    }
    REQUIRE( l_nz == 24 );
  }
}
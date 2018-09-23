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
 * Unit tests for operations on series.
 **/

#include <catch.hpp>
#include "Series.hpp"

TEST_CASE( "Series: Integrate", "[ts][integrate]" ) {
  /*
   *
   *   vals:    0.5       1.2       -0.3      0.0     -0.1
   *   diff:        0.7      -1.5        0.3     -0.1
   *   diff/2:      0.35     -0.75       0.15    -0.05
   *
   *   sum first four vals:  1.4
   *   sum diff:            -0.6
   *   sum diff/2:          -0.3
   */

  // input values
  double l_vals[5][8] = {
    {  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 },
    {  1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2 },
    { -0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3 },
    {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1 }
  };

  // result
  double l_re[8];

  // all points, dt=1.0
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.0,
                                        4.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(1.1) );
  }

  // all points, dt=0.1
  edge::linalg::Series< 8 >::integrate( 0.1,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.0,
                                        4.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx( (1.1)*0.1) );
  }

  // all points, except first
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        1.0,
                                        4.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.25) );
  }

  // all points, except last
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.0,
                                        3.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
   REQUIRE( l_re[l_se] == Approx(1.15) );
  }

  // only middle points
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        1.0,
                                        3.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.3) );
  }

  // only first interval
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.0,
                                        1.0,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.85) );
  }

  // subinterval in 2nd interval
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.2,
                                        1.3,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(1.0285) );
  }

  // subinterval in 2nd interval, scaled
  edge::linalg::Series< 8 >::integrate( 0.1,
                                        0.0,
                                        5,
                                        l_vals,
                                        0.02,
                                        0.13,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.10285) );
  }

  // start before series
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        -0.5,
                                        1.3,
                                        l_re );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(1.1425) );
  }

  // start and end before series
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        -7.0,
                                        -3.0,
                                        l_re,
                                        0.5 );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(2.0) );
  }

  // start before series
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        1.0,
                                        5,
                                        l_vals,
                                        3.9,
                                        8.0,
                                        l_re,
                                        0.1 );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.2485) );
  }

  // start and end after series
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        6.0,
                                        9.0,
                                        l_re,
                                        0.1 );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(0.3) );
  }

  // start before end end after series
  edge::linalg::Series< 8 >::integrate( 1.0,
                                        0.0,
                                        5,
                                        l_vals,
                                        -5.0,
                                        9.0,
                                        l_re,
                                        0.3 );
  for( unsigned short l_se = 0; l_se < 8; l_se++ ) {
    REQUIRE( l_re[l_se] == Approx(4.1) );
  }
}

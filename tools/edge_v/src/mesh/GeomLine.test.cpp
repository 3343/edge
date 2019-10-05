/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Tests the geometry computations for line elements.
 **/
#include <catch.hpp>
#define private public
#include "GeomLine.h"
#undef private

TEST_CASE( "Tests the volume computation for lines.", "[volume][line]" ) {
  double l_ves[2][3] = { {1, 0, 2}, {2, 1, 3} };

  double l_vol = edge_v::geom::Line::length( l_ves );

  REQUIRE( l_vol == Approx( std::sqrt(3.0) ) );
}

TEST_CASE( "Tests the normal computation for lines.", "[normal][line]" ) {
  double l_ves0[2][3] = { {0.0, 0.0, 0.0},
                          {7.0, 0.0, 0.0} };
  double l_np0[3] = {0.0, 10.0, 0.0};

  double l_normal[3] = {0};

  edge_v::geom::Line::normal( l_ves0, l_np0, l_normal );

  REQUIRE( l_normal[0] == Approx( 0.0) );
  REQUIRE( l_normal[1] == Approx(-1.0) );
  REQUIRE( l_normal[2] == Approx( 0.0) );

  double l_np1[3] = {0.0, -10.0, 0.0};

  edge_v::geom::Line::normal( l_ves0, l_np1, l_normal );

  REQUIRE( l_normal[0] == Approx( 0.0) );
  REQUIRE( l_normal[1] == Approx( 1.0) );
  REQUIRE( l_normal[2] == Approx( 0.0) );
}

TEST_CASE( "Tests the tangent computation for lines.", "[tangent][line]" ) {
  double l_ves0[2][3] = { {0.0, 0.0, 0.0},
                          {7.0, 0.0, 0.0} };
  double l_np0[3] = {0.0, 10.0, 0.0};

  double l_tangent[3] = {0};

  edge_v::geom::Line::tangent( l_ves0, l_np0, l_tangent );

  REQUIRE( l_tangent[0] == Approx( 1.0) );
  REQUIRE( l_tangent[1] == Approx( 0.0) );
  REQUIRE( l_tangent[2] == Approx( 0.0) );

  double l_np1[3] = {0.0, -10.0, 0.0};
  edge_v::geom::Line::tangent( l_ves0, l_np1, l_tangent );

  REQUIRE( l_tangent[0] == Approx(-1.0) );
  REQUIRE( l_tangent[1] == Approx( 0.0) );
  REQUIRE( l_tangent[2] == Approx( 0.0) );
}
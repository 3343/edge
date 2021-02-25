/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Tests the geometry computations for quad4r elements.
 **/
#include <catch.hpp>
#define private public
#include "Quad4r.h"
#undef private

TEST_CASE( "Tests the area computation for quad4r elements.", "[quad4r][area]" ) {
  double l_ves[4][3] = { {1, 0, 0}, {5, 0, 0}, {5, 10, 0}, {1, 10, 0} };

  double l_area = edge_v::geom::Quad4r::area( l_ves );

  REQUIRE( l_area == Approx( 40 ) );
}

TEST_CASE( "Tests the incircle diameter computation for quad4r elements.", "[quad4r][inDiameter]" ) {
  double l_ves[4][3] = { {1, 0, 0}, {5, 0, 0}, {5, 10, 0}, {1, 10, 0} };

  double l_dia = edge_v::geom::Quad4r::inDiameter( l_ves );

  REQUIRE( l_dia == Approx( 4 ) );
}

TEST_CASE( "Tests the normalization of the elements' vertex-orders.", "[quad4r][normVe]" ) {
  /*
   *        0
   *  |0---------2|
   *  |           |
   * 1|           |2
   *  |           |
   *  |3_________1|
   *        3
   */

  double l_ves[4][3] = { {1, 10, 0}, {5, 0, 0}, {5, 10, 0}, {1, 0, 0} };
  edge_v::t_idx l_elVe[4]   = { 0,  1,  2,  3 };
  edge_v::t_idx l_elFa[4]   = { 4,  5,  6,  7 };
  edge_v::t_idx l_elFaEl[4] = { 8,  9, 10, 11 };

  edge_v::geom::Quad4r::normVesFas( l_ves,
                                    l_elVe,
                                    l_elFa,
                                    l_elFaEl );

  REQUIRE( l_elVe[0]   ==  3 );
  REQUIRE( l_elVe[1]   ==  1 );
  REQUIRE( l_elVe[2]   ==  2 );
  REQUIRE( l_elVe[3]   ==  0 );

  REQUIRE( l_elFa[0]   ==  7 );
  REQUIRE( l_elFa[1]   ==  6 );
  REQUIRE( l_elFa[2]   ==  4 );
  REQUIRE( l_elFa[3]   ==  5 );

  REQUIRE( l_elFaEl[0] == 11 );
  REQUIRE( l_elFaEl[1] == 10 );
  REQUIRE( l_elFaEl[2] ==  8 );
  REQUIRE( l_elFaEl[3] ==  9 );
}
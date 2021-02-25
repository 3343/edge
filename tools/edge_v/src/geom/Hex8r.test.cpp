/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Tests the geometry computations for hex8r elements.
 **/
#include <catch.hpp>
#define private public
#include "Hex8r.h"
#undef private

TEST_CASE( "Tests the volume computation for hex8r elements.", "[hex8r][volume]" ) {
  double l_ves[8][3] = { {1, 0, 0}, {5, 0, 0},
                         {5, 10, 0}, {1, 10, 0},

                         {1, 0, 3}, {5, 0, 3},
                         {5, 10, 3}, {1, 10, 3} };

  double l_volume = edge_v::geom::Hex8r::volume( l_ves );

  REQUIRE( l_volume == Approx( 120 ) );
}

TEST_CASE( "Tests the insphere diameter computation for hex8r elements.", "[hex8r][inDiameter]" ) {
  double l_ves[8][3] = { {1, 0, 0}, {5, 0, 0},
                         {5, 10, 0}, {1, 10, 0},

                         {1, 0, 3}, {5, 0, 3},
                         {5, 10, 3}, {1, 10, 3} };

  double l_dia = edge_v::geom::Hex8r::inDiameter( l_ves );

  REQUIRE( l_dia == Approx( 3 ) );
}

TEST_CASE( "Tests the normalization of the hex8r element's vertex-orders.", "[hex8r][normVesFas]" ) {
  double l_ves[8][3] = { {1, 0, 3}, {5, 0, 3},
                         {1, 10, 3}, {5, 10, 3},

                         {1, 0, 0}, {5, 0, 0},
                         {5, 10, 0}, {1, 10, 0} };

  edge_v::t_idx l_elVe[8]   = { 2, 3, 4, 5,
                                6, 7, 8, 9 };
  //                              0   1   2   3   4   5
  edge_v::t_idx l_elFa[6]   = {   8,  9, 10, 11, 12, 13 };
  edge_v::t_idx l_elFaEl[6] = {  14, 15, 16, 17, 18, 19 };

  edge_v::geom::Hex8r::normVesFas( l_ves,
                                   l_elVe,
                                   l_elFa,
                                   l_elFaEl );

  // bottom
  REQUIRE( l_elVe[0] == 6 );
  REQUIRE( l_elVe[1] == 7 );
  REQUIRE( l_elVe[2] == 8 );
  REQUIRE( l_elVe[3] == 9 );

  // top
  REQUIRE( l_elVe[4] == 2 );
  REQUIRE( l_elVe[5] == 3 );
  REQUIRE( l_elVe[6] == 5 );
  REQUIRE( l_elVe[7] == 4 );

  /*
   *           ref-ids   in-ids    lsort     id 
   *   face 0: 0-3-2-1 - 6-7-8-9 - 6-7-8-9 - 5
   *   face 1: 0-1-5-4 - 6-7-3-2 - 2-3-6-7 - 1
   *   face 2: 1-2-6-5 - 7-8-5-3 - 3-5-7-8 - 3
   *   face 3: 3-7-6-2 - 9-4-5-8 - 4-5-8-9 - 4
   *   face 4: 0-4-7-3 - 6-2-4-9 - 2-4-6-9 - 2
   *   face 5: 4-5-6-7 - 2-3-4-5 - 2-3-4-5 - 0
   */
  REQUIRE( l_elFa[0] == 13 );
  REQUIRE( l_elFaEl[0] == 19 );

  REQUIRE( l_elFa[1] == 9 );
  REQUIRE( l_elFaEl[1] == 15 );

  REQUIRE( l_elFa[2] == 11 );
  REQUIRE( l_elFaEl[2] == 17 );

  REQUIRE( l_elFa[3] == 12 );
  REQUIRE( l_elFaEl[3] == 18 );

  REQUIRE( l_elFa[4] == 10 );
  REQUIRE( l_elFaEl[4] == 16 );

  REQUIRE( l_elFa[5] == 8 );
  REQUIRE( l_elFaEl[5] == 14 );
}
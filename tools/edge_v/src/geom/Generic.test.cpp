/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Tests the generic geometry functions.
 **/
#include <catch.hpp>
#define private public
#include "Generic.h"
#undef private

TEST_CASE( "Tests the centroid computation", "[generic][centroid]" ) {
  double l_veCrds[4][3] = { { 0,  1,  2 },
                            { 3,  4,  5 },
                            { 6,  7,  8 },
                            { 9, 10, 11 } };

  double l_cen[3] = {0};

  edge_v::geom::Generic::centroid( edge_v::t_entityType::TET4,
                                   l_veCrds,
                                   l_cen );

  REQUIRE( l_cen[0] == Approx( 4.5 ) );
  REQUIRE( l_cen[1] == Approx( 5.5 ) );
  REQUIRE( l_cen[2] == Approx( 6.5 ) );
}

TEST_CASE( "Tests the lexicographic sorting.", "[generic][sortLex]" ) {
  edge_v::t_idx l_data[7*3] = { 2, 3, 1,   // 3
                                1, 4, 5,   // 1
                                3, 5, 1,   // 5
                                0, 3, 4,   // 0
                                2, 4, 0,   // 4
                                1, 4, 6,   // 2
                                6, 7, 9 }; // 6

  edge_v::t_idx l_ref[7*3] = { 0, 3, 4,   // 0
                               1, 4, 5,   // 1
                               1, 4, 6,   // 2
                               2, 3, 1,   // 3
                               2, 4, 0,   // 4
                               3, 5, 1,   // 5
                               6, 7, 9 }; // 6

  edge_v::geom::Generic::sortLex( 7,
                                  3,
                                  l_data,
                                 nullptr  );

  for( unsigned short l_en = 0; l_en < 7*3; l_en++ ) {
    REQUIRE( l_data[l_en] == l_ref[l_en] );
  }
}
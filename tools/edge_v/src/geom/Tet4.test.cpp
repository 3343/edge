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
 * Tests the geometry computations on tet4-elements.
 **/
#include <catch.hpp>
#define private public
#include "Tet4.h"
#undef private

TEST_CASE( "Tests the volume computation of 4-node tetrahedrons.", "[volume][tet4]" ) {
    // construct tet
  double l_ves[4][3] = { {0.3, 0.1, 0.2},
                         {1.1, 0.4, 0.3},
                         {0.5, 1.1, 0.4},
                         {0.5, 0.5, 1.2} };

  // compute diameter
  double l_vol = edge_v::geom::Tet4::volume( l_ves );

  // check results
  REQUIRE( l_vol == Approx(0.112667) );
}

TEST_CASE( "Tests the insphere diameter computation for 4-node tetrahedrons.", "[diameter][tet4]" ) { 
  // construct tet
  double l_ves[4][3] = { {0.0, 0.0, 0.0},
                         {1.0, 0.0, 0.0},
                         {0.5, 1.0, 0.0},
                         {0.5, 0.5, 1.0} };

  // compute diameter
  double l_dia = edge_v::geom::Tet4::inDiameter( l_ves );

  // check result
  REQUIRE( l_dia == Approx(0.45358449083204) );
}

TEST_CASE( "Tests the derivation of vertex-ids for 4-node tetrahedrons.", "[getVeIdsAd][tet4]" ) {
  // test derivation through ids
  edge_v::t_idx l_el0[2] = {0, 1};
  unsigned short l_fa0[2] = {3, 1};
  edge_v::t_idx l_elVe0[8] = { 0, 1, 2, 3,
                               3, 1, 5, 2 };
  edge_v::t_idx l_elFaEl0[8] = { 3343, 3343, 3343, 1,
                                 3343, 0,    3343, 3343 };
  unsigned short l_veIdsAd0[2] = { 3343, 3343 };

  edge_v::geom::Tet4::getVeIdsAd( 2,
                                  0,
                                  l_el0,
                                  l_fa0,
                                  l_elVe0,
                                  l_elFaEl0,
                                  nullptr,
                                  l_veIdsAd0 );

  REQUIRE( l_veIdsAd0[0] == 1 );
  REQUIRE( l_veIdsAd0[1] == 2 );

  // test derivation through coordinates
  edge_v::t_idx l_el1[2] = {0, 1};
  unsigned short l_fa1[2] = {3, 1};
  edge_v::t_idx l_elVe1[8] = { 0, 1, 2, 3,
                               3, 6, 5, 2 };
  edge_v::t_idx l_elFaEl1[8] = { 3343, 3343, 3343, 1,
                                 3343, 0,    3343, 3343 };
  unsigned short l_veIdsAd1[2] = { 3343, 3343 };

  double l_veCrds1[7][3] = { { 0,  1, 3343 },   // 0
                             { 7, -5, 3343 },   // 1
                             { 2,  3, 3343 },   // 2
                             { 4,  5, 3343 },   // 3
                             { 6,  7, 6686 },   // 4
                             { 8,  9, 6686 },   // 5
                             { 7, -5, 6686 } }; // 6

  edge_v::geom::Tet4::getVeIdsAd( 2,
                                  0,
                                  l_el1,
                                  l_fa1,
                                  l_elVe1,
                                  l_elFaEl1,
                                  l_veCrds1,
                                  l_veIdsAd1 );

  REQUIRE( l_veIdsAd1[0] == 1 );
  REQUIRE( l_veIdsAd1[1] == 2 );
}
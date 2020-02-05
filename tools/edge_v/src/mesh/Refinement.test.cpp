/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Tests the velocity-based mesh refinement.
 **/
#include <catch.hpp>
#include "models/Constant.h"
#define private public
#include "Refinement.h"
#undef private

TEST_CASE( "Tests the mesh refinement.", "[meshRefinement][init]" ) {
  double l_crds[4][3] = { { 1, 2,   3},
                          {-2, 3,   4},
                          {-4, -2,  1},
                          { 5, 10, 20} };

  std::size_t l_elVe[4][4] = { {0, 0, 0, 0},
                               {1, 1, 1, 1},
                               {2, 2, 2, 2},
                               {3, 3, 3, 3} };

  std::string l_exp0 = "edges_per_wave_length := 4.3;\
                        frequency := 2.7;";

  edge_v::models::Constant l_mod(3);

  edge_v::mesh::Refinement l_ref;
  l_ref.init( 4,
              4,
              4,
              l_elVe[0],
              l_crds,
              l_exp0,
              l_mod );

  REQUIRE( l_ref.m_refVe[0] == Approx( 3.0 / (2.7*4.3) ) );
  REQUIRE( l_ref.m_refVe[1] == Approx( 3.0 / (2.7*4.3) ) );
  REQUIRE( l_ref.m_refVe[2] == Approx( 3.0 / (2.7*4.3) ) );
  REQUIRE( l_ref.m_refVe[3] == Approx( 3.0 / (2.7*4.3) ) );

  std::string l_exp1 = "if( x > 0 ) {\
                          edges_per_wave_length := 4.3;\
                        }\
                        else {\
                          edges_per_wave_length := 1.3;\
                        }\
                        if( y < 3 ) {\
                          frequency := 2.7;\
                        }\
                        else {\
                          frequency := 9.3;\
                        }";

  l_ref.init( 4,
              4,
              4,
              l_elVe[0],
              l_crds,
              l_exp1,
              l_mod );

  REQUIRE( l_ref.m_refVe[0] == Approx( 3.0 / (2.7*4.3) ) );
  REQUIRE( l_ref.m_refVe[1] == Approx( 3.0 / (9.3*1.3) ) );
  REQUIRE( l_ref.m_refVe[2] == Approx( 3.0 / (2.7*1.3) ) );
  REQUIRE( l_ref.m_refVe[3] == Approx( 3.0 / (9.3*4.3) ) );
}
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
 * Tests the expression-based seismic velocity model.
 **/
#include <catch.hpp>
#define private public
#include "Expression.h"
#undef private

#include "io/logging.h"

TEST_CASE( "Tests the initialization.", "[expression][init]" ) {
  std::string l_exp0 = "vp := x;\
                        vs := y;\
                        qp := z + vs;\
                        qs := vp - vs;";
  edge_v::models::seismic::Expression l_veMod0( l_exp0 );

  double l_pts0[2][3] = { {1.0, 2.0, 3.0},
                          {4.0, 5.0, 6.0} };

  l_veMod0.init( 2,
                 l_pts0 );

  // check the results
  REQUIRE( l_veMod0.m_vp[0] == Approx(1.0) );
  REQUIRE( l_veMod0.m_vs[0] == Approx(2.0) );
  REQUIRE( l_veMod0.m_qp[0] == Approx(3.0+2.0) );
  REQUIRE( l_veMod0.m_qs[0] == Approx(1.0-2.0) );

  REQUIRE( l_veMod0.m_vp[1] == Approx(4.0) );
  REQUIRE( l_veMod0.m_vs[1] == Approx(5.0) );
  REQUIRE( l_veMod0.m_qp[1] == Approx(6.0+5.0) );
  REQUIRE( l_veMod0.m_qs[1] == Approx(4.0-5.0) );

  std::string l_exp1 = "if( z > -3.5 ) {\
                          vp := 0.2;\
                          vs := 0.3;\
                          qp := 0.4;\
                          qs := 0.5;\
                        }\
                        else {\
                          vp := 0.6;\
                          vs := 0.7;\
                          qp := 0.8;\
                          qs := 0.9;\
                        }";

  double l_pts1[5][3] = { {0.0, 1.0, -3.0},
                          {2.0, 4.0, -4.0},
                          {9.0, 2.0, -5.0},
                          {1.0, 3.0, -2.0},
                          {3.0, 4.0, -4.0} };

  edge_v::models::seismic::Expression l_veMod1( l_exp1 );
  l_veMod1.init( 5,
                 l_pts1 );  

  // check the results
  REQUIRE( l_veMod1.m_vp[0] == Approx(0.2) );
  REQUIRE( l_veMod1.m_vs[0] == Approx(0.3) );
  REQUIRE( l_veMod1.m_qp[0] == Approx(0.4) );
  REQUIRE( l_veMod1.m_qs[0] == Approx(0.5) );

  REQUIRE( l_veMod1.m_vp[1] == Approx(0.6) );
  REQUIRE( l_veMod1.m_vs[1] == Approx(0.7) );
  REQUIRE( l_veMod1.m_qp[1] == Approx(0.8) );
  REQUIRE( l_veMod1.m_qs[1] == Approx(0.9) );

  REQUIRE( l_veMod1.m_vp[2] == Approx(0.6) );
  REQUIRE( l_veMod1.m_vs[2] == Approx(0.7) );
  REQUIRE( l_veMod1.m_qp[2] == Approx(0.8) );
  REQUIRE( l_veMod1.m_qs[2] == Approx(0.9) );

  REQUIRE( l_veMod1.m_vp[3] == Approx(0.2) );
  REQUIRE( l_veMod1.m_vs[3] == Approx(0.3) );
  REQUIRE( l_veMod1.m_qp[3] == Approx(0.4) );
  REQUIRE( l_veMod1.m_qs[3] == Approx(0.5) );

  REQUIRE( l_veMod1.m_vp[4] == Approx(0.6) );
  REQUIRE( l_veMod1.m_vs[4] == Approx(0.7) );
  REQUIRE( l_veMod1.m_qp[4] == Approx(0.8) );
  REQUIRE( l_veMod1.m_qs[4] == Approx(0.9) );
}
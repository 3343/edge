/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
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
 * Tests the UCVM-derived seismic velocity model.
 **/
#include <catch.hpp>
#define private public
#include "Ucvm.h"
#undef private

TEST_CASE( "Tests the UCVM model abstraction.", "[ucvm][model]" ) {
  std::string l_cfg = PP_UCVM_CONF;
  std::string l_models = "cvms";
  std::string l_mode = "UCVM_COORD_GEO_DEPTH";

  edge_v::io::Ucvm l_ucvmReader( l_cfg,
                                 l_models,
                                 l_mode );

  double l_trafoSrc[3][3] = { {1, 0, 0},
                              {0, 1, 0},
                              {0, 0, 1} };
  std::string l_projSrc = "+proj=tmerc +units=m +axis=eun +no_defs +datum=WGS84 +k=0.9996 +lon_0=-118 +lat_0=34";
  std::string l_projDes = "+proj=longlat +datum=WGS84";
  std::string l_ucvmType = "crust";

  edge_v::models::seismic::Ucvm l_ucvm( l_ucvmReader,
                                        l_trafoSrc,
                                        l_projSrc,
                                        l_projDes,
                                        l_ucvmType );

  double l_pts[5][3] = { {0, 0,     0},
                         {0, 0,    50},
                         {0, 0,   100},
                         {0, 0,   500},
                         {0, 0,  1000} };

  l_ucvm.init( 5,
               l_pts );

  // check results (analogue to io unit tests)
  REQUIRE( l_ucvm.m_velP[0] == Approx(  696.491 ) );
  REQUIRE( l_ucvm.m_velP[1] == Approx( 1669.540 ) );
  REQUIRE( l_ucvm.m_velP[2] == Approx( 1683.174 ) );
  REQUIRE( l_ucvm.m_velP[3] == Approx( 2701.216 ) );
  REQUIRE( l_ucvm.m_velP[4] == Approx( 3330.909 ) );

  REQUIRE( l_ucvm.m_velS[0] == Approx(  213.000 ) );
  REQUIRE( l_ucvm.m_velS[1] == Approx(  548.000 ) );
  REQUIRE( l_ucvm.m_velS[2] == Approx(  603.470 ) );
  REQUIRE( l_ucvm.m_velS[3] == Approx( 1475.608 ) );
  REQUIRE( l_ucvm.m_velS[4] == Approx( 1945.594 ) );

  REQUIRE( l_ucvm.m_rho[0] == Approx( 1974.976 ) );
  REQUIRE( l_ucvm.m_rho[1] == Approx( 2128.620 ) );
  REQUIRE( l_ucvm.m_rho[2] == Approx( 2130.773 ) );
  REQUIRE( l_ucvm.m_rho[3] == Approx( 2354.105 ) );
  REQUIRE( l_ucvm.m_rho[4] == Approx( 2443.042 ) );

  REQUIRE( l_ucvm.getMaxSpeed(0) == Approx(  696.491 ) );
  REQUIRE( l_ucvm.getMaxSpeed(1) == Approx( 1669.540 ) );
  REQUIRE( l_ucvm.getMaxSpeed(2) == Approx( 1683.174 ) );
  REQUIRE( l_ucvm.getMaxSpeed(3) == Approx( 2701.216 ) );
  REQUIRE( l_ucvm.getMaxSpeed(4) == Approx( 3330.909 ) );
}
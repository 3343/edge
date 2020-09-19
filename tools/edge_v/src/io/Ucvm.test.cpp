/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Tests the UCVM interface.
 **/
#include <catch.hpp>
#define private public
#include "Ucvm.h"
#undef private

TEST_CASE( "Tests CVM-S4 queries via UCVM.", "[ucvm][cvms]" ) {
  // derivation of reference data
#if 0
[alex@jenalex ucvm-19.4.0]$ echo -e "-118 34 0\n-118 34 50\n-118 34 100\n-118 34 500\n-118 34 1000" |  ./bin/ucvm_query -f conf/ucvm.conf -m cvmsi
Using Geo Depth coordinates as default mode.
 -118.0000    34.0000      0.000    281.668    468.400      cvmsi    696.491    213.000   1974.976       none      0.000      0.000      0.000      crust    696.491    213.000   1974.976
 -118.0000    34.0000     50.000    281.668    468.400      cvmsi   1669.540    548.000   2128.620       none      0.000      0.000      0.000      crust   1669.540    548.000   2128.620
 -118.0000    34.0000    100.000    281.668    468.400      cvmsi   1683.174    603.470   2130.773       none      0.000      0.000      0.000      crust   1683.174    603.470   2130.773
 -118.0000    34.0000    500.000    281.668    468.400      cvmsi   2701.216   1475.608   2354.105       none      0.000      0.000      0.000      crust   2701.216   1475.608   2354.105
 -118.0000    34.0000   1000.000    281.668    468.400      cvmsi   3330.909   1945.594   2443.042       none      0.000      0.000      0.000      crust   3330.909   1945.594   2443.042
#endif

  std::string l_cfg = PP_UCVM_CONF;
  std::string l_models = "cvmsi";
  std::string l_mode = "UCVM_COORD_GEO_DEPTH";

  double l_trafoSrc[3][3] = { {1, 0, 0},
                              {0, 1, 0},
                              {0, 0, 1} };

  std::string l_projSrc = "+proj=tmerc +units=m +axis=eun +no_defs +datum=WGS84 +k=0.9996 +lon_0=-118 +lat_0=34";
  std::string l_projDes = "+proj=longlat +datum=WGS84";
  std::string l_ucvmType = "crust";

  double l_pts[5][3] = { {0, 0,     0},
                         {0, 0,    50},
                         {0, 0,   100},
                         {0, 0,   500},
                         {0, 0,  1000} };
  float l_vps[5];
  float l_vss[5];
  float l_rhos[5];

  edge_v::io::Ucvm l_ucvm( l_cfg,
                           l_models,
                           l_mode );

  l_ucvm.getVels( 5,
                  l_trafoSrc,
                  l_projSrc,
                  l_projDes,
                  l_ucvmType,
                  l_pts,
                  l_vps,
                  l_vss,
                  l_rhos );

  REQUIRE( l_vps[0] == Approx(  696.491 ) );
  REQUIRE( l_vps[1] == Approx( 1669.540 ) );
  REQUIRE( l_vps[2] == Approx( 1683.174 ) );
  REQUIRE( l_vps[3] == Approx( 2701.216 ) );
  REQUIRE( l_vps[4] == Approx( 3330.909 ) );

  REQUIRE( l_vss[0] == Approx(  213.000 ) );
  REQUIRE( l_vss[1] == Approx(  548.000 ) );
  REQUIRE( l_vss[2] == Approx(  603.470 ) );
  REQUIRE( l_vss[3] == Approx( 1475.608 ) );
  REQUIRE( l_vss[4] == Approx( 1945.594 ) );

  REQUIRE( l_rhos[0] == Approx( 1974.976 ) );
  REQUIRE( l_rhos[1] == Approx( 2128.620 ) );
  REQUIRE( l_rhos[2] == Approx( 2130.773 ) );
  REQUIRE( l_rhos[3] == Approx( 2354.105 ) );
  REQUIRE( l_rhos[4] == Approx( 2443.042 ) );
}
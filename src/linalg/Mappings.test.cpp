/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Unit tests for the mappings.
 **/
#include <catch.hpp>

#include "constants.hpp"
#include "Mappings.hpp"

TEST_CASE( "Mappings: Physical to reference coords, tet4", "[mappings][phyToRefTet4]" ) {
  // physical tet
  double l_ves[3][4] = { {-3,-1,-3,-3},
                         {-3,-3,-1,-3},
                         {-3,-3,-3,-1} };
  double l_pt[3];
  double l_res[3];

  l_pt[0] = -2.5;
  l_pt[1] = -2.5;
  l_pt[2] = -2.5;
  edge::linalg::Mappings::phyToRefTet4( l_ves, l_pt, l_res );
  REQUIRE( l_res[0] == Approx(0.25) );
  REQUIRE( l_res[1] == Approx(0.25) );
  REQUIRE( l_res[2] == Approx(0.25) );

  l_pt[0] = -3;
  l_pt[1] = -1;
  l_pt[2] = -6;
  edge::linalg::Mappings::phyToRefTet4( l_ves, l_pt, l_res );
  REQUIRE( l_res[0] == Approx( 0.0) );
  REQUIRE( l_res[1] == Approx( 1.0) );
  REQUIRE( l_res[2] == Approx(-1.5) );

  l_pt[0] = -3;
  l_pt[1] = -2.5;
  l_pt[2] = -6;
  edge::linalg::Mappings::phyToRef( TET4, l_ves[0], l_pt, l_res );
  REQUIRE( l_res[0] == Approx( 0.0) );
  REQUIRE( l_res[1] == Approx( 0.25) );
  REQUIRE( l_res[2] == Approx(-1.5) );
}

TEST_CASE( "Mappings: Face to volume coordinates (reference element)", "[mappings][faToVolRef]" ) {
  // init vars
  unsigned short l_fa;
  t_entityType l_elType;
  double l_faPt[2];
  double l_volPt[3];

  // line elements
  l_fa = 0;
  l_elType = LINE;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.0 ) );

  l_fa = 1;
  l_elType = LINE;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 1.0 ) );

  // quad4r elements
  l_fa = 0;
  l_faPt[0] = 0.45;
  l_elType = QUAD4R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.45 ) );
  REQUIRE( l_volPt[1] == Approx( 0.0  ) );

  l_fa = 1;
  l_faPt[0] = 0.25;
  l_elType = QUAD4R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 1.0  ) );
  REQUIRE( l_volPt[1] == Approx( 0.25 ) );

  l_fa = 2;
  l_faPt[0] = 0.3;
  l_elType = QUAD4R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.7 ) );
  REQUIRE( l_volPt[1] == Approx( 1.0 ) );

  l_fa = 3;
  l_faPt[0] = 0.7;
  l_elType = QUAD4R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.0 ) );
  REQUIRE( l_volPt[1] == Approx( 0.3 ) );

  // tria3 elements
  l_fa = 0;
  l_faPt[0] = 0.2;
  l_elType = TRIA3;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.2 ) );
  REQUIRE( l_volPt[1] == Approx( 0.0 ) );

  l_fa = 1;
  l_faPt[0] = 0.9;
  l_elType = TRIA3;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.1 ) );
  REQUIRE( l_volPt[1] == Approx( 0.9 ) );

  l_fa = 2;
  l_faPt[0] = 0.4;
  l_elType = TRIA3;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.0 ) );
  REQUIRE( l_volPt[1] == Approx( 0.6 ) );

  // hex8r elements
  l_fa = 0;
  l_faPt[0] = 0.2;
  l_faPt[1] = 0.7;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.7 ) );
  REQUIRE( l_volPt[1] == Approx( 0.2 ) );
  REQUIRE( l_volPt[2] == Approx( 0.0 ) );

  l_fa = 1;
  l_faPt[0] = 0.3;
  l_faPt[1] = 0.4;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.3 ) );
  REQUIRE( l_volPt[1] == Approx( 0.0 ) );
  REQUIRE( l_volPt[2] == Approx( 0.4 ) );

  l_fa = 2;
  l_faPt[0] = 0.9;
  l_faPt[1] = 0.8;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 1.0 ) );
  REQUIRE( l_volPt[1] == Approx( 0.9 ) );
  REQUIRE( l_volPt[2] == Approx( 0.8 ) );

  l_fa = 3;
  l_faPt[0] = 0.4;
  l_faPt[1] = 0.7;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.7 ) );
  REQUIRE( l_volPt[1] == Approx( 1.0 ) );
  REQUIRE( l_volPt[2] == Approx( 0.4 ) );

  l_fa = 4;
  l_faPt[0] = 0.5;
  l_faPt[1] = 0.6;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.0 ) );
  REQUIRE( l_volPt[1] == Approx( 0.6 ) );
  REQUIRE( l_volPt[2] == Approx( 0.5 ) );

  l_fa = 5;
  l_faPt[0] = 0.1;
  l_faPt[1] = 0.6;
  l_elType = HEX8R;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.1 ) );
  REQUIRE( l_volPt[1] == Approx( 0.6 ) );
  REQUIRE( l_volPt[2] == Approx( 1.0 ) );

  // tet4 elements
  l_fa = 0;
  l_faPt[0] = 0.2;
  l_faPt[1] = 0.3;
  l_elType = TET4;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.3 ) );
  REQUIRE( l_volPt[1] == Approx( 0.2 ) );
  REQUIRE( l_volPt[2] == Approx( 0.0 ) );

  l_fa = 1;
  l_faPt[0] = 0.7;
  l_faPt[1] = 0.2;
  l_elType = TET4;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.7 ) );
  REQUIRE( l_volPt[1] == Approx( 0.0 ) );
  REQUIRE( l_volPt[2] == Approx( 0.2 ) );

  l_fa = 2;
  l_faPt[0] = 0.4;
  l_faPt[1] = 0.9;
  l_elType = TET4;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.0 ) );
  REQUIRE( l_volPt[1] == Approx( 0.9 ) );
  REQUIRE( l_volPt[2] == Approx( 0.4 ) );

  l_fa = 3;
  l_faPt[0] = 0.3;
  l_faPt[1] = 0.2;
  l_elType = TET4;
  edge::linalg::Mappings::faToVolRef( l_fa, l_elType, l_faPt, l_volPt );
  REQUIRE( l_volPt[0] == Approx( 0.5 ) );
  REQUIRE( l_volPt[1] == Approx( 0.3 ) );
  REQUIRE( l_volPt[2] == Approx( 0.2 ) );
}

TEST_CASE( "Mappings: Local face coordinates to neighboring face coordinates", "[mappings][faLocToFaNei]" ) {
  unsigned short l_ve;
  real_mesh l_faPtL[2];
  real_mesh l_faPtN[2];

  /*
   * Tria3
   */
  l_faPtL[0] = 0.2;
  edge::linalg::Mappings::faLocToFaNei( TRIA3, l_faPtL, l_faPtN );
  REQUIRE( l_faPtN[0] == Approx(0.8) );

  l_faPtL[0] = 0.7;
  edge::linalg::Mappings::faLocToFaNei( TRIA3, l_faPtL, l_faPtN );
  REQUIRE( l_faPtN[0] == Approx(0.3) );

  /*
   * Hex8r
   */
  l_faPtL[0] = 0.4;
  l_faPtL[1] = 0.3;
  edge::linalg::Mappings::faLocToFaNei( HEX8R, l_faPtL, l_faPtN );
  REQUIRE( l_faPtN[0] == Approx(0.3) );
  REQUIRE( l_faPtN[1] == Approx(0.4) );

  /*
   * Tet4
   */
  l_ve = 0;
  l_faPtL[0] = 0.2;
  l_faPtL[1] = 0.7;
  edge::linalg::Mappings::faLocToFaNei( TET4, l_faPtL, l_faPtN, l_ve );
  REQUIRE( l_faPtN[0] == Approx(0.7) );
  REQUIRE( l_faPtN[1] == Approx(0.2) );

  l_ve = 1;
  l_faPtL[0] = 0.3;
  l_faPtL[1] = 0.6;
  edge::linalg::Mappings::faLocToFaNei( TET4, l_faPtL, l_faPtN, l_ve );
  REQUIRE( l_faPtN[0] == Approx(0.1) );
  REQUIRE( l_faPtN[1] == Approx(0.6) );

  l_ve = 2;
  l_faPtL[0] = 0.2;
  l_faPtL[1] = 0.1;
  edge::linalg::Mappings::faLocToFaNei( TET4, l_faPtL, l_faPtN, l_ve );
  REQUIRE( l_faPtN[0] == Approx(0.2) );
  REQUIRE( l_faPtN[1] == Approx(0.7) );
}

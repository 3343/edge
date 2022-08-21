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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONsTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Unit tests of geometry module.
 **/
#include <catch.hpp>
#include "constants.hpp"

#include "Geom.hpp"
#include "Matrix.h"

TEST_CASE( "Inside/outside derivation: Points with respect to lines", "[geom][insideLine]" ) {
  real_mesh l_veCrds[2] = {1,2}, l_pt;

  l_pt = 0.5;
  REQUIRE( edge::linalg::Geom::inside( LINE, l_veCrds, &l_pt ) == 0 );

  l_pt = 1.5;
  REQUIRE( edge::linalg::Geom::inside( LINE, l_veCrds, &l_pt ) == 1 );

  l_pt = 2.5;
  REQUIRE( edge::linalg::Geom::inside( LINE, l_veCrds, &l_pt ) == 0 );

  l_pt = 2;
  REQUIRE( edge::linalg::Geom::inside( LINE, l_veCrds, &l_pt ) == 2 );
}

TEST_CASE( "Projection of a point outside of an entity to the surface of the entity.", "[geom][outProjSurf]" ) {
  double l_pt[3] = { 99, 99, 99 };

  double l_line[2] = {0, 1};

  l_pt[0] = -2;
  edge::linalg::Geom::closestPoint( LINE, l_line, l_pt );
  REQUIRE( l_pt[0] == Approx(0) );

  l_pt[0] = 3;
  edge::linalg::Geom::closestPoint( LINE, l_line, l_pt );
  REQUIRE( l_pt[0] == Approx(1) );

  double l_quad[2][4] = { {0,1,1,0}, {0,0,1,1} };

  l_pt[0] = -2; l_pt[1] = -2;
  edge::linalg::Geom::closestPoint( QUAD4R, l_quad[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0) );
  REQUIRE( l_pt[1] == Approx(0) );

  l_pt[0] = -2; l_pt[1] = 0.5;
  edge::linalg::Geom::closestPoint( QUAD4R, l_quad[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0  ) );
  REQUIRE( l_pt[1] == Approx(0.5) );

  l_pt[0] = -2; l_pt[1] = 2;
  edge::linalg::Geom::closestPoint( QUAD4R, l_quad[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0) );
  REQUIRE( l_pt[1] == Approx(1) );


  double l_hex[3][8] = { {0,1,1,0,0,1,1,0}, {0,0,1,1,0,0,1,1}, {0,0,0,0,1,1,1,1} };

  l_pt[0] = -2; l_pt[1] = 2; l_pt[2] = 3;
  edge::linalg::Geom::closestPoint( HEX8R, l_hex[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0) );
  REQUIRE( l_pt[1] == Approx(1) );
  REQUIRE( l_pt[2] == Approx(1) );

  l_pt[0] = -2; l_pt[1] = 2; l_pt[2] = -3;
  edge::linalg::Geom::closestPoint( HEX8R, l_hex[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0) );
  REQUIRE( l_pt[1] == Approx(1) );
  REQUIRE( l_pt[2] == Approx(0) );

  l_pt[0] = -2; l_pt[1] = 2; l_pt[2] = 0.25;
  edge::linalg::Geom::closestPoint( HEX8R, l_hex[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0)    );
  REQUIRE( l_pt[1] == Approx(1)    );
  REQUIRE( l_pt[2] == Approx(0.25) );

  double l_tria[2][3] = { {0,1,0}, {0,0,1} };

  // TODO: Fails because of a bug in wykobi
  // l_pt[0] = -1; l_pt[1] = -1;
  // edge::linalg::Geom::closestPoint( TRIA3, l_tria[0], l_pt );
  // REQUIRE( l_pt[0] == Approx(0) );
  // REQUIRE( l_pt[1] == Approx(0) );

  l_pt[0] = 0.25; l_pt[1] = -1;
  edge::linalg::Geom::closestPoint( TRIA3, l_tria[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0.25) );
  REQUIRE( l_pt[1] == Approx(0) );

  l_pt[0] = 10; l_pt[1] = 10;
  edge::linalg::Geom::closestPoint( TRIA3, l_tria[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0.5) );
  REQUIRE( l_pt[1] == Approx(0.5) );

  l_pt[0] = 1.0 - 0.33 + 500; l_pt[1] = 0.33 + 500;
  edge::linalg::Geom::closestPoint( TRIA3, l_tria[0], l_pt );
  REQUIRE( l_pt[0] == Approx(1 - 0.33) );
  REQUIRE( l_pt[1] == Approx(0.33) );

  l_pt[0] = 0.25; l_pt[1] = 0.25;
  edge::linalg::Geom::closestPoint( TRIA3, l_tria[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0.25) );
  REQUIRE( l_pt[1] == Approx(0.25) );

  double l_tet[3][4] = { {0,1,0,0}, {0,0,1,0}, {0,0,0,1} };

  l_pt[0] = 0; l_pt[1] = 0; l_pt[2] = 0;
  edge::linalg::Geom::closestPoint( TET4, l_tet[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0) );
  REQUIRE( l_pt[1] == Approx(0) );
  REQUIRE( l_pt[2] == Approx(0) );

  l_pt[0] = 0.5; l_pt[1] = 0; l_pt[2] = 0;
  edge::linalg::Geom::closestPoint( TET4, l_tet[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0.5) );
  REQUIRE( l_pt[1] == Approx(0) );
  REQUIRE( l_pt[2] == Approx(0) );

  l_pt[0] = 0.5; l_pt[1] = 0.5; l_pt[2] = 0;
  edge::linalg::Geom::closestPoint( TET4, l_tet[0], l_pt );
  REQUIRE( l_pt[0] == Approx(0.5) );
  REQUIRE( l_pt[1] == Approx(0.5) );
  REQUIRE( l_pt[2] == Approx(0) );

  l_pt[0] = 1; l_pt[1] = 1; l_pt[2] = 1;
  edge::linalg::Geom::closestPoint( TET4, l_tet[0], l_pt );
  REQUIRE( l_pt[0] == Approx(1.0/3.0) );
  REQUIRE( l_pt[1] == Approx(1.0/3.0) );
  REQUIRE( l_pt[2] == Approx(1.0/3.0) );
}

TEST_CASE( "Inside/outside derivation: Points with respect to quad4 elements", "[geom][nsideQuad4]" ) {
  real_mesh l_veCrds[2][4];
  real_mesh l_pt[2];

  /*
   *   ********* y1=15.3
   *   *       *
   *   *       *
   *   ********* y0=-1.5
   *  x0=-0.25 x1=3.5
   */
  l_veCrds[0][0] = -0.25; l_veCrds[1][0] = -1.5;
  l_veCrds[0][1] =  3.5;  l_veCrds[1][1] = -1.5;
  l_veCrds[0][2] =  3.5;  l_veCrds[1][2] = 15.3;
  l_veCrds[0][3] = -0.25; l_veCrds[1][3] = 15.3;

  l_pt[0] = 0.5; l_pt[1] = 1;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -0.1; l_pt[1] = 14;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 2.3; l_pt[1] = -1;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 2.3; l_pt[1] = 15.3;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = -1.5; l_pt[1] = 15.3;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -2.3; l_pt[1] = -1;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] =  4; l_pt[1] = -1;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 2.3; l_pt[1] = -3;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 2.3; l_pt[1] = 5000;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -1000; l_pt[1] = 10000;
  REQUIRE( edge::linalg::Geom::inside( QUAD4R, (real_mesh*) l_veCrds, l_pt ) == 0 );
}

TEST_CASE( "Inside/outside derivation: Points with respect to tria3 elements", "[geom][insideTria3]" ) {
  real_mesh l_veCrds[2][3];
  real_mesh l_pt[2];

  /*
   *        * 0.3, 7.8
   *       *
   *      *  *
   *     *
   *    *     *
   *   *
   *    *      *
   *      *
   * -1, 2  *   *
   *           *
   *             * 3.3, -4
   */
  l_veCrds[0][0] = -1.0; l_veCrds[1][0] =  2.0;
  l_veCrds[0][1] =  3.3; l_veCrds[1][1] = -4.0;
  l_veCrds[0][2] =  0.3; l_veCrds[1][2] =  7.8;

  l_pt[0] = 0.5; l_pt[1] = 0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 1.0; l_pt[1] = 1.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 1.0; l_pt[1] = 2.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -0.5; l_pt[1] = 1.5;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0; l_pt[1] = 2;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0; l_pt[1] = 6;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0.5; l_pt[1] = 6;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -1; l_pt[1] = 2;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 3.3; l_pt[1] = -4;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 0.5; l_pt[1] = -0.5;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 1.0; l_pt[1] = -2;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 2.0; l_pt[1] = 4.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -2.0; l_pt[1] = 0.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.5; l_pt[1] = 5.5;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 1.5; l_pt[1] = 6.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0; l_pt[1] = 7.0;
  REQUIRE( edge::linalg::Geom::inside( TRIA3, (real_mesh*) l_veCrds, l_pt ) == 0 );
}


TEST_CASE( "Inside/outside derivation: Points with respect to hex8r elements", "[geom][insideHex8r]" ) {
  real_mesh l_veCrds[3][8];
  real_mesh l_pt[3];

  l_veCrds[0][0] = -0.25; l_veCrds[1][0] = -1.5; l_veCrds[2][0] = 3;
  l_veCrds[0][1] =  3.5;  l_veCrds[1][1] = -1.5; l_veCrds[2][1] = 3;
  l_veCrds[0][2] =  3.5;  l_veCrds[1][2] = 15.3; l_veCrds[2][2] = 3;
  l_veCrds[0][3] = -0.25; l_veCrds[1][3] = 15.3; l_veCrds[2][3] = 3;

  l_veCrds[0][4] = -0.25; l_veCrds[1][4] = -1.5; l_veCrds[2][4] = 7;
  l_veCrds[0][5] =  3.5;  l_veCrds[1][5] = -1.5; l_veCrds[2][5] = 7;
  l_veCrds[0][6] =  3.5;  l_veCrds[1][6] = 15.3; l_veCrds[2][6] = 7;
  l_veCrds[0][7] = -0.25; l_veCrds[1][7] = 15.3; l_veCrds[2][7] = 7;

  l_pt[0] = 0.5; l_pt[1] = 6; l_pt[2] = 5;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -0.1; l_pt[1] = 14; l_pt[2] = 6;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 3.0; l_pt[1] = -1; l_pt[2] = 3.01;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -0.25; l_pt[1] = -1.5; l_pt[2] = 6;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 3.5; l_pt[1] = -1; l_pt[2] = 7;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 3.0; l_pt[1] = -1; l_pt[2] = 2.99;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 3.6; l_pt[1] = -1.6; l_pt[2] = 2.99;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0.5; l_pt[1] = 6; l_pt[2] = 7.5;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.1; l_pt[1] = 15.31; l_pt[2] = 6;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.2501; l_pt[1] = 15.31; l_pt[2] = 6;
  REQUIRE( edge::linalg::Geom::inside( HEX8R, (real_mesh*) l_veCrds, l_pt ) == 0 );
}

TEST_CASE( "Inside/outside derivation: Points with respect to tet4 elements", "[geom][insideTet4]" ) {
  real_mesh l_veCrds[3][4];
  real_mesh l_pt[3];

  l_veCrds[0][0] = 0; l_veCrds[1][0] = 0; l_veCrds[2][0] = 0;
  l_veCrds[0][1] = 1; l_veCrds[1][1] = 0; l_veCrds[2][1] = 0;
  l_veCrds[0][2] = 0; l_veCrds[1][2] = 1; l_veCrds[2][2] = 0;
  l_veCrds[0][3] = 0; l_veCrds[1][3] = 0; l_veCrds[2][3] = 1;

  l_pt[0] = 0.25; l_pt[1] = 0.25; l_pt[2] = 0.25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0.1; l_pt[1] = 0.25; l_pt[2] = 0.5;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0.45; l_pt[1] = 0.45; l_pt[2] = 0.001;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = 0.001; l_pt[1] = 0.99; l_pt[2] = 0.0001;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 1 );

  l_pt[0] = -0.0001; l_pt[1] = 0.99; l_pt[2] = 0.00001;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.0001; l_pt[1] = 0.1; l_pt[2] = 0.1;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0.5; l_pt[1] = 0.5; l_pt[2] = 0.5;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0.25; l_pt[1] = -0.25; l_pt[2] = 0.25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0.25; l_pt[1] = 0.25; l_pt[2] = -0.25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.25; l_pt[1] = 0.25; l_pt[2] = -0.25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = -0.25; l_pt[1] = -0.25; l_pt[2] = -0.25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 25; l_pt[1] = -25; l_pt[2] = 25;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 0 );

  l_pt[0] = 0; l_pt[1] = 0; l_pt[2] = 0;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 0; l_pt[1] = 1; l_pt[2] = 0;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 2 );

  l_pt[0] = 0.5; l_pt[1] = 0.5; l_pt[2] = 0;
  REQUIRE( edge::linalg::Geom::inside( TET4, (real_mesh*) l_veCrds, l_pt ) == 2 );
}

TEST_CASE( "Volume: Line 2D", "[GeomVolLine2D]" ) {
  double l_ves[2][2] = { {0.25, 0.3}, {-3.5, 15.12} };

  double l_vol = edge::linalg::Geom::volume( LINE, l_ves[0], 2 );
  REQUIRE( l_vol == Approx( 18.620067132 ) );
}

TEST_CASE( "Volume: Triangle", "[Geom][VolTria]" ) {
  double l_ves[3][3] = { {0.0, 2.0, 2.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0} };
  double l_vol = edge::linalg::Geom::volume( TRIA3, l_ves[0] );
  REQUIRE( l_vol == (2.0 * 3.0) / 2.0 );
}

TEST_CASE( "Diameter: Triangle 2D", "[Geom][DiaTria2D]" ) {
  // construct equilateral triangle
  double l_ves[2][3] = { {0.0, 1.0, 0.5}, {0.0, 0.0, std::sqrt(3.0)/2.0 } };
  double l_inDia = edge::linalg::Geom::inDia( TRIA3, l_ves[0] );
  REQUIRE( l_inDia == Approx(std::sqrt(3)/3.0) );
}

TEST_CASE( "Diameter: Tetrahedron 3D", "[Geom][DiaTet3D]" ) {
  double l_ves[3][4] = { {0.0, 1.0, 0.5, 0.5},
                         {0.0, 0.0, 1.0, 0.5},
                         {0.0, 0.0, 0.0, 1.0} };

  REQUIRE( edge::linalg::Geom::inDia( TET4, l_ves[0] ) == Approx(0.45358449083204) );
}

TEST_CASE( "Geom: Matrix rotating a 3D vector onto another vector", "[Geom][rotMat3D]" ) {
  double l_v0[3];
  double l_v1[3];
  double l_v2[3];
  double l_rt[3];
  double l_rm[3][3];

  l_v0[0] = 1; l_v0[1] = 0; l_v0[2] = 0;
  l_v1[0] = 0; l_v1[1] = 1; l_v1[2] = 0;

  // get rotation matrix
  edge::linalg::GeomTs<3>::rotMat( l_v0, l_v1, l_rm );

  // apply rotation matrix
  edge::linalg::Matrix::matMulB0( 3, 1, 3,
                                  l_rm[0], l_v0, l_rt );

  REQUIRE( l_rt[0] == Approx(l_v1[0]) );
  REQUIRE( l_rt[1] == Approx(l_v1[1]) );
  REQUIRE( l_rt[2] == Approx(l_v1[2]) );

  // apply rotation matrix to "second" coordinate direction
  l_v2[0] = 0; l_v2[1] = 1; l_v2[2] = 0;
  edge::linalg::Matrix::matMulB0( 3, 1, 3,
                                  l_rm[0], l_v2, l_rt );
  REQUIRE( l_rt[0] == Approx(-1.0) );
  REQUIRE( l_rt[1] == Approx( 0.0) );
  REQUIRE( l_rt[2] == Approx( 0.0) );

  // apply rotation matrix to "third" coordinate direction
  l_v2[0] = 0; l_v2[1] = 0; l_v2[2] = 1;
  edge::linalg::Matrix::matMulB0( 3, 1, 3,
                                  l_rm[0], l_v2, l_rt );
  REQUIRE( l_rt[0] == Approx(0.0) );
  REQUIRE( l_rt[1] == Approx(0.0) );
  REQUIRE( l_rt[2] == Approx(1.0) );
}

TEST_CASE( "Geom: Angle between two vectors.", "[Geom][angle]" ) {
  double l_angle = std::numeric_limits< double >::max();

  double l_v1[2] = { 1.0, 0.0 };
  double l_v2[2] = { 0.0, 1.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v1, l_v2 );
  REQUIRE( l_angle == Approx(3.14159265/2) );

  double l_v3[2] = { 10.0, 0.0 };
  double l_v4[2] = {  0.0,  9.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v3, l_v4 );
  REQUIRE( l_angle == Approx(3.14159265/2) );

  double l_v5[2] = { 10.0, 0.0 };
  double l_v6[2] = { 12.0, 0.0} ;
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v5, l_v6 );
  REQUIRE( l_angle == Approx(0) );

  double l_v7[2] = {  10.0, 0.0 };
  double l_v8[2] = { -12.0, 0.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v7, l_v8 );
  REQUIRE( l_angle == Approx(3.14159265) );

  double l_v9[2]  = { 10.0,   0.0 };
  double l_v10[2] = {  0,   -30.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v9, l_v10 );
  REQUIRE( l_angle == Approx(1.5*3.14159265) );

  double l_v11[2] = { 0,    3.0 };
  double l_v12[2] = { 1.0,  1.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v11, l_v12 );
  REQUIRE( l_angle == Approx(1.75*3.14159265) );

  double l_v13[2] = { -1.0, -1.0  };
  double l_v14[2] = {  1.0, -1.0 };
  l_angle = edge::linalg::GeomTs< 2 >::angle( l_v13, l_v14 );
  REQUIRE( l_angle == Approx(0.5*3.14159265) );
}
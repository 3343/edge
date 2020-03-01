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
 * Tests the geometry computations on tria3 elements.
 **/
#include <catch.hpp>
#define private public
#include "Tria3.h"
#undef private

TEST_CASE( "Tests the volume computation for 3-node triangles.", "[volume][tria3]" ) {
  double l_ves1[3][3] = { {0.0, 2.0, 2.0}, {0.0, 0.0, 3.0}, {0.0, 0.0, 0.0} };
  double l_vol1 = edge_v::geom::Tria3::area( l_ves1 );
  REQUIRE( l_vol1 == Approx((2.0 * 3.0) / 2.0) );


  double l_ves2[3][3] = { {1.2, 3.0, 4.0}, {1.0, 5.0, 2.0}, {2.0, 4.0, 3.0} };
  double l_vol2 = edge_v::geom::Tria3::area( l_ves2 );
  REQUIRE( l_vol2 == Approx(1.27279) );
}

TEST_CASE( "Tests the tangents computation for 3-node triangles.", "[tangents][tria3]" ) {
  double l_ves0[3][3] = { {3.0, 2.0, 2.0}, {5.0, 4.0, 2.0}, {3.0, 4.0, 2.0} };
  double l_np0[3] = { 4.5, 2.5, -3.0 };

  double l_tangents[2][3] = {{0}};

  edge_v::geom::Tria3::tangents( l_ves0,
                                 l_np0,
                                 l_tangents );

  REQUIRE( l_tangents[0][0] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[0][1] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[0][2] == Approx(  0              ) );

  REQUIRE( l_tangents[1][0] == Approx( -std::sqrt(0.5) ) );
  REQUIRE( l_tangents[1][1] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[1][2] == Approx(  0              ) );

  double l_np1[3] = { 4.5, 2.5, 3.0 };

  edge_v::geom::Tria3::tangents( l_ves0,
                                 l_np1,
                                 l_tangents );

  REQUIRE( l_tangents[0][0] == Approx( -std::sqrt(0.5) ) );
  REQUIRE( l_tangents[0][1] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[0][2] == Approx(  0              ) );

  REQUIRE( l_tangents[1][0] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[1][1] == Approx(  std::sqrt(0.5) ) );
  REQUIRE( l_tangents[1][2] == Approx(  0              ) );
}

TEST_CASE( "Tests the normal computation for 3-node triangles.", "[normal][tria3]" ) {
  double l_ves0[3][3] = { {3.0, 2.0, 2.0}, {5.0, 4.0, 2.0}, {3.0, 4.0, 2.0} };
  double l_np0[3] = { 4.5, 2.5, -3.0 };

  double l_normal[3] = {0};

  edge_v::geom::Tria3::normal( l_ves0,
                               l_np0,
                               l_normal );

  REQUIRE( l_normal[0] == Approx(0.0) );
  REQUIRE( l_normal[1] == Approx(0.0) );
  REQUIRE( l_normal[2] == Approx(1.0) );

  double l_np1[3] = { 4.5, 2.5, 5.0 };

  edge_v::geom::Tria3::normal( l_ves0,
                               l_np1,
                               l_normal );

  REQUIRE( l_normal[0] == Approx( 0.0) );
  REQUIRE( l_normal[1] == Approx( 0.0) );
  REQUIRE( l_normal[2] == Approx(-1.0) );
}

TEST_CASE( "Tests the incircle diameter computation for 3-node triangles.", "[inDiameter][tria3]" ) { 
  // construct equilateral triangle
  double l_ves[3][3] = { {0.0, 0.0, 0.0} , {1.0, 0.0, 0.0}, {0.5, std::sqrt(3.0)/2.0, 0.0 } };

  // compute diameter
  double l_dia = edge_v::geom::Tria3::inDiameter( l_ves );

  // check result
  REQUIRE( l_dia == Approx( std::sqrt(3)/3.0 ) );
}


TEST_CASE( "Tests the vertex and face reordering for triangles.", "[normVesFas][tria3]" ) {
  /*
   * counterclockwise triangle:
   *
   *   3
   *   *  *
   *   *    *
   *   *      *
   *   *        *
   *   1 ********* 2
   *
   * faces (ascending): (1, 2) < (1, 3) < (2, 3)
   * faces (counter-clockwise): (1, 2), (2, 3), (1, 3)
   */
  double l_ves0[3][3] = { {0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.5, 1.0, 0.0} };

  edge_v::t_idx l_elVe0[3]   = { 1, 2, 3 };
  edge_v::t_idx l_elFa0[3]   = { 4, 5, 6 };
  edge_v::t_idx l_elFaEl0[3] = { 7, 8, 9 };

  edge_v::geom::Tria3::normVesFas( l_ves0,
                                   l_elVe0,
                                   l_elFa0,
                                   l_elFaEl0 );

  REQUIRE( l_elVe0[0] == 1 );
  REQUIRE( l_elVe0[1] == 2 );
  REQUIRE( l_elVe0[2] == 3 );

  REQUIRE( l_elFa0[0] == 4 );
  REQUIRE( l_elFa0[1] == 6 );
  REQUIRE( l_elFa0[2] == 5 );

  REQUIRE( l_elFaEl0[0] == 7 );
  REQUIRE( l_elFaEl0[1] == 9 );
  REQUIRE( l_elFaEl0[2] == 8 );

  /*
   * counterclockwise triangle:
   *
   *   3                  2
   *   *  *               * *
   *   *    *         ->  *   *
   *   *      *           *     *
   *   *        *         *       *
   *   2 ********* 1      3 ******** 1
   *
   * faces (ascending): (1, 2) < (1, 3) < (2, 3)
   * faces (counter-clockwise): (1, 2), (2, 3), (1, 3)
   *
   * counterclockwise ordering exchanges vertices 2 and 3:
   * faces ordering: (1, 3), (2, 3), (1, 2)
   */
  double l_ves1[3][3] = { {1.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0},
                          {0.5, 1.0, 0.0} };


  edge_v::t_idx l_elVe1[3]   = { 1, 2, 3 };
  edge_v::t_idx l_elFa1[3]   = { 4, 5, 6 };
  edge_v::t_idx l_elFaEl1[3] = { 7, 8, 9 };

  edge_v::geom::Tria3::normVesFas( l_ves1,
                                   l_elVe1,
                                   l_elFa1,
                                   l_elFaEl1 );

  REQUIRE( l_elVe1[0] == 1 );
  REQUIRE( l_elVe1[1] == 3 );
  REQUIRE( l_elVe1[2] == 2 );

  REQUIRE( l_elFa1[0] == 5 );
  REQUIRE( l_elFa1[1] == 6 );
  REQUIRE( l_elFa1[2] == 4 );

  REQUIRE( l_elFaEl1[0] == 8 );
  REQUIRE( l_elFaEl1[1] == 9 );
  REQUIRE( l_elFaEl1[2] == 7 );
}
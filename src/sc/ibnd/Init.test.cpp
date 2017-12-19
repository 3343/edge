/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Unit tests for initialization of internal boundary data.
 **/
#include <catch.hpp>
#define private public
#include "Init.hpp"
#undef private

TEST_CASE( "Ibnd init: Connectivity.", "[ibndInit][conn]") {
  /*
   * Our mini-example. Unit tests only / not consistent.
   *
   * 1: limited
   * 2: internal boundary
   *
   *   fa | sp | el0 | el1 | ibndId |
   *    0 |  0 |   0 |   3 |      x |
   *    1 |  0 |   3 |   0 |      x |
   *    2 |  2 |   2 |   1 |      0 |
   *    3 |  2 |   7 |   5 |      1 |
   *    4 |  2 |   4 |   6 |      2 |
   *    5 |  2 |   1 |   2 |      3 |
   *    6 |  0 |   9 |   8 |      x |
   *    7 |  2 |   4 |   3 |      4 |
   *
   *   el  | sp | fa0 | fa1 | fa2 | limId | ibndId |
   *    0  |  0 |   0 |   1 |   6 |     x |      x |
   *    1  |  1 |   1 |   0 |   6 |     0 |      x |
   *    2  |  3 |   2 |   7 |   5 |     1 |      0 |
   *    3  |  1 |   0 |   0 |   1 |     2 |      x |
   *    4  |  3 |   4 |   0 |   7 |     3 |      1 |
   *    5  |  3 |   2 |   4 |   3 |     4 |      2 |
   *    6  |  0 |   6 |   0 |   1 |     x |      x |
   *    7  |  1 |   0 |   1 |   6 |     5 |      x |
   *    8  |  0 |   1 |   0 |   6 |     x |      x |
   *    9  |  0 |   6 |   1 |   0 |     x |      x |
   */

  // set up example
  typedef struct { unsigned short spType; } t_enChars;

  t_enChars l_charsFa1[8]  = { {0}, {0}, {2}, {2}, {2},
                               {2}, {0}, {2} };
  t_enChars l_charsEl1[10] = { {0}, {1}, {3}, {1}, {3},
                               {3}, {0}, {1}, {0}, {0} };

  int l_faEl1[8][2] = { {0, 3}, {3, 0}, {2, 1}, {7, 5}, {4, 6},
                        {1, 2}, {9, 8}, {4, 3} };

  int l_elFa1[10][3] = { {0, 1, 6}, {1, 0, 6}, {2, 7, 5}, {0, 0, 1}, {4, 0, 7},
                         {2, 4, 3}, {6, 0, 1}, {0, 1, 6}, {1, 0, 6}, {6, 1, 0 } };

  // "allocate" memore for ibnd connectivity
  int l_bfLe1[5][2];
  int l_bfBe1[5][2] = { {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1} };
  int l_beBf1[3][3];
  int l_beEl1[3];

  edge::sc::ibnd::t_connect< int, TRIA3 > l_conn1;
  l_conn1.bfLe = l_bfLe1;
  l_conn1.bfBe = l_bfBe1;
  l_conn1.beBf = l_beBf1;
  l_conn1.beEl = l_beEl1;

  // call the initialization
  edge::sc::ibnd::Init::connect( 8,
                                 10,
                                 2,
                                 1,
                                 l_faEl1,
                                 l_elFa1,
                                 l_charsFa1,
                                 l_charsEl1,
                                 l_conn1 );

  // check link between ibnd faces and ibnd elements, remark: max() would be invalid in application
  REQUIRE( l_conn1.bfBe[0][0] == 0 );
  REQUIRE( l_conn1.bfBe[0][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_conn1.bfBe[1][0] == std::numeric_limits< int >::max() );
  REQUIRE( l_conn1.bfBe[1][1] == 2 );

  REQUIRE( l_conn1.bfBe[2][0] == 1 );
  REQUIRE( l_conn1.bfBe[2][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_conn1.bfBe[3][0] == std::numeric_limits< int >::max() );
  REQUIRE( l_conn1.bfBe[3][1] == 0 );

  REQUIRE( l_conn1.bfBe[4][0] == 1 );
  REQUIRE( l_conn1.bfBe[4][1] == std::numeric_limits< int >::max() );

  // check link between ibnd elements and ibnd faces
  REQUIRE( l_conn1.beBf[0][0] == 0 );
  REQUIRE( l_conn1.beBf[0][1] == 4 );
  REQUIRE( l_conn1.beBf[0][2] == 3 );

  REQUIRE( l_conn1.beBf[1][0] == 2 );
  REQUIRE( l_conn1.beBf[1][1] == std::numeric_limits< int >::max() );
  REQUIRE( l_conn1.beBf[1][2] == 4 );

  REQUIRE( l_conn1.beBf[2][0] == 0 );
  REQUIRE( l_conn1.beBf[2][1] == 2 );
  REQUIRE( l_conn1.beBf[2][2] == 1 );

  // check link between ibnd faces and limited faces, remark: max() would be invalid in application
  REQUIRE( l_conn1.bfLe[0][0] == 1 );
  REQUIRE( l_conn1.bfLe[0][1] == 0 );

  REQUIRE( l_conn1.bfLe[1][0] == 5 );
  REQUIRE( l_conn1.bfLe[1][1] == 4 );

  REQUIRE( l_conn1.bfLe[2][0] == 3 );
  REQUIRE( l_conn1.bfLe[2][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_conn1.bfLe[3][0] == 0 );
  REQUIRE( l_conn1.bfLe[3][1] == 1 );

  REQUIRE( l_conn1.bfLe[4][0] == 3 );
  REQUIRE( l_conn1.bfLe[4][1] == 2 );

  // check link between internal boundary elements and dense elements
  REQUIRE( l_conn1.beEl[0] == 2 );
  REQUIRE( l_conn1.beEl[1] == 4 );
  REQUIRE( l_conn1.beEl[2] == 5 );
}

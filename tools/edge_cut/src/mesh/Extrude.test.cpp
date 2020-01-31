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
 * Tests the 2D surface mesh extrusion.
 **/

#include <catch.hpp>
#define private public
#include "Extrude.h"
#undef private
#include <regex>

TEST_CASE( "Tests the extrusion of surface meshes.", "[extrude][Extrude]" ) {
  // coordinates on a regular x-y grid
  double l_veCrds0[12][3] = {  { 1, 5, 2},
                               { 3, 5, 1},
                               { 6, 5, 3},
                               {10, 5, 3},

                               { 1, 8, 3},
                               { 3, 8, 2},
                               { 6, 8, 4},
                               {10, 8, 5},

                               { 1, 9, 2},
                               { 3, 9, 4},
                               { 6, 9, 2},
                               {10, 9, 1} };

  // reorder version above
  double l_veCrds1[12][3] = {  { 1, 5, 2},
                               { 3, 5, 1},
                               {10, 8, 5},
                               { 6, 5, 3},
                               { 3, 8, 2},
                               {10, 9, 1},
                               {10, 5, 3},
                               { 1, 8, 3},
                               { 6, 8, 4},
                               { 6, 9, 2},
                               { 1, 9, 2},
                               { 3, 9, 4} };

  edge_cut::mesh::Extrude l_ext0( 12,
                                  l_veCrds0,
                                  -4 );

  REQUIRE( l_ext0.m_veMinSurf[0] == 1 );
  REQUIRE( l_ext0.m_veMinSurf[1] == 5 );
  REQUIRE( l_ext0.m_veMinSurf[2] == 1 );

  REQUIRE( l_ext0.m_veMaxSurf[0] == 10 );
  REQUIRE( l_ext0.m_veMaxSurf[1] ==  9 );
  REQUIRE( l_ext0.m_veMaxSurf[2] ==  5 );

  REQUIRE( l_ext0.m_nx == 4 );
  REQUIRE( l_ext0.m_ny == 3 );

  REQUIRE( l_ext0.m_dMin[0] == 2 );
  REQUIRE( l_ext0.m_dMin[1] == 1 );

  // test output writer
  edge_cut::mesh::Extrude l_ext1( 12,
                                  l_veCrds1,
                                  -5 );

  std::stringstream l_left, l_right, l_front, l_back, l_bottom, l_top;

  l_ext1.writeOff( l_left,
                   l_right,
                   l_front,
                   l_back,
                   l_bottom,
                   l_top );

  std::string l_refLeft = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
36 44 0
1 5 2
1 8 3
1 9 2
1 5 1.36364
1 8 2.27273
1 9 1.36364
1 5 0.727273
1 8 1.54545
1 9 0.727273
1 5 0.0909091
1 8 0.818182
1 9 0.0909091
1 5 -0.545455
1 8 0.0909091
1 9 -0.545455
1 5 -1.18182
1 8 -0.636364
1 9 -1.18182
1 5 -1.81818
1 8 -1.36364
1 9 -1.81818
1 5 -2.45455
1 8 -2.09091
1 9 -2.45455
1 5 -3.09091
1 8 -2.81818
1 9 -3.09091
1 5 -3.72727
1 8 -3.54545
1 9 -3.72727
1 5 -4.36364
1 8 -4.27273
1 9 -4.36364
1 5 -5
1 8 -5
1 9 -5
3 0 1 4
3 0 4 3
3 3 4 7
3 3 7 6
3 6 7 10
3 6 10 9
3 9 10 13
3 9 13 12
3 12 13 16
3 12 16 15
3 15 16 19
3 15 19 18
3 18 19 22
3 18 22 21
3 21 22 25
3 21 25 24
3 24 25 28
3 24 28 27
3 27 28 31
3 27 31 30
3 30 31 34
3 30 34 33
3 1 2 5
3 1 5 4
3 4 5 8
3 4 8 7
3 7 8 11
3 7 11 10
3 10 11 14
3 10 14 13
3 13 14 17
3 13 17 16
3 16 17 20
3 16 20 19
3 19 20 23
3 19 23 22
3 22 23 26
3 22 26 25
3 25 26 29
3 25 29 28
3 28 29 32
3 28 32 31
3 31 32 35
3 31 35 34
)V0G0N";

  std::string l_refRight = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
36 44 0
10 5 3
10 8 5
10 9 1
10 5 2.27273
10 8 4.09091
10 9 0.454545
10 5 1.54545
10 8 3.18182
10 9 -0.0909091
10 5 0.818182
10 8 2.27273
10 9 -0.636364
10 5 0.0909091
10 8 1.36364
10 9 -1.18182
10 5 -0.636364
10 8 0.454545
10 9 -1.72727
10 5 -1.36364
10 8 -0.454545
10 9 -2.27273
10 5 -2.09091
10 8 -1.36364
10 9 -2.81818
10 5 -2.81818
10 8 -2.27273
10 9 -3.36364
10 5 -3.54545
10 8 -3.18182
10 9 -3.90909
10 5 -4.27273
10 8 -4.09091
10 9 -4.45455
10 5 -5
10 8 -5
10 9 -5
3 0 1 4
3 0 4 3
3 3 4 7
3 3 7 6
3 6 7 10
3 6 10 9
3 9 10 13
3 9 13 12
3 12 13 16
3 12 16 15
3 15 16 19
3 15 19 18
3 18 19 22
3 18 22 21
3 21 22 25
3 21 25 24
3 24 25 28
3 24 28 27
3 27 28 31
3 27 31 30
3 30 31 34
3 30 34 33
3 1 2 5
3 1 5 4
3 4 5 8
3 4 8 7
3 7 8 11
3 7 11 10
3 10 11 14
3 10 14 13
3 13 14 17
3 13 17 16
3 16 17 20
3 16 20 19
3 19 20 23
3 19 23 22
3 22 23 26
3 22 26 25
3 25 26 29
3 25 29 28
3 28 29 32
3 28 32 31
3 31 32 35
3 31 35 34
)V0G0N";

  std::string l_refFront = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
48 66 0
1 5 2
3 5 1
6 5 3
10 5 3
1 5 1.36364
3 5 0.454545
6 5 2.27273
10 5 2.27273
1 5 0.727273
3 5 -0.0909091
6 5 1.54545
10 5 1.54545
1 5 0.0909091
3 5 -0.636364
6 5 0.818182
10 5 0.818182
1 5 -0.545455
3 5 -1.18182
6 5 0.0909091
10 5 0.0909091
1 5 -1.18182
3 5 -1.72727
6 5 -0.636364
10 5 -0.636364
1 5 -1.81818
3 5 -2.27273
6 5 -1.36364
10 5 -1.36364
1 5 -2.45455
3 5 -2.81818
6 5 -2.09091
10 5 -2.09091
1 5 -3.09091
3 5 -3.36364
6 5 -2.81818
10 5 -2.81818
1 5 -3.72727
3 5 -3.90909
6 5 -3.54545
10 5 -3.54545
1 5 -4.36364
3 5 -4.45455
6 5 -4.27273
10 5 -4.27273
1 5 -5
3 5 -5
6 5 -5
10 5 -5
3 0 1 5
3 0 5 4
3 4 5 9
3 4 9 8
3 8 9 13
3 8 13 12
3 12 13 17
3 12 17 16
3 16 17 21
3 16 21 20
3 20 21 25
3 20 25 24
3 24 25 29
3 24 29 28
3 28 29 33
3 28 33 32
3 32 33 37
3 32 37 36
3 36 37 41
3 36 41 40
3 40 41 45
3 40 45 44
3 1 2 6
3 1 6 5
3 5 6 10
3 5 10 9
3 9 10 14
3 9 14 13
3 13 14 18
3 13 18 17
3 17 18 22
3 17 22 21
3 21 22 26
3 21 26 25
3 25 26 30
3 25 30 29
3 29 30 34
3 29 34 33
3 33 34 38
3 33 38 37
3 37 38 42
3 37 42 41
3 41 42 46
3 41 46 45
3 2 3 7
3 2 7 6
3 6 7 11
3 6 11 10
3 10 11 15
3 10 15 14
3 14 15 19
3 14 19 18
3 18 19 23
3 18 23 22
3 22 23 27
3 22 27 26
3 26 27 31
3 26 31 30
3 30 31 35
3 30 35 34
3 34 35 39
3 34 39 38
3 38 39 43
3 38 43 42
3 42 43 47
3 42 47 46
)V0G0N";

  std::string l_refBack = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
48 66 0
1 9 2
3 9 4
6 9 2
10 9 1
1 9 1.36364
3 9 3.18182
6 9 1.36364
10 9 0.454545
1 9 0.727273
3 9 2.36364
6 9 0.727273
10 9 -0.0909091
1 9 0.0909091
3 9 1.54545
6 9 0.0909091
10 9 -0.636364
1 9 -0.545455
3 9 0.727273
6 9 -0.545455
10 9 -1.18182
1 9 -1.18182
3 9 -0.0909091
6 9 -1.18182
10 9 -1.72727
1 9 -1.81818
3 9 -0.909091
6 9 -1.81818
10 9 -2.27273
1 9 -2.45455
3 9 -1.72727
6 9 -2.45455
10 9 -2.81818
1 9 -3.09091
3 9 -2.54545
6 9 -3.09091
10 9 -3.36364
1 9 -3.72727
3 9 -3.36364
6 9 -3.72727
10 9 -3.90909
1 9 -4.36364
3 9 -4.18182
6 9 -4.36364
10 9 -4.45455
1 9 -5
3 9 -5
6 9 -5
10 9 -5
3 0 1 5
3 0 5 4
3 4 5 9
3 4 9 8
3 8 9 13
3 8 13 12
3 12 13 17
3 12 17 16
3 16 17 21
3 16 21 20
3 20 21 25
3 20 25 24
3 24 25 29
3 24 29 28
3 28 29 33
3 28 33 32
3 32 33 37
3 32 37 36
3 36 37 41
3 36 41 40
3 40 41 45
3 40 45 44
3 1 2 6
3 1 6 5
3 5 6 10
3 5 10 9
3 9 10 14
3 9 14 13
3 13 14 18
3 13 18 17
3 17 18 22
3 17 22 21
3 21 22 26
3 21 26 25
3 25 26 30
3 25 30 29
3 29 30 34
3 29 34 33
3 33 34 38
3 33 38 37
3 37 38 42
3 37 42 41
3 41 42 46
3 41 46 45
3 2 3 7
3 2 7 6
3 6 7 11
3 6 11 10
3 10 11 15
3 10 15 14
3 14 15 19
3 14 19 18
3 18 19 23
3 18 23 22
3 22 23 27
3 22 27 26
3 26 27 31
3 26 31 30
3 30 31 35
3 30 35 34
3 34 35 39
3 34 39 38
3 38 39 43
3 38 43 42
3 42 43 47
3 42 47 46
)V0G0N";

  std::string l_refBottom = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
12 12 0
1 5 -5
1 8 -5
1 9 -5
3 5 -5
3 8 -5
3 9 -5
6 5 -5
6 8 -5
6 9 -5
10 5 -5
10 8 -5
10 9 -5
3 0 1 4
3 0 4 3
3 3 4 7
3 3 7 6
3 6 7 10
3 6 10 9
3 1 2 5
3 1 5 4
3 4 5 8
3 4 8 7
3 7 8 11
3 7 11 10
)V0G0N";

  std::string l_refTop = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
12 12 0
1 5 2
1 8 3
1 9 2
3 5 1
3 8 2
3 9 4
6 5 3
6 8 4
6 9 2
10 5 3
10 8 5
10 9 1
3 0 1 4
3 0 4 3
3 3 4 7
3 3 7 6
3 6 7 10
3 6 10 9
3 1 2 5
3 1 5 4
3 4 5 8
3 4 8 7
3 7 8 11
3 7 11 10
)V0G0N";

  // replace version strings with compile time constant
  l_refLeft   = std::regex_replace( l_refLeft,   std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );
  l_refRight  = std::regex_replace( l_refRight,  std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );
  l_refFront  = std::regex_replace( l_refFront,  std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );
  l_refBack   = std::regex_replace( l_refBack,   std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );
  l_refBottom = std::regex_replace( l_refBottom, std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );
  l_refTop    = std::regex_replace( l_refTop,    std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );


  REQUIRE( l_left.str()   == l_refLeft   );
  REQUIRE( l_right.str()  == l_refRight  );
  REQUIRE( l_front.str()  == l_refFront  );
  REQUIRE( l_back.str()   == l_refBack   );
  REQUIRE( l_bottom.str() == l_refBottom );
  REQUIRE( l_top.str()    == l_refTop    );
}
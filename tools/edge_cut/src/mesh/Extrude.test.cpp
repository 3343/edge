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

  REQUIRE( l_ext0.m_nx[0] == 4 );
  REQUIRE( l_ext0.m_ny[0] == 3 );

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
1 9 2
1 8 3
1 5 2
1 9 1.36364
1 8 2.27273
1 5 1.36364
1 9 0.727273
1 8 1.54545
1 5 0.727273
1 9 0.0909091
1 8 0.818182
1 5 0.0909091
1 9 -0.545455
1 8 0.0909091
1 5 -0.545455
1 9 -1.18182
1 8 -0.636364
1 5 -1.18182
1 9 -1.81818
1 8 -1.36364
1 5 -1.81818
1 9 -2.45455
1 8 -2.09091
1 5 -2.45455
1 9 -3.09091
1 8 -2.81818
1 5 -3.09091
1 9 -3.72727
1 8 -3.54545
1 5 -3.72727
1 9 -4.36364
1 8 -4.27273
1 5 -4.36364
1 9 -5
1 8 -5
1 5 -5
3 0 3 4
3 0 4 1
3 1 4 5
3 1 5 2
3 3 6 7
3 3 7 4
3 4 7 8
3 4 8 5
3 6 9 10
3 6 10 7
3 7 10 11
3 7 11 8
3 9 12 13
3 9 13 10
3 10 13 14
3 10 14 11
3 12 15 16
3 12 16 13
3 13 16 17
3 13 17 14
3 15 18 19
3 15 19 16
3 16 19 20
3 16 20 17
3 18 21 22
3 18 22 19
3 19 22 23
3 19 23 20
3 21 24 25
3 21 25 22
3 22 25 26
3 22 26 23
3 24 27 28
3 24 28 25
3 25 28 29
3 25 29 26
3 27 30 31
3 27 31 28
3 28 31 32
3 28 32 29
3 30 33 34
3 30 34 31
3 31 34 35
3 31 35 32
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
3 0 3 4
3 0 4 1
3 1 4 5
3 1 5 2
3 3 6 7
3 3 7 4
3 4 7 8
3 4 8 5
3 6 9 10
3 6 10 7
3 7 10 11
3 7 11 8
3 9 12 13
3 9 13 10
3 10 13 14
3 10 14 11
3 12 15 16
3 12 16 13
3 13 16 17
3 13 17 14
3 15 18 19
3 15 19 16
3 16 19 20
3 16 20 17
3 18 21 22
3 18 22 19
3 19 22 23
3 19 23 20
3 21 24 25
3 21 25 22
3 22 25 26
3 22 26 23
3 24 27 28
3 24 28 25
3 25 28 29
3 25 29 26
3 27 30 31
3 27 31 28
3 28 31 32
3 28 32 29
3 30 33 34
3 30 34 31
3 31 34 35
3 31 35 32
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
3 0 4 5
3 0 5 1
3 1 5 6
3 1 6 2
3 2 6 7
3 2 7 3
3 4 8 9
3 4 9 5
3 5 9 10
3 5 10 6
3 6 10 11
3 6 11 7
3 8 12 13
3 8 13 9
3 9 13 14
3 9 14 10
3 10 14 15
3 10 15 11
3 12 16 17
3 12 17 13
3 13 17 18
3 13 18 14
3 14 18 19
3 14 19 15
3 16 20 21
3 16 21 17
3 17 21 22
3 17 22 18
3 18 22 23
3 18 23 19
3 20 24 25
3 20 25 21
3 21 25 26
3 21 26 22
3 22 26 27
3 22 27 23
3 24 28 29
3 24 29 25
3 25 29 30
3 25 30 26
3 26 30 31
3 26 31 27
3 28 32 33
3 28 33 29
3 29 33 34
3 29 34 30
3 30 34 35
3 30 35 31
3 32 36 37
3 32 37 33
3 33 37 38
3 33 38 34
3 34 38 39
3 34 39 35
3 36 40 41
3 36 41 37
3 37 41 42
3 37 42 38
3 38 42 43
3 38 43 39
3 40 44 45
3 40 45 41
3 41 45 46
3 41 46 42
3 42 46 47
3 42 47 43
)V0G0N";

  std::string l_refBack = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
48 66 0
10 9 1
6 9 2
3 9 4
1 9 2
10 9 0.454545
6 9 1.36364
3 9 3.18182
1 9 1.36364
10 9 -0.0909091
6 9 0.727273
3 9 2.36364
1 9 0.727273
10 9 -0.636364
6 9 0.0909091
3 9 1.54545
1 9 0.0909091
10 9 -1.18182
6 9 -0.545455
3 9 0.727273
1 9 -0.545455
10 9 -1.72727
6 9 -1.18182
3 9 -0.0909091
1 9 -1.18182
10 9 -2.27273
6 9 -1.81818
3 9 -0.909091
1 9 -1.81818
10 9 -2.81818
6 9 -2.45455
3 9 -1.72727
1 9 -2.45455
10 9 -3.36364
6 9 -3.09091
3 9 -2.54545
1 9 -3.09091
10 9 -3.90909
6 9 -3.72727
3 9 -3.36364
1 9 -3.72727
10 9 -4.45455
6 9 -4.36364
3 9 -4.18182
1 9 -4.36364
10 9 -5
6 9 -5
3 9 -5
1 9 -5
3 0 4 5
3 0 5 1
3 1 5 6
3 1 6 2
3 2 6 7
3 2 7 3
3 4 8 9
3 4 9 5
3 5 9 10
3 5 10 6
3 6 10 11
3 6 11 7
3 8 12 13
3 8 13 9
3 9 13 14
3 9 14 10
3 10 14 15
3 10 15 11
3 12 16 17
3 12 17 13
3 13 17 18
3 13 18 14
3 14 18 19
3 14 19 15
3 16 20 21
3 16 21 17
3 17 21 22
3 17 22 18
3 18 22 23
3 18 23 19
3 20 24 25
3 20 25 21
3 21 25 26
3 21 26 22
3 22 26 27
3 22 27 23
3 24 28 29
3 24 29 25
3 25 29 30
3 25 30 26
3 26 30 31
3 26 31 27
3 28 32 33
3 28 33 29
3 29 33 34
3 29 34 30
3 30 34 35
3 30 35 31
3 32 36 37
3 32 37 33
3 33 37 38
3 33 38 34
3 34 38 39
3 34 39 35
3 36 40 41
3 36 41 37
3 37 41 42
3 37 42 38
3 38 42 43
3 38 43 39
3 40 44 45
3 40 45 41
3 41 45 46
3 41 46 42
3 42 46 47
3 42 47 43
)V0G0N";

  std::string l_refBottom = R"V0G0N(OFF
#
# EDGEcut generated surface mesh (version: TEMPLATE_EDGE_VERSION)
# EDGEcut is available from https://dial3343.org
#
12 12 0
1 5 -5
3 5 -5
6 5 -5
10 5 -5
1 8 -5
3 8 -5
6 8 -5
10 8 -5
1 9 -5
3 9 -5
6 9 -5
10 9 -5
3 0 4 5
3 0 5 1
3 1 5 6
3 1 6 2
3 2 6 7
3 2 7 3
3 4 8 9
3 4 9 5
3 5 9 10
3 5 10 6
3 6 10 11
3 6 11 7
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
3 0 3 4
3 0 4 1
3 1 4 5
3 1 5 2
3 3 6 7
3 3 7 4
3 4 7 8
3 4 8 5
3 6 9 10
3 6 10 7
3 7 10 11
3 7 11 8
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

TEST_CASE( "Tests splitting quads bottom left to top right.", "[extrude][writeTriasBlTr]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasBlTr( 5,
                                           7,
                                           1,
                                           l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 12 13
3 5 13 6
)V0G0N" );
}

TEST_CASE( "Tests splitting quads bottom right to top left.", "[extrude][writeTriasBrTl]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasBrTl( 5,
                                           7,
                                           1,
                                           l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 12 6
3 12 13 6
)V0G0N" );
}

TEST_CASE( "Tests splitting quads at the center right.", "[extrude][writeTriasCr]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasCr( 5,
                                         7,
                                         1,
                                         l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 19 20
3 5 20 7
3 20 21 7
)V0G0N" );
}

TEST_CASE( "Tests splitting quads at the center left.", "[extrude][writeTriasCl]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasCl( 5,
                                         7,
                                         1,
                                         l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 19 6
3 19 21 6
3 6 21 7
)V0G0N" );
}

TEST_CASE( "Tests splitting quads at the center top.", "[extrude][writeTriasCt]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasCt( 5,
                                         7,
                                         1,
                                         l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 14 7
3 5 19 14
3 19 21 14
)V0G0N" );
}

TEST_CASE( "Tests splitting quads at the center bottom.", "[extrude][writeTriasCb]" ) {
  std::stringstream l_stream;
  edge_cut::mesh::Extrude::writeTriasCb( 5,
                                         7,
                                         1,
                                         l_stream );
  REQUIRE( l_stream.str() == R"V0G0N(3 5 12 7
3 12 21 7
3 12 19 21
)V0G0N" );
}
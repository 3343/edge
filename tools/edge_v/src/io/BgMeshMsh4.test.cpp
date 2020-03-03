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
 * Tests the background mesh writer in msh4 format.
 **/
#include <catch.hpp>
#define private public
#include "BgMeshMsh4.h"
#undef private
#include <regex>

TEST_CASE( "Tests the bg mesh writer in msh4 format.", "[bgMeshMsh4][write]" ) {
  double l_veCrds[6][3] = { {0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {2, 0, 0},
                            {2, 1, 0} };

  edge_v::t_idx l_elVe[4][3] = { {0, 1, 2},
                                 {0, 2, 3},
                                 {1, 5, 2},
                                 {1, 4, 5} };

  float l_targetLengths[6] = { 0.1, 0.15, 0.2, 0.4, 0.8, 1 };

  std::stringstream l_stream;
  edge_v::io::BgMeshMsh4::write( edge_v::TRIA3,
                                 6,
                                 4,
                                 l_elVe[0],
                                 l_veCrds,
                                 l_targetLengths,
                                 l_stream );
  std::string l_ref = R"V0G0N($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 6 0 5
2 1 0 6
0
1
2
3
4
5
0 0 0
1 0 0
1 1 0
0 1 0
2 0 0
2 1 0
$EndNodes
$Elements
1 4 0 3
2 1 2 4
0 0 1 2
1 0 2 3
2 1 5 2
3 1 4 5
$EndElements
$NodeData
1
"EDGE-V Background Mesh (git-version: TEMPLATE_EDGE_VERSION)"
1
0.0
3
0
1
6
0 0.1
1 0.15
2 0.2
3 0.4
4 0.8
5 1
$EndNodeData
)V0G0N";

  // replace version string with compile time constant
  l_ref = std::regex_replace( l_ref, std::regex("TEMPLATE_EDGE_VERSION"), PP_EDGE_VERSION );

  REQUIRE( l_stream.str().size() == l_ref.size() );
  REQUIRE( l_stream.str() == l_ref );
}
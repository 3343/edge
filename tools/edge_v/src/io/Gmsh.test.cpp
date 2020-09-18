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
 * Tests the Gmsh interface.
 **/
#include <catch.hpp>
#define private public
#include "Gmsh.h"
#undef private

namespace edge_v {
  namespace test {
    extern std::string g_files;
  }
}

TEST_CASE( "Tests the derivation of entity types.", "[gmsh][entityTypes]" ) {
  edge_v::io::Gmsh l_gmsh;
  int l_gmshPoint  = l_gmsh.getGmshType( edge_v::POINT );
  int l_gmshLine   = l_gmsh.getGmshType( edge_v::LINE );
  int l_gmshQuad4r = l_gmsh.getGmshType( edge_v::QUAD4R );
  int l_gmshTria3  = l_gmsh.getGmshType( edge_v::TRIA3 );
  int l_gmshHex8r  = l_gmsh.getGmshType( edge_v::HEX8R );
  int l_gmshTet4   = l_gmsh.getGmshType( edge_v::TET4 );

  edge_v::t_entityType l_edgePoint  = l_gmsh.getEntityType( l_gmshPoint );
  edge_v::t_entityType l_edgeLine   = l_gmsh.getEntityType( l_gmshLine );
  edge_v::t_entityType l_edgeQuad4r = l_gmsh.getEntityType( l_gmshQuad4r );
  edge_v::t_entityType l_edgeTria3  = l_gmsh.getEntityType( l_gmshTria3 );
  edge_v::t_entityType l_edgeHex8r  = l_gmsh.getEntityType( l_gmshHex8r );
  edge_v::t_entityType l_edgeTet4   = l_gmsh.getEntityType( l_gmshTet4 );

  REQUIRE( l_edgePoint  == edge_v::POINT );
  REQUIRE( l_edgeLine   == edge_v::LINE );
  REQUIRE( l_edgeQuad4r == edge_v::QUAD4R );
  REQUIRE( l_edgeTria3  == edge_v::TRIA3 );
  REQUIRE( l_edgeHex8r  == edge_v::HEX8R );
  REQUIRE( l_edgeTet4   == edge_v::TET4 );
}
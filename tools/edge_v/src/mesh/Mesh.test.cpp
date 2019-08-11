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
 * Tests the mesh functions.
 **/
#include <catch.hpp>
#define private public
#include "Mesh.h"
#undef private

namespace edge_v {
  namespace test {
    extern std::string g_files;
  }
}

TEST_CASE( "Tests the mesh interface for triangular meshes.", "[mesh][Tria3]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the mesh-interface
    edge_v::mesh::Mesh l_mesh( l_path );

    REQUIRE( l_mesh.nVes() == 21 );
    REQUIRE( l_mesh.nFas() == 48 );
    REQUIRE( l_mesh.nEls() == 28 );
  }
}

TEST_CASE( "Tests the mesh interface in the case of periodic boundaries.", "[mesh][periodic]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the mesh-interface without periodic boundaries
    edge_v::mesh::Mesh l_mesh0( l_path );

    // check three faces
    REQUIRE( l_mesh0.getFaEl()[0*2 + 0] == 5 );
    REQUIRE( l_mesh0.getFaEl()[0*2 + 1] == std::numeric_limits< std::size_t >::max() );

    REQUIRE( l_mesh0.getFaEl()[1*2 + 0] == 18 );
    REQUIRE( l_mesh0.getFaEl()[1*2 + 1] == std::numeric_limits< std::size_t >::max() );

    REQUIRE( l_mesh0.getFaEl()[2*2 + 0] == 1 );
    REQUIRE( l_mesh0.getFaEl()[2*2 + 1] == std::numeric_limits< std::size_t >::max() );

    REQUIRE( l_mesh0.getFaEl()[3*2 + 0] == 3 );
    REQUIRE( l_mesh0.getFaEl()[3*2 + 1] == std::numeric_limits< std::size_t >::max() );

    // construct the mesh-interface with periodic boundaries
    edge_v::mesh::Mesh l_mesh1( l_path,
                                true );

    REQUIRE( l_mesh1.getFaEl()[0*2 + 0] == 2 );
    REQUIRE( l_mesh1.getFaEl()[0*2 + 1] == 5 );

    REQUIRE( l_mesh1.getFaEl()[1*2 + 0] == 17 );
    REQUIRE( l_mesh1.getFaEl()[1*2 + 1] == 18 );

    REQUIRE( l_mesh1.getFaEl()[2*2 + 0] == 1 );
    REQUIRE( l_mesh1.getFaEl()[2*2 + 1] == 7 );

    REQUIRE( l_mesh1.getFaEl()[3*2 + 0] == 3 );
    REQUIRE( l_mesh1.getFaEl()[3*2 + 1] == 4 );
  }
}
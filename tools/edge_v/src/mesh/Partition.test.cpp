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
 * Tests the mesh partitioning.
 **/
#include <catch.hpp>
#define private public
#include "Partition.h"
#undef private

namespace edge_v {
  namespace test {
    extern std::string g_files;
  }
}

TEST_CASE( "Tests the PartGraphKway call of the partitioner.", "[partition][kWay]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the mesh-interface
    edge_v::mesh::Mesh l_mesh( l_path );

    // construct partitioner
    edge_v::mesh::Partition l_part( l_mesh );
 
    // test partitioner
    l_part.kWay( 3 );
    std::size_t const * l_elPa = l_part.getElPa();

    for( std::size_t l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      REQUIRE( l_elPa[l_el] < 3 );
    }

    // assign time groups
    unsigned short * l_elTg;
    l_elTg = new unsigned short[ l_mesh.nEls() ];

    for( std::size_t l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      l_elTg[l_el] = l_el%4;
    }

    // test with assigned time groups
    l_part.kWay( 4, l_elTg );
    l_elPa = l_part.getElPa();
    for( std::size_t l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      REQUIRE( l_elPa[l_el] < 4 );
    }

    delete[] l_elTg;
  }
}
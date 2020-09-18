/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
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
    std::string l_path = edge_v::test::g_files + "/la_habra_small.msh";

    // construct the mesh-interface
    edge_v::io::Gmsh l_gmsh;
    l_gmsh.open( l_path );
    l_gmsh.readMesh();

    edge_v::mesh::Mesh l_mesh( l_gmsh );

    // construct partitioner
    edge_v::mesh::Partition l_part0( l_mesh, nullptr );
 
    // test partitioner
    l_part0.kWay( 3 );
    edge_v::t_idx const * l_elPa = l_part0.getElPa();

    for( edge_v::t_idx l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      REQUIRE( l_elPa[l_el] < 3 );
    }

    // assign time groups
    unsigned short * l_elTg;
    l_elTg = new unsigned short[ l_mesh.nEls() ];

    for( edge_v::t_idx l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      l_elTg[l_el] = l_el%4;
    }

    // construct partitioner
    edge_v::mesh::Partition l_part1( l_mesh, l_elTg );

    // test with assigned time groups
    l_part1.kWay( 4 );
    l_elPa = l_part1.getElPa();
    for( edge_v::t_idx l_el = 0; l_el < l_mesh.nEls(); l_el++ ) {
      REQUIRE( l_elPa[l_el] < 4 );
    }

    delete[] l_elTg;
  }
}

TEST_CASE( "Tests the computation of element priorities.", "[partition][elPr]" ) {
  /*
   * Our test case.
   * The adjacency info is not consistent.
   *
   * adjusted ordering:
   *
   * id | rank | tg  | send? | prio | ini id
   *  0 |   0  |  0  |  0    | 0    | 6
   *  1 |   0  |  1  |  0    | 1    | 3
   *  2 |   0  |  0  |  1    | 3    | 2
   *  3 |   1  |  1  |  0    | 7    | 9
   *  4 |   1  |  2  |  0    | 8    | 4
   *  5 |   1  |  0  |  1    | 9    | 5
   *  6 |   1  |  1  |  1    | 10   | 1
   *  7 |   1  |  1  |  1    | 10   | 8
   *  8 |   2  |  1  |  1    | 16   | 7
   *  9 |   3  |  2  |  1    | 23   | 0
   */
  //                            0  1  2  3  4  5  6  7  8  9
  edge_v::t_idx  l_elPa[10] = { 3, 1, 0, 0, 1, 1, 0, 2, 1, 1 };
  unsigned short l_elTg[10] = { 2, 1, 0, 1, 2, 0, 0, 1, 1, 1 };
  edge_v::t_idx l_elFaEl[10][3] = {
    // 0
    { 7,
      std::numeric_limits< edge_v::t_idx >::max(),
      2 },
    // 1
    { 0,
      4,
      2 },
    // 2
    { 7,
      6,
      8 },
    // 3
    { std::numeric_limits< edge_v::t_idx >::max(),
      6,
      2 },
    // 4
    { 9,
      5,
      1 },
    // 5
    { 4,
      0,
      1 },
    // 6
    { std::numeric_limits< edge_v::t_idx >::max(),
      3,
      std::numeric_limits< edge_v::t_idx >::max() },
    // 7
    { 0,
      1,
      2 },
    // 8
    { 1,
      3,
      5 },
    // 9
    { 5,
      1,
      4 }
  };

  edge_v::t_idx l_elPr[10] = {0};

  edge_v::mesh::Partition::getElPr( edge_v::TRIA3,
                                    10,
                                    l_elFaEl[0],
                                    l_elPa,
                                    l_elTg,
                                    l_elPr );

  // check the priorities
  REQUIRE( l_elPr[6] ==  0 );
  REQUIRE( l_elPr[3] ==  1 );
  REQUIRE( l_elPr[2] ==  3 );
  REQUIRE( l_elPr[9] ==  7 );
  REQUIRE( l_elPr[4] ==  8 );
  REQUIRE( l_elPr[5] ==  9 );
  REQUIRE( l_elPr[1] == 10 );
  REQUIRE( l_elPr[8] == 10 );
  REQUIRE( l_elPr[7] == 16 );
  REQUIRE( l_elPr[0] == 23 );
}
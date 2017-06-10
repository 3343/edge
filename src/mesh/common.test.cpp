/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Unit tests of common mesh functions.
 **/

#include <catch.hpp>
#include "constants.hpp"

#define private public
#include "common.hpp"
#undef private

TEST_CASE( "Find minimum global id of a point", "[findMinGid]" ) {
  /*
   * Our dummy quad4r-example:
   *
   * time group 0: 0-2
   * time group 1: 3-4
   *
   * 0, 3, 4 are inner elements
   * 1 is send-eleemnt
   * 2 is receive-element
   *
   * global ids are given by:
   *
   * 0: 3
   * 1: 2
   * 2: 4
   * 3: 1
   * 4: 0
   *
   *     0.0  0.5 1.0 1.5
   *      |   |   |   |
   *  2.0-0***3***6   |  
   *      * 0 * 2 *   |
   *  1.0-1***4***7**10
   *      * 1 * 3 * 4 *
   *  0.0-2***5***8***9
   */

  // setup dummy element layout
  t_enLayout l_elLayout;
  l_elLayout.nEnts = 5;
  l_elLayout.timeGroups.resize(2);

  l_elLayout.timeGroups[0].nEntsOwn    = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[0].inner.first = 0;
  l_elLayout.timeGroups[0].inner.size  = 1;

  l_elLayout.timeGroups[0].send.resize(1);
  l_elLayout.timeGroups[0].send[0].first = 1;
  l_elLayout.timeGroups[0].send[0].size  = 1;

  l_elLayout.timeGroups[0].receive.resize(1);
  l_elLayout.timeGroups[0].receive[0].first = 2;
  l_elLayout.timeGroups[0].receive[0].size  = 1;

  l_elLayout.timeGroups[1].nEntsOwn    = 2;
  l_elLayout.timeGroups[1].nEntsNotOwn = 0;
  l_elLayout.timeGroups[1].inner.first = 3;
  l_elLayout.timeGroups[1].inner.size  = 2;

  // assign vertex chars
  t_vertexChars l_veChars[11];
  l_veChars[ 0].coords[0] = 0.0; l_veChars[ 0].coords[1] = 2.0; l_veChars[ 0].coords[2] = 0.0;
  l_veChars[ 1].coords[0] = 0.0; l_veChars[ 1].coords[1] = 1.0; l_veChars[ 1].coords[2] = 0.0;
  l_veChars[ 2].coords[0] = 0.0; l_veChars[ 2].coords[1] = 0.0; l_veChars[ 2].coords[2] = 0.0;
  l_veChars[ 3].coords[0] = 0.5; l_veChars[ 3].coords[1] = 2.0; l_veChars[ 3].coords[2] = 0.0;
  l_veChars[ 4].coords[0] = 0.5; l_veChars[ 4].coords[1] = 1.0; l_veChars[ 4].coords[2] = 0.0;
  l_veChars[ 5].coords[0] = 0.5; l_veChars[ 5].coords[1] = 0.0; l_veChars[ 5].coords[2] = 0.0;
  l_veChars[ 6].coords[0] = 1.0; l_veChars[ 6].coords[1] = 2.0; l_veChars[ 6].coords[2] = 0.0;
  l_veChars[ 7].coords[0] = 1.0; l_veChars[ 7].coords[1] = 1.0; l_veChars[ 7].coords[2] = 0.0;
  l_veChars[ 8].coords[0] = 1.0; l_veChars[ 8].coords[1] = 0.0; l_veChars[ 8].coords[2] = 0.0;
  l_veChars[ 9].coords[0] = 1.5; l_veChars[ 9].coords[1] = 0.0; l_veChars[ 9].coords[2] = 0.0;
  l_veChars[10].coords[0] = 1.5; l_veChars[10].coords[1] = 1.0; l_veChars[10].coords[2] = 0.0;

  // assign el-ve adjacency
  int_el l_elVe[5][4];
  l_elVe[0][0] = 0; l_elVe[0][1] = 1; l_elVe[0][2] =  4; l_elVe[0][3] = 3;
  l_elVe[1][0] = 5; l_elVe[1][1] = 4; l_elVe[1][2] =  1; l_elVe[1][3] = 2;
  l_elVe[2][0] = 7; l_elVe[2][1] = 6; l_elVe[2][2] =  3; l_elVe[2][3] = 4;
  l_elVe[3][0] = 4; l_elVe[3][1] = 5; l_elVe[3][2] =  8; l_elVe[3][3] = 7;
  l_elVe[4][0] = 8; l_elVe[4][1] = 9; l_elVe[4][2] = 10; l_elVe[4][3] = 7;

  // assign global ids
  int_gid l_gIdsEl[5];
  l_gIdsEl[0] = 3; l_gIdsEl[1] = 2; l_gIdsEl[2] = 4; l_gIdsEl[3] = 1; l_gIdsEl[4] = 0;

  // do the tests
  real_mesh l_crds[3]; l_crds[2] = 0;
  bool l_ownedT; int_el l_lIdT; int_gid l_gIdT; int_tg l_tgT;

  /*
   * point outside mesh
   */
  l_crds[0] = -1.0; l_crds[1] = 0.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == std::numeric_limits< int_el  >::max() );
  REQUIRE( l_gIdT == std::numeric_limits< int_gid >::max() );
  REQUIRE( l_tgT  == std::numeric_limits< int_tg  >::max() );

  l_crds[0] = 0.25; l_crds[1] = 2.01;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == std::numeric_limits< int_el  >::max() );
  REQUIRE( l_gIdT == std::numeric_limits< int_gid >::max() );
  REQUIRE( l_tgT  == std::numeric_limits< int_tg  >::max() );

  l_crds[0] = 1.25; l_crds[1] = 1.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == std::numeric_limits< int_el  >::max() );
  REQUIRE( l_gIdT == std::numeric_limits< int_gid >::max() );
  REQUIRE( l_tgT  == std::numeric_limits< int_tg  >::max() );


  /*
   * point inside quad
   */
  l_crds[0] = 0.25; l_crds[1] = 0.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 1 );
  REQUIRE( l_gIdT == 2 );
  REQUIRE( l_tgT  == 0 );

  l_crds[0] = 0.25; l_crds[1] = 1.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 0 );
  REQUIRE( l_gIdT == 3 );
  REQUIRE( l_tgT  == 0 );

  l_crds[0] = 0.75; l_crds[1] = 0.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 3 );
  REQUIRE( l_gIdT == 1 );
  REQUIRE( l_tgT  == 1 );

  l_crds[0] = 0.75; l_crds[1] = 1.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == 2 );
  REQUIRE( l_gIdT == 4 );
  REQUIRE( l_tgT  == 0 );

  l_crds[0] = 1.25; l_crds[1] = 0.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 4 );
  REQUIRE( l_gIdT == 0 );
  REQUIRE( l_tgT  == 1 );

  /*
   * point on face
   */
  l_crds[0] = 0.5; l_crds[1] = 0.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 3 );
  REQUIRE( l_gIdT == 1 );
  REQUIRE( l_tgT  == 1 );

  l_crds[0] = 0.6; l_crds[1] = 1.0;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                     l_elLayout,
                                                     l_elVe[0],
                                                     l_veChars,
                                                     l_gIdsEl,
                                                     l_lIdT,
                                                     l_gIdT,
                                                     l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 3 );
  REQUIRE( l_gIdT == 1 );
  REQUIRE( l_tgT  == 1 );

  l_crds[0] = 1.0; l_crds[1] = 1.5;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == 2 );
  REQUIRE( l_gIdT == 4 );
  REQUIRE( l_tgT  == 0 );

  /*
   * point on vertex
   */
  l_crds[0] = 1.0; l_crds[1] = 1.0;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 4 );
  REQUIRE( l_gIdT == 0 );
  REQUIRE( l_tgT  == 1 );

  l_crds[0] = 0.5; l_crds[1] = 1.0;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 3 );
  REQUIRE( l_gIdT == 1 );
  REQUIRE( l_tgT  == 1 );

  l_crds[0] = 0.5; l_crds[1] = 2.0;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == true );
  REQUIRE( l_lIdT == 0 );
  REQUIRE( l_gIdT == 3 );
  REQUIRE( l_tgT  == 0 );

  l_crds[0] = 1.; l_crds[1] = 2.0;
  l_ownedT = edge::mesh::common< QUAD4R >::findMinGid( l_crds,
                                                       l_elLayout,
                                                       l_elVe[0],
                                                       l_veChars,
                                                       l_gIdsEl,
                                                       l_lIdT,
                                                       l_gIdT,
                                                       l_tgT );
  REQUIRE( l_ownedT == false );
  REQUIRE( l_lIdT == 2 );
  REQUIRE( l_gIdT == 4 );
  REQUIRE( l_tgT  == 0 );
}

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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Tests base discretization using quads or hexes.
 **/

#include <catch.hpp>

#define private public
#include "Base.h"
#undef private

TEST_CASE( "Quadrilateral base discretization", "[regBase][quad]" ) {
  unsigned int l_nX[2], l_nPartX[2], l_partX[2];

  // check the partitioning
  unsigned int l_nXAll[2];

  l_nXAll[0] = 50; l_nXAll[1] = 50;
  edge::mesh::regular::Base::getSetup2d( 4, 0, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  2 );
  REQUIRE( l_nPartX[1] ==  2 );
  REQUIRE( l_partX[0]  ==  0 );
  REQUIRE( l_partX[1]  ==  0 );
  REQUIRE( l_nX[0]     == 25 );
  REQUIRE( l_nX[1]     == 25 );

  l_nXAll[0] = 50; l_nXAll[1] = 25;
  edge::mesh::regular::Base::getSetup2d( 2, 1, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  2 );
  REQUIRE( l_nPartX[1] ==  1 );
  REQUIRE( l_partX[0]  ==  1 );
  REQUIRE( l_partX[1]  ==  0 );
  REQUIRE( l_nX[0]     == 25 );
  REQUIRE( l_nX[1]     == 25 );

  l_nXAll[0] = 25; l_nXAll[1] = 50;
  edge::mesh::regular::Base::getSetup2d( 2, 1, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  1 );
  REQUIRE( l_nPartX[1] ==  2 );
  REQUIRE( l_partX[0]  ==  0 );
  REQUIRE( l_partX[1]  ==  1 );
  REQUIRE( l_nX[0]     == 25 );
  REQUIRE( l_nX[1]     == 25 );

  l_nXAll[0] = 250; l_nXAll[1] = 500;
  edge::mesh::regular::Base::getSetup2d( 24, 13, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  4 );
  REQUIRE( l_nPartX[1] ==  6 );
  REQUIRE( l_partX[0]  ==  1 );
  REQUIRE( l_partX[1]  ==  3 );
  REQUIRE( l_nX[0]     == 63 );
  REQUIRE( l_nX[1]     == 83 );

  l_nXAll[0] = 5000; l_nXAll[1] = 2500;
  edge::mesh::regular::Base::getSetup2d( 235, 234, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  47 );
  REQUIRE( l_nPartX[1] ==   5 );
  REQUIRE( l_partX[0]  ==  46 );
  REQUIRE( l_partX[1]  ==   4 );
  REQUIRE( l_nX[0]     == 106 );
  REQUIRE( l_nX[1]     == 500 );

  l_nXAll[0] = 5000; l_nXAll[1] = 2500;
  edge::mesh::regular::Base::getSetup2d( 2048, 31, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] == 64 );
  REQUIRE( l_nPartX[1] == 32 );
  REQUIRE( l_partX[0]  == 31 );
  REQUIRE( l_partX[1]  ==  0 );
  REQUIRE( l_nX[0]     == 78 );
  REQUIRE( l_nX[1]     == 79 );
}

TEST_CASE( "Hexahedral base discretization", "[regBase][hex]" ) {
  unsigned int l_nX[3], l_nPartX[3], l_partX[3];
  unsigned int l_nXAll[3];

  l_nXAll[0] = 50; l_nXAll[1] = 50; l_nXAll[2] = 50;
  edge::mesh::regular::Base::getSetup3d( 4, 0, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  1 );
  REQUIRE( l_nPartX[1] ==  1 );
  REQUIRE( l_nPartX[2] ==  4 );
  REQUIRE( l_partX[0]  ==  0 );
  REQUIRE( l_partX[1]  ==  0 );
  REQUIRE( l_partX[2]  ==  0 );
  REQUIRE( l_nX[0]     == 50 );
  REQUIRE( l_nX[1]     == 50 );
  REQUIRE( l_nX[2]     == 13 );

  l_nXAll[0] = 50; l_nXAll[1] = 50; l_nXAll[2] = 50;
  edge::mesh::regular::Base::getSetup3d( 4, 2, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==  1 );
  REQUIRE( l_nPartX[1] ==  1 );
  REQUIRE( l_nPartX[2] ==  4 );
  REQUIRE( l_partX[0]  ==  0 );
  REQUIRE( l_partX[1]  ==  0 );
  REQUIRE( l_partX[2]  ==  2 );
  REQUIRE( l_nX[0]     == 50 );
  REQUIRE( l_nX[1]     == 50 );
  REQUIRE( l_nX[2]     == 12 );

  l_nXAll[0] = 300; l_nXAll[1] = 500; l_nXAll[2] = 250;
  edge::mesh::regular::Base::getSetup3d( 16, 0, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==   2 );
  REQUIRE( l_nPartX[1] ==   4 );
  REQUIRE( l_nPartX[2] ==   2 );
  REQUIRE( l_partX[0]  ==   0 );
  REQUIRE( l_partX[1]  ==   0 );
  REQUIRE( l_partX[2]  ==   0 );
  REQUIRE( l_nX[0]     == 150 );
  REQUIRE( l_nX[1]     == 125 );
  REQUIRE( l_nX[2]     == 125 );

  l_nXAll[0] = 500; l_nXAll[1] = 250; l_nXAll[2] = 300;
  edge::mesh::regular::Base::getSetup3d( 24, 13, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==   6 );
  REQUIRE( l_nPartX[1] ==   2 );
  REQUIRE( l_nPartX[2] ==   2 );
  REQUIRE( l_partX[0]  ==   1 );
  REQUIRE( l_partX[1]  ==   0 );
  REQUIRE( l_partX[2]  ==   1 );
  REQUIRE( l_nX[0]     ==  84 );
  REQUIRE( l_nX[1]     == 125 );
  REQUIRE( l_nX[2]     == 150 );

  l_nXAll[0] = 5000; l_nXAll[1] = 2500; l_nXAll[2] = 3121;
  edge::mesh::regular::Base::getSetup3d( 27, 11, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==    3 );
  REQUIRE( l_nPartX[1] ==    3 );
  REQUIRE( l_nPartX[2] ==    3 );
  REQUIRE( l_partX[0]  ==    2 );
  REQUIRE( l_partX[1]  ==    0 );
  REQUIRE( l_partX[2]  ==    1 );
  REQUIRE( l_nX[0]     == 1666 );
  REQUIRE( l_nX[1]     ==  834 );
  REQUIRE( l_nX[2]     == 1040 );

  l_nXAll[0] = 5000; l_nXAll[1] = 6000; l_nXAll[2] = 3121;
  edge::mesh::regular::Base::getSetup3d( 108, 61, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==    3 );
  REQUIRE( l_nPartX[1] ==   12 );
  REQUIRE( l_nPartX[2] ==    3 );
  REQUIRE( l_partX[0]  ==    1 );
  REQUIRE( l_partX[1]  ==    8 );
  REQUIRE( l_partX[2]  ==    1 );
  REQUIRE( l_nX[0]     == 1667 );
  REQUIRE( l_nX[1]     ==  500 );
  REQUIRE( l_nX[2]     == 1040 );

  l_nXAll[0] = 5000; l_nXAll[1] = 2500; l_nXAll[2] = 7100;
  edge::mesh::regular::Base::getSetup3d( 2048, 31, true, l_nXAll, l_nX, l_nPartX, l_partX );
  REQUIRE( l_nPartX[0] ==   8 );
  REQUIRE( l_nPartX[1] ==   8 );
  REQUIRE( l_nPartX[2] ==  32 );
  REQUIRE( l_partX[0]  ==   7 );
  REQUIRE( l_partX[1]  ==   3 );
  REQUIRE( l_partX[2]  ==   0 );
  REQUIRE( l_nX[0]     == 625 );
  REQUIRE( l_nX[1]     == 313 );
  REQUIRE( l_nX[2]     == 222 );
}

TEST_CASE( "Vertex ids of the hexahedral base discretization", "[regBase][veIdsHex]" ) {
  edge::mesh::regular::Base l_base;

  int_el l_ves[8];
  l_base.getVesHex( -1, -1, -1,
                    17, 13,  7,
                    l_ves );

  REQUIRE( l_ves[0] == 0  );
  REQUIRE( l_ves[1] == 1  );
  REQUIRE( l_ves[3] == 20 );
  REQUIRE( l_ves[2] == 21 );
  REQUIRE( l_ves[4] == 320 );
  REQUIRE( l_ves[5] == 321 );
  REQUIRE( l_ves[7] == 340 );
  REQUIRE( l_ves[6] == 341 );
}

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
 * Unit tests for the entity layout.
 **/

#include <catch.hpp>
#include "EntityLayout.h"

TEST_CASE( "Entity layout: Completion of a partial layout.", "[sizesToLayout][EntityLayout]" ) {
  t_enLayout l_layout;

  // assign the partial layout
  l_layout.timeGroups.resize( 2 );

  l_layout.timeGroups[0].inner.size  =  5;

  l_layout.timeGroups[0].neRanks.resize( 3 );
  l_layout.timeGroups[0].neRanks[0] = 13;
  l_layout.timeGroups[0].neRanks[1] = 13;
  l_layout.timeGroups[0].neRanks[2] = 2;

  l_layout.timeGroups[0].neTgs.resize( 3 );
  l_layout.timeGroups[0].neTgs[0] = 0;
  l_layout.timeGroups[0].neTgs[1] = 1;
  l_layout.timeGroups[0].neTgs[2] = 2;

  l_layout.timeGroups[0].send.resize(    3 );
  l_layout.timeGroups[0].receive.resize( 3 );

  l_layout.timeGroups[0].send[0].size     = 3;
  l_layout.timeGroups[0].send[1].size     = 2;
  l_layout.timeGroups[0].send[2].size     = 3;
  l_layout.timeGroups[0].receive[0].size  = 3;
  l_layout.timeGroups[0].receive[1].size  = 2;
  l_layout.timeGroups[0].receive[2].size  = 2;

  l_layout.timeGroups[1].inner.size  = 9;

  l_layout.timeGroups[1].neRanks.resize( 3 );
  l_layout.timeGroups[1].neRanks[0] = 4;
  l_layout.timeGroups[1].neRanks[1] = 5;
  l_layout.timeGroups[1].neRanks[2] = 5;

  l_layout.timeGroups[1].neTgs.resize( 3 );
  l_layout.timeGroups[1].neTgs[0] = 3;
  l_layout.timeGroups[1].neTgs[1] = 1;
  l_layout.timeGroups[1].neTgs[2] = 2;

  l_layout.timeGroups[1].send.resize(    3 );
  l_layout.timeGroups[1].receive.resize( 3 );

  l_layout.timeGroups[1].send[0].size     = 2;
  l_layout.timeGroups[1].send[1].size     = 1;
  l_layout.timeGroups[1].send[2].size     = 1;
  l_layout.timeGroups[1].receive[0].size  = 2;
  l_layout.timeGroups[1].receive[1].size  = 2;
  l_layout.timeGroups[1].receive[2].size  = 1;

  // complete the layout
  edge::data::EntityLayout::sizesToLayout( l_layout );

  // check the resulting layout
  REQUIRE( l_layout.nEnts == 38 );

  REQUIRE( l_layout.timeGroups[0].nEntsOwn    == 13 );
  REQUIRE( l_layout.timeGroups[0].nEntsNotOwn ==  7 );
  REQUIRE( l_layout.timeGroups[0].inner.first ==  0 );
  REQUIRE( l_layout.timeGroups[0].inner.size  ==  5 );

  REQUIRE( l_layout.timeGroups[0].neRanks[0] == 13 );
  REQUIRE( l_layout.timeGroups[0].neRanks[1] == 13 );
  REQUIRE( l_layout.timeGroups[0].neRanks[2] == 2 );

  REQUIRE( l_layout.timeGroups[0].neTgs[0] == 0 );
  REQUIRE( l_layout.timeGroups[0].neTgs[1] == 1 );
  REQUIRE( l_layout.timeGroups[0].neTgs[2] == 2 );

  REQUIRE( l_layout.timeGroups[0].send[0].first    == 5 );
  REQUIRE( l_layout.timeGroups[0].send[0].size     == 3 );

  REQUIRE( l_layout.timeGroups[0].send[1].first    == 8 );
  REQUIRE( l_layout.timeGroups[0].send[1].size     == 2 );

  REQUIRE( l_layout.timeGroups[0].send[2].first    == 10 );
  REQUIRE( l_layout.timeGroups[0].send[2].size     == 3  );

  REQUIRE( l_layout.timeGroups[0].receive[0].first == 13 );
  REQUIRE( l_layout.timeGroups[0].receive[0].size  == 3  );

  REQUIRE( l_layout.timeGroups[0].receive[1].first == 16 );
  REQUIRE( l_layout.timeGroups[0].receive[1].size  == 2  );

  REQUIRE( l_layout.timeGroups[0].receive[2].first == 18 );
  REQUIRE( l_layout.timeGroups[0].receive[2].size  == 2  );

  REQUIRE( l_layout.timeGroups[1].nEntsOwn    == 13 );
  REQUIRE( l_layout.timeGroups[1].nEntsNotOwn == 5  );
  REQUIRE( l_layout.timeGroups[1].inner.first == 20 );
  REQUIRE( l_layout.timeGroups[1].inner.size  == 9  );

  REQUIRE( l_layout.timeGroups[1].neRanks[0] == 4 );
  REQUIRE( l_layout.timeGroups[1].neRanks[1] == 5 );
  REQUIRE( l_layout.timeGroups[1].neRanks[2] == 5 );

  REQUIRE( l_layout.timeGroups[1].neTgs[0] == 3 );
  REQUIRE( l_layout.timeGroups[1].neTgs[1] == 1 );
  REQUIRE( l_layout.timeGroups[1].neTgs[2] == 2 );

  REQUIRE( l_layout.timeGroups[1].send[0].first    == 29 );
  REQUIRE( l_layout.timeGroups[1].send[0].size     == 2  );

  REQUIRE( l_layout.timeGroups[1].send[1].first    == 31 );
  REQUIRE( l_layout.timeGroups[1].send[1].size     == 1  );

  REQUIRE( l_layout.timeGroups[1].send[2].first    == 32 );
  REQUIRE( l_layout.timeGroups[1].send[2].size     == 1  );

  REQUIRE( l_layout.timeGroups[1].receive[0].first == 33 );
  REQUIRE( l_layout.timeGroups[1].receive[0].size  == 2  );

  REQUIRE( l_layout.timeGroups[1].receive[1].first == 35 );
  REQUIRE( l_layout.timeGroups[1].receive[1].size  == 2  );

  REQUIRE( l_layout.timeGroups[1].receive[2].first == 37 );
  REQUIRE( l_layout.timeGroups[1].receive[2].size  == 1  );
}

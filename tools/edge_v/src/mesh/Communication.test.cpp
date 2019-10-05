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
 * Tests the derivation of communication structures.
 **/
#include <catch.hpp>
#define private public
#include "Communication.h"
#undef private

TEST_CASE( "Tests the derivation of the partition-time group pairs.", "[communication][getPaTgPairs]" ) {
  /*
   * Our test case:
   *
   *  id | elFaEl   | partition | time group | pairs
   *   0 |          | 0         | 2          |
   *   1 |          | 0         | 2          |
   *   2 |          | 1         | 0          |
   *  ==============|===========|============|
   *   3 |  6  4  2 | 2         | 1          | (1, 0)
   *   4 |  1  5  9 | 2         | 1          | (0, 2), (4, 2)
   *   5 |  7  2 12 | 2         | 1          | (1, 0), (5, 0)
   *   6 | 10  5  4 | 2         | 1          | (4, 1)
   *   7 |  3 11  0 | 2         | 1          | (5, 2), (0, 2)
   *  ==============|===========|============|
   *   8 |          | 3         | 1          |
   *   9 |          | 4         | 2          |
   *  10 |          | 4         | 1          |
   *  11 |          | 5         | 2          |
   *  12 |          | 5         | 0
   *
   * Sorted, unique pairs: (0, 2), (1, 0), (4, 1), (4, 2), (5, 0), (5, 2)
   */

  // assemble test structures
  std::size_t l_elFaEl[13][3];
  for( unsigned short l_el = 0; l_el < 13; l_el++ )
    for( unsigned short l_fa = 0; l_fa < 3; l_fa++ )
      l_elFaEl[l_el][l_fa] = std::numeric_limits< std::size_t >::max();

  l_elFaEl[3][0] =  6; l_elFaEl[3][1] =  4; l_elFaEl[3][2] =  2;
  l_elFaEl[4][0] =  1; l_elFaEl[4][1] =  5; l_elFaEl[4][2] =  9;
  l_elFaEl[5][0] =  7; l_elFaEl[5][1] =  2; l_elFaEl[5][2] = 12;
  l_elFaEl[6][0] = 10; l_elFaEl[6][1] =  5; l_elFaEl[6][2] =  4;
  l_elFaEl[7][0] =  3; l_elFaEl[7][1] = 11; l_elFaEl[7][2] =  0;

  //                            0  1  2  3  4  5  6  7  8  9 10 11 12
  std::size_t l_elPa[13]    = { 0, 0, 1, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5 };
  unsigned short l_elTg[13] = { 2, 2, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 0 };

  // get the unique and sort pairs
  std::set< std::pair< std::size_t, unsigned short > > l_pairs;
  edge_v::mesh::Communication::getPaTgPairs( edge_v::TRIA3,
                                             3,
                                             5,
                                             l_elFaEl[0],
                                             l_elTg,
                                             l_elPa,
                                             l_pairs );


  // check the results
  auto l_it = l_pairs.begin();

  REQUIRE( std::get<0>(*l_it) == 0 );
  REQUIRE( std::get<1>(*l_it) == 2 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 1 );
  REQUIRE( std::get<1>(*l_it) == 0 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 4 );
  REQUIRE( std::get<1>(*l_it) == 1 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 4 );
  REQUIRE( std::get<1>(*l_it) == 2 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 5 );
  REQUIRE( std::get<1>(*l_it) == 0 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 5 );
  REQUIRE( std::get<1>(*l_it) == 2 );
}
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
 * Unit tests for time step groups.
 **/
#include <catch.hpp>

#define private public
#include "Groups.h"
#undef private

TEST_CASE( "Tests the derivation of time step groups.", "[time][groups]" ) {
  /**
   * First example:
   * 
   *  
   *   unordered:
   *     0     1     2     3     4     5     6     7     8    9     10
   *    [12.3, 3.7,  10.4, 3.5,  5.0,  1.0 , 3.9,  3.6,  3.2, 1.5,  2.1]
   * 
   *   ordered:
   *     5     9     10    3     8     7     1     6     4    2     0
   *    [1.0,  1.5,  2.1,  3.5,  3.2,  3.6,  3.7,  3.9,  5.0, 10.4, 12.3] 
   *
   *    time groups:
   *      0 [1.0,  1.6 ] 1.6
   *      1 [1.6,  3.04] 1.9
   *      2 [3.04, 6.08] 2.0
   *      3 [6.08, 9.12] 1.5
   *      4 [9.12, inf]
   **/

  double l_rates1[4] = {1.6, 1.9, 2.0, 1.5};
  std::size_t l_nEls1 = 11;
  double l_ts1[11] = {12.3, 3.7,  10.4, 3.5,  5.0,  1.0 , 3.9,  3.6,  3.2, 1.5,  2.1};

  std::size_t l_elFaEl1[11*3];
  for( unsigned short l_en = 0; l_en < 11*3; l_en++ ) l_elFaEl1[l_en] = std::numeric_limits< std::size_t >::max();

  edge_v::time::Groups l_groups1( edge_v::TRIA3,
                                  11,
                                  l_elFaEl1,
                                  4,
                                  l_rates1,
                                  l_ts1 );

  REQUIRE( l_groups1.m_elTg[ 0] == 4 );
  REQUIRE( l_groups1.m_elTg[ 1] == 2 );
  REQUIRE( l_groups1.m_elTg[ 2] == 4 );
  REQUIRE( l_groups1.m_elTg[ 3] == 2 );
  REQUIRE( l_groups1.m_elTg[ 4] == 2 );
  REQUIRE( l_groups1.m_elTg[ 5] == 0 );
  REQUIRE( l_groups1.m_elTg[ 6] == 2 );
  REQUIRE( l_groups1.m_elTg[ 7] == 2 );
  REQUIRE( l_groups1.m_elTg[ 8] == 2 );
  REQUIRE( l_groups1.m_elTg[ 9] == 0 );
  REQUIRE( l_groups1.m_elTg[10] == 1 );

  double l_load[3] = {};
  edge_v::time::Groups::getLoads( 11,
                                  l_ts1,
                                  l_groups1.m_tsIntervals,
                                  l_groups1.m_elTg,
                                  l_load[0],
                                  l_load[1],
                                  l_load[2] );

  REQUIRE( l_load[0] == Approx(11) );
  REQUIRE( l_load[1] == Approx(2 / 1.0 + 1 / 1.6 + 6 / 3.04 + 0 / 6.08 + 2 / 9.12 ) );
  REQUIRE( l_load[2] == Approx(1 / 12.3 + 1 / 3.7 + 1 / 10.4 + 1 / 3.5 + 1 / 5.0 + 1 / 1.0 + 1 / 3.9 + 1 / 3.6 + 1 / 3.2 + 1 / 1.5 + 1 /  2.1) );

  std::size_t l_elFaEl2[11*3];
  for( unsigned short l_en = 0; l_en < 11*3; l_en++ ) l_elFaEl2[l_en] = std::numeric_limits< std::size_t >::max();
  l_elFaEl2[0*3 + 0] =  2;
  l_elFaEl2[4*3 + 0] =  5;
  l_elFaEl2[2*3 + 1] = 10;
  edge_v::time::Groups l_groups2( edge_v::TRIA3,
                                  11,
                                  l_elFaEl2,
                                  4,
                                  l_rates1,
                                  l_ts1 );
  REQUIRE( l_groups2.m_elTg[ 0] == 3 );
  REQUIRE( l_groups2.m_elTg[ 1] == 2 );
  REQUIRE( l_groups2.m_elTg[ 2] == 2 );
  REQUIRE( l_groups2.m_elTg[ 3] == 2 );
  REQUIRE( l_groups2.m_elTg[ 4] == 1 );
  REQUIRE( l_groups2.m_elTg[ 5] == 0 );
  REQUIRE( l_groups2.m_elTg[ 6] == 2 );
  REQUIRE( l_groups2.m_elTg[ 7] == 2 );
  REQUIRE( l_groups2.m_elTg[ 8] == 2 );
  REQUIRE( l_groups2.m_elTg[ 9] == 0 );
  REQUIRE( l_groups2.m_elTg[10] == 1 );
}
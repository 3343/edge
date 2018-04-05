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
 * Unit tests for the sub-grid operations.
 **/
#include <catch.hpp>
#define private public
#include "SubGrid.hpp"
#undef private

TEST_CASE( "Sub-grid: ptSc.", "[subGrid][ptSc]" ) {
  /*
   * Order: 4
   * #sub-cells: 7
   *
   *   0   1/7  2/7  3/7  4/7  5/7  6/7   1
   *   |----x----x----x----x----x----x----|
   *  sv0  sv1  sv2  sv3  sv4  sv5  sv6  sv7
   *
   *    sc5  sc0  sc1  sc2  sc3   sc4  sc6
   *
   *  sc |  sv
   *   0 | 1 2
   *   1 | 2 3
   *   2 | 3 4
   *   3 | 4 5
   *   4 | 5 6
   *   5 | 0 1
   *   6 | 6 7
   */
  unsigned short l_scSv1[7][2] = { {1,2},
                                   {2,3},
                                   {3,4},
                                   {4,5},
                                   {5,6},
                                   {0,1},
                                   {6,7} };

  typedef struct {
    double coords[1];
  } t_vcChars1;

  t_vcChars1 l_vcChars1[8] = { {{    0.0}},
                               {{1.0/7.0}},
                               {{2.0/7.0}},
                               {{3.0/7.0}},
                               {{4.0/7.0}},
                               {{5.0/7.0}},
                               {{6.0/7.0}},
                               {{    1.0}} };

  // point for which the sub-cell is queried
  double l_pt1[1];
  // resulting subcell
  unsigned short l_sc;

  // check results
  l_pt1[0] = -1.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 5 );

  l_pt1[0] = 0.25/7.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 5 );

  l_pt1[0] = 2.75/7.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 1 );

  l_pt1[0] = 4.1/7.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 3 );

  l_pt1[0] = 5.9/7.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 4 );

  l_pt1[0] = 6.3/7.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 6 );

  l_pt1[0] = 10.0;
  l_sc = edge::sc::SubGrid< LINE, 4 >::ptSc( l_pt1, l_scSv1, l_vcChars1 );
  REQUIRE( l_sc == 6 );
}

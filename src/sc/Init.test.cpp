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
 * Unit tests for initialization of sub-cell data.
 **/
#include <catch.hpp>
#include "Memory.hpp"
#define private public
#include "Init.hpp"
#undef private

TEST_CASE( "Sub-cell init: Scatter operation.", "[subCellInit][scatter]" ) {
  // check constants
  unsigned short l_const;
  l_const = edge::sc::Init< QUAD4R, 2 >::TL_N_SCS; REQUIRE( l_const == 9 );
  l_const = edge::sc::Init< QUAD4R, 2 >::TL_N_MDS; REQUIRE( l_const == 4 );

  // sub-cell data
  edge::sc::t_SubCell< int,
                       float,
                       T_SDISC.ELEMENT,
                       ORDER,
                       N_QUANTITIES,
                       N_CRUNS > l_sc1;

  // dynamic memory allocations
  edge::data::Dynamic l_dynMem;

  // init sc operators
  edge::sc::Init< T_SDISC.ELEMENT, ORDER >::ops( l_sc1.ops );

  // compilation specific (matrices are linkes) sanity checks
#if defined PP_T_ELEMENTS_QUAD4R
#if PP_ORDER==2
  // first row
  REQUIRE( l_sc1.ops.scatter[0][0] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][1] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][2] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][3] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][4] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][5] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][6] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][7] == Approx(      1.0 ) );
  REQUIRE( l_sc1.ops.scatter[0][8] == Approx(      1.0 ) );

  // second row
  REQUIRE( l_sc1.ops.scatter[1][0] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][1] == Approx( -2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][2] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][3] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][4] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][5] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][6] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][7] == Approx( -2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[1][8] == Approx( -2.0/3.0 ) );

  // third row
  REQUIRE( l_sc1.ops.scatter[2][0] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][1] == Approx( -2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][2] == Approx( -2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][3] == Approx( -2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][4] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][5] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][6] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][7] == Approx(  2.0/3.0 ) );
  REQUIRE( l_sc1.ops.scatter[2][8] == Approx(      0.0 ) );

  // third row
  REQUIRE( l_sc1.ops.scatter[3][0] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][1] == Approx(  4.0/9.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][2] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][3] == Approx( -4.0/9.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][4] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][5] == Approx(  4.0/9.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][6] == Approx(      0.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][7] == Approx( -4.0/9.0 ) );
  REQUIRE( l_sc1.ops.scatter[3][8] == Approx(      0.0 ) );

  // first column
  REQUIRE( l_sc1.ops.gather[0][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[1][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[2][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[3][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[4][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[5][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[6][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[7][0] == Approx(  1.0/9.0 ) );
  REQUIRE( l_sc1.ops.gather[8][0] == Approx(  1.0/9.0 ) );

#endif
#endif

  // check inverse in a generic way
  float l_inv[ N_ELEMENT_MODES ][ N_ELEMENT_MODES ];

  // init matrix
  for( unsigned short l_m1 = 0; l_m1 < N_ELEMENT_MODES; l_m1++ )
    for( unsigned short l_m2 = 0; l_m2 < N_ELEMENT_MODES; l_m2++ )
      l_inv[l_m1][l_m2] = 0;

  // computer matrix-matrix product
  for( unsigned short l_ro = 0; l_ro < N_ELEMENT_MODES; l_ro++ )
    for( unsigned short l_co = 0; l_co < N_ELEMENT_MODES; l_co++ )
      for( unsigned int l_sc = 0; l_sc < edge::sc::Init< T_SDISC.ELEMENT, ORDER >::TL_N_SCS; l_sc++ )
        l_inv[l_ro][l_co] += l_sc1.ops.scatter[l_ro][l_sc] * l_sc1.ops.gather[l_sc][l_co];

   // check for indentity
  for( unsigned short l_ro = 0; l_ro < N_ELEMENT_MODES; l_ro++ )
    for( unsigned short l_co = 0; l_co < N_ELEMENT_MODES; l_co++ ) {
      if( l_ro == l_co ) REQUIRE( l_inv[l_ro][l_co] == Approx(1.0) );
      else               REQUIRE( l_inv[l_ro][l_co] == Approx(0.0) );
    }
}

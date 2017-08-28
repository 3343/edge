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
 * Unit test of sub-cell resolution kernels.
 **/
#include <catch.hpp>
#define private public
#include "Kernels.hpp"
#undef private

TEST_CASE( "Sub-cell kernels: Scatter operation.", "[subCellKernels][scatter]" ) {
  /*
   * Dofs (2nd fused run opp sign):
   *        1.0 -0.1 -2.0 3.0
   *       -0.3  0.4  2.2 8.0
   *
   * Scatter:
   *    1.0  0.1  1.0  2.0  5.0  2.0  3.0  1.0  0.0
   *   -1.0 -0.1 -1.0 -2.0 -5.0 -2.0 -3.0 -1.0  0.0
   *    0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0
   *   -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2
   */

  // dg dofs
  float l_dgDofs1[2][4][2] = { { { 1.0, -1.0}, {-0.1,  0.1}, {-2.0,  2.0}, {3.0, -3.0} },
                               { {-0.3,  0.3}, { 0.4, -0.4}, { 2.2, -2.2}, {8.0, -8.0} } };

  // scatter operator
  float l_scatter1[4][9] = { { 1.0,  0.1,  1.0,   2.0,  5.0,  2.0,   3.0,  1.0,  0.0},
                             {-1.0, -0.1, -1.0,  -2.0, -5.0, -2.0,  -3.0, -1.0, -0.0},
                             { 0.2,  0.3,  0.4,   0.5,  0.6,  0.7,   0.8,  0.9,  1.0},
                             {-1.0, -0.9, -0.8,  -0.7, -0.6, -0.5,  -0.4, -0.3, -0.2} };

  // sub-cell DOFs
  float l_scDofs1[2][9][2];

  // perform scatter operation
  edge::sc::Kernels< QUAD4R, 2, 2, 2 >::scatterVanilla( l_dgDofs1, l_scatter1, l_scDofs1 );

  // result (1st fused run)
  float l_ut1[2][9] = { {  -2.3, -3.19,  -2.1,  -0.9,  2.5,   -0.7,    0.5,  -1.6, -2.6 },
                        { -8.26, -6.61, -6.22,  -5.9, -6.98, -3.86,  -3.54, -1.12,  0.6 } };

  // check the result
  for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
    for( unsigned short l_sc = 0; l_sc < 9; l_sc++ ) {
      REQUIRE( l_scDofs1[l_qt][l_sc][0] == Approx(  l_ut1[l_qt][l_sc] ) );
      REQUIRE( l_scDofs1[l_qt][l_sc][1] == Approx( -l_ut1[l_qt][l_sc] ) );
    }
  }
}

TEST_CASE( "Sub-cell kernels: Gather operation.", "[subCellKernels][gather]" ) {
  float l_scDofs1[3][9][2] = { { { -2.3,   2.3}, {-3.19,  3.19}, { -2.1,   2.1}, {-0.9,  0.9}, { 2.5,   -2.5}, { -0.7,   0.7}, {  0.5,  -0.5}, { -1.6,   1.6}, {-2.6,  2.6} },
                               { {-8.26,  8.26}, {-6.61,  6.61}, {-6.22,  6.22}, {-5.9,  5.9}, {-6.98,  6.98}, {-3.86,  3.86}, {-3.54,  3.54}, {-1.12,  1.12}, { 0.6, -0.6} },
                               { { 8.26, -8.26}, {-6.61,  6.61}, { 6.22, -6.22}, {-5.9,  5.9}, { 6.98, -6.98}, {-3.86,  3.86}, { 3.54, -3.54}, { 1.12, -1.12}, { 0.6, -0.6} } };

  float l_gather1[9][4] = { { 1.0,  2.0,  3.0,  4.0},
                            {-0.1, -0.2, -0.3, -0.4},
                            { 0.7, -0.2, -0.3,  0.8},

                            { 2.3, -4.2,  5.9,  9.0},
                            {-8.5, -9.0, -2.6,  1.0}, 
                            { 5.0,  2.1, -3.9,  5.1},

                            {-9.0,  5.2, -9.2,  8.1},
                            {-4.0, -5.0, -2.0,  0.0},
                            { 0.1,  0.2,  0.3, -0.4} };

  float l_dgDofs1[3][4][2];

  // perform gather operation
  edge::sc::Kernels< QUAD4R, 2, 3, 2 >::gatherVanilla( l_scDofs1, l_gather1, l_dgDofs1 );

  float l_ut1[3][4] = { { -28.631, -13.652, -16.573,  -13.684},
                        {  50.907,  52.852,  12.449, -144.052},
                        {-115.205, -16.62,  -47.635,    3.288} };

  // check result
  for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 4; l_md++ ) {
      REQUIRE( l_dgDofs1[l_qt][l_md][0] == Approx(  l_ut1[l_qt][l_md] ) );
      REQUIRE( l_dgDofs1[l_qt][l_md][1] == Approx( -l_ut1[l_qt][l_md] ) );
    }
  }
}

TEST_CASE( "Sub-cell kernels: Sub-cell extrema.", "[subCellKernels][scExtrema]" ) {
  float l_scDofs1[3][9][2] = { { { -2.3,   2.3}, {-3.19,  3.19}, { -2.1,   2.1}, {-0.9,  0.9}, { 2.5,   -2.5}, { -0.7,   0.7}, {  0.5,  -0.5}, { -1.6,   1.6}, {-2.6,  2.6} },
                               { {-8.26,  8.26}, {-6.61,  6.61}, {-6.22,  6.22}, {-5.9,  5.9}, {-6.98,  6.98}, {-3.86,  3.86}, {-3.54,  3.54}, {-1.12,  1.12}, { 0.6, -0.6} },
                               { { 8.26, -8.26}, {-6.61,  6.61}, { 6.22, -6.22}, {-5.9,  5.9}, { 6.98, -6.98}, {-3.86,  3.86}, { 3.54, -3.54}, { 1.12, -1.12}, { 0.6, -0.6} } };

  float l_min1[3][2], l_max1[3][2];

  // perform reduction
  edge::sc::Kernels< QUAD4R, 2, 3, 2 >::scExtrema( l_scDofs1, l_min1, l_max1 );

  // check results
  REQUIRE( l_min1[0][0] == Approx(-3.19) );
  REQUIRE( l_min1[1][0] == Approx(-8.26) );
  REQUIRE( l_min1[2][0] == Approx(-6.61) );

  REQUIRE( l_min1[0][1] == Approx( -2.5) );
  REQUIRE( l_min1[1][1] == Approx( -0.6) );
  REQUIRE( l_min1[2][1] == Approx(-8.26) );

  REQUIRE( l_max1[0][0] == Approx(  2.5) );
  REQUIRE( l_max1[1][0] == Approx(  0.6) );
  REQUIRE( l_max1[2][0] == Approx( 8.26) );

  REQUIRE( l_max1[0][1] == Approx( 3.19) );
  REQUIRE( l_max1[1][1] == Approx( 8.26) );
  REQUIRE( l_max1[2][1] == Approx( 6.61) );
}

TEST_CASE( "Sub-cell kernels: DG extrema.", "[subCellKernels][dgExtrema]" ) {
  // dg dofs
  double l_dgDofs1[2][4][2] = { { { 1.0, -1.0}, {-0.1,  0.1}, {-2.0,  2.0}, {3.0, -3.0} },
                                { {-0.3,  0.3}, { 0.4, -0.4}, { 2.2, -2.2}, {8.0, -8.0} } };

  // scatter operator
  double l_scatter1[4][9] = { { 1.0,  0.1,  1.0,   2.0,  5.0,  2.0,   3.0,  1.0,  0.0},
                              {-1.0, -0.1, -1.0,  -2.0, -5.0, -2.0,  -3.0, -1.0, -0.0},
                              { 0.2,  0.3,  0.4,   0.5,  0.6,  0.7,   0.8,  0.9,  1.0},
                              {-1.0, -0.9, -0.8,  -0.7, -0.6, -0.5,  -0.4, -0.3, -0.2} };

  double l_min1[2][2], l_max1[2][2];

  double l_scratch1[2][9][2];

  // determine extrema through sub-cell resolution
  edge::sc::Kernels< QUAD4R, 2, 2, 2 >::dgExtremaVanilla( l_dgDofs1, l_scatter1, l_scratch1, l_min1, l_max1 );

  // check results
  REQUIRE( l_min1[0][0] == Approx(-3.19) );
  REQUIRE( l_min1[1][0] == Approx(-8.26) );

  REQUIRE( l_min1[0][1] == Approx( -2.5) );
  REQUIRE( l_min1[1][1] == Approx( -0.6) );

  REQUIRE( l_max1[0][0] == Approx(  2.5) );
  REQUIRE( l_max1[1][0] == Approx(  0.6) );

  REQUIRE( l_max1[0][1] == Approx( 3.19) );
  REQUIRE( l_max1[1][1] == Approx( 8.26) );
}

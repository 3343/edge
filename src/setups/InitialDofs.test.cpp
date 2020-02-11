/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Unit tests for the initial DOFs.
 **/
#include <catch.hpp>
#include "dg/Basis.h"
#define private public
#include "InitialDofs.hpp"
#undef private
#include <cmath>

// hardcoded due to EDGEpre-dependency
#if PP_PRECISION == 64 && defined(PP_T_ELEMENTS_TRIA3) && PP_ORDER == 4
TEST_CASE( "Initial Dofs: DG modes.", "[initialDofs][dg]" ) {
  /*
   *        v0             v1
   *   -2 - *****************
   *      | * *     e0    * *
   *      | *   *       *   *
   *      | * e2  * v2*  e1 *
   *      | *       *       *
   *   -3 - *****************
   *        v3             v4
   *
   *        |---------------|
   *        1               2
   *
   *    ve | x   |  y
   *     0 | 1   | -2
   *     1 | 2   | -2
   *     2 | 1.5 | -3
   *     3 | 1   | -3
   *     4 | 2   | -3
   *
   *    el  | ve0 | ve1 | ve2
   *     0  |  2  |  1  |  0
   *     1  |  2  |  4  |  1
   *     2  |  0  |  3  |  2
   *
   *  ref pt (0.25, 0.25):
   *
   *    el |   x    |   y
   *     0 |  1.5   | -2.5
   *     1 |  1.75  | -2.75
   *     2 |  1.125 | -2.5
   */

  struct {
    double coords[2];
  } l_veChars1[5] = { {{ 1.0, -2.0}},
                      {{ 2.0, -2.0}},
                      {{ 1.5, -3.0}},
                      {{ 1.0, -3.0}},
                      {{ 2.0, -3.0}} };

  int l_elVe[3][3] = { {2, 1, 0},
                       {2, 4, 1},
                       {0, 3, 2} };

  // basis
  edge::dg::Basis l_basis1( TRIA3, 4);

  std::string l_exprs1[2] = { "q[0] := x;\
                               q[1] := y;\
                               q[2] := x+y;",

                              "q[0] := x*y;\
                               // missing init -> zero\n\
                               q[2] := x*x*y + 2.5;" };

  // compute dofs
  double l_dofs1[3][3][10][2];
  edge::setups::InitialDofs< TRIA3, 4, 3, 2>::dg( 0,
                                                  2,
                                                  l_exprs1,
                                                  l_basis1,
                                                  l_elVe,
                                                  l_veChars1,
                                                  l_dofs1 );

  edge::setups::InitialDofs< TRIA3, 4, 3, 2>::dg( 2,
                                                  1,
                                                  l_exprs1,
                                                  l_basis1,
                                                  l_elVe,
                                                  l_veChars1,
                                                  l_dofs1 );

  // reordered DOFs, modes as fastest dimension
  double l_dofsR1[3][3][2][10];
  for( unsigned short l_el = 0; l_el < 3; l_el++ )
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ )
      for( unsigned short l_md = 0; l_md < 10; l_md++ )
        for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
          l_dofsR1[l_el][l_qt][l_cr][l_md] = l_dofs1[l_el][l_qt][l_md][l_cr];

  // check DOFs
  double l_ref1[3] = {0.25,0.25,0};
  double l_out1;


  double l_x1[3] = {1.5,   1.75, 1.125};
  double l_y1[3] = {-2.5, -2.75, -2.5};

  for( unsigned short l_el = 0; l_el < 3; l_el++ ) {
    // q0, c0
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][0][0] );
    REQUIRE( l_out1 == Approx( l_x1[l_el] ) );

    // q1, c0
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][1][0] );
    REQUIRE( l_out1 == Approx( l_y1[l_el] ) );

    // q2, c0
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][2][0] );
    REQUIRE( l_out1 == Approx( l_x1[l_el] + l_y1[l_el] ) );


    // q0, c1
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][0][1] );
    REQUIRE( l_out1 == Approx( l_x1[l_el] * l_y1[l_el] ) );

    // q1, c1
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][1][1] );
    REQUIRE( l_out1 == Approx(0.0) );

    // q2, c1
    l_out1  = l_basis1.modal2ptval( l_ref1,
                                    l_dofsR1[l_el][2][1] );
    REQUIRE( l_out1 == Approx( l_x1[l_el] * l_x1[l_el] * l_y1[l_el] + 2.5 ) );
  }
}
#endif
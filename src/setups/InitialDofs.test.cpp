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
 * Unit tests for the initial DOFs.
 **/

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64

#include <catch.hpp>
#include "dg/Basis.h"
#define private public
#include "InitialDofs.hpp"
#undef private
#include <cmath>

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

TEST_CASE( "Initial Dofs: sub-cells.", "[initialDofs][sc]" ) {
  /* DG order: 4
   * #sub-vertices: 8
   * #sub-cells: 7
   *
   *  ve2        ve0   ve3   ve1        ve4   ve5
   *   |----------|-----|-----|----------|-----|
   *  -2    el1   0 el0 1 el2 2    el3   4 el4 5
   *
   *    ve |  x
   *     0 |  0
   *     1 |  2
   *     2 | -2
   *     3 |  1
   *     4 |  4
   *     5 |  5
   *
   *    el  | ve0 | ve1 | limited?
   *     0  |  0  |  3  | -
   *     1  |  0  |  2  | x
   *     2  |  3  |  1  | x
   *     3  |  4  |  2  | -
   *     4  |  4  |  5  | x
   */

  short l_elVe1[5][2] = { {0, 3},
                          {0, 2},
                          {3, 1},
                          {4, 1},
                          {4, 5} };

  short l_scSv1[7][2] = { { 1,2 },
                          { 2,3 },
                          { 3,4 },
                          { 4,5 },
                          { 5,6 },
                          { 0,1 },
                          { 6,7 } };

  typedef struct {
    double coords[1];
  } t_veChars1;

  t_veChars1 l_veChars1[6] = { {{ 0.0}},
                               {{ 2.0}},
                               {{-2.0}},
                               {{ 1.0}},
                               {{ 4.0}},
                               {{ 5.0}} };

  t_veChars1 l_svChars1[8] = { {{ 0.0                  }},
                               {{ 0.14285714285714285  }},
                               {{ 0.2857142857142857   }},
                               {{ 0.42857142857142855  }},
                               {{ 0.5714285714285714   }},
                               {{ 0.7142857142857143   }},
                               {{ 0.8571428571428571   }},
                               {{ 1.0                  }}
                             };

  typedef struct {
    int spType;
  } t_elChars1;

  t_elChars1 l_elChars[5] = { { 1},
                              {17},
                              {17},
                              { 1},
                              {17} };

  std::string l_exprs1[2] = { "if( x < -0.75 ) {\
                                q[0] := x + 1;\
                                q[1] := 2*x + 2;\
                                q[2] := x*x;\
                               }\
                               else if( x < 1.5 ) {\
                                q[0] := 0;\
                                q[1] := x*x;\
                                // missing init -> zero\n\
                               }\
                               else {\
                                q[0] := x*x*x;\
                                q[1] := 2.5;\
                                q[2] := 1.0/x;\
                               };",

                              "if( x < -1.9 ) {\
                                 q[0] := 4.2*x;\
                                 // missing init -> zero\n\
                                 q[2] := x*x*x;\
                               }\
                               else if( x < 4.2 ) {\
                                 q[0] := 1;\
                                 q[1] := sin(x);\
                                 q[2] := sin(x)+abs( cos(x) );\
                               }\
                               else {\
                                 // missing init -> zero\n\
                                 // missing init -> zero\n\
                                 q[2] := sqrt(x);\
                               }" };

  double l_dofsSc[3][3][7][2];

  edge::setups::InitialDofs< LINE, 4, 3, 2>::sc( (short) 0,
                                                 (short) 5,
                                                 (short) 0,
                                                 17,
                                                 l_exprs1,
                                                 l_elVe1,
                                                 l_scSv1,
                                                 l_veChars1,
                                                 l_svChars1,
                                                 l_elChars,
                                                 l_dofsSc );

   /*
    * check solution by checking some examples going from left to right
    * Remark: sub-cell order might be right-to-left due to DG-vertices
    */
   // vertex coordinates
   double l_vc[2];

   // check el1's left sub-cell DOFs (discontinuous in fused sim #2)
   double l_sh = 2.0 / 7.0;

   l_vc[0] = -2.0;
   l_vc[1] = -2.0 + l_sh;

   REQUIRE( l_dofsSc[0][0][6][0] == Approx( (   (l_vc[0]+1)
                                              + (l_vc[1]+1)
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][1][6][0] == Approx( (   (2.0*l_vc[0]+2)
                                              + (2.0*l_vc[1]+2)
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][2][6][0] == Approx( (   (l_vc[0]*l_vc[0])
                                              + (l_vc[1]*l_vc[1])
                                            ) *0.5
                                          ) );


   REQUIRE( l_dofsSc[0][0][6][1] == Approx( (   (4.2*l_vc[0])
                                              + (1)
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][1][6][1] == Approx( (   (0.0)
                                              + (sin(l_vc[1]))
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][2][6][1] == Approx( (   (l_vc[0]*l_vc[0]*l_vc[0])
                                              + (sin(l_vc[1])+std::abs( cos(l_vc[1]) ) )
                                            ) *0.5
                                          ) );

   // check el1's fifth sub-cell from the left (discontinuous in fused sim #1)
   //
   //   6    4    3    2     1    0   5
   // |----|----|----|----|xxxxx|----|----|
   //
   l_vc[0] = -2.0 + 4*l_sh;
   l_vc[1] = -2.0 + 5*l_sh;

   REQUIRE( l_dofsSc[0][0][1][0] == Approx( (   (l_vc[0]+1)
                                              + (0)
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][1][1][0] == Approx( (   (2.0*l_vc[0]+2)
                                              + (l_vc[1]*l_vc[1])
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][2][1][0] == Approx( (   (l_vc[0]*l_vc[0])
                                              + (0)
                                            ) *0.5
                                          ) );


   REQUIRE( l_dofsSc[0][0][1][1] == Approx( (   1.0
                                              + 1.0
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][1][1][1] == Approx( (   sin(l_vc[0])
                                              + sin(l_vc[1])
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[0][2][1][1] == Approx( (   (sin(l_vc[0])+std::abs( cos(l_vc[0]) ) )
                                              + (sin(l_vc[1])+std::abs( cos(l_vc[1]) ) )
                                            ) *0.5
                                          ) );

   // check el2's fourth sub-cell from the left (discontinuous in fused sim#1)
   //
   //   6    4    3    2     1    0   5
   // |----|----|----|xxxx|----|----|----|
   //
   l_sh = 1.0 / 7.0;

   l_vc[0] = 1.0 + l_sh * 3;
   l_vc[1] = 1.0 + l_sh * 4;

   REQUIRE( l_dofsSc[1][0][2][0] == Approx( (   (0)
                                              + (l_vc[1]*l_vc[1]*l_vc[1])
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[1][1][2][0] == Approx( (   (l_vc[0]*l_vc[0])
                                              + (2.5)
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[1][2][2][0] == Approx( (   (0)
                                              + (1.0/l_vc[1])
                                            ) *0.5
                                          ) );


   REQUIRE( l_dofsSc[1][0][2][1] == Approx( (   1
                                              + 1
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[1][1][2][1] == Approx( (   sin(l_vc[0])
                                              + sin(l_vc[1])
                                            ) *0.5
                                          ) );

   REQUIRE( l_dofsSc[1][2][2][1] == Approx( (   (sin(l_vc[0])+std::abs( cos(l_vc[0]) ) )
                                              + (sin(l_vc[1])+std::abs( cos(l_vc[1]) ) )
                                            ) *0.5
                                          ) );
}

TEST_CASE( "Initial Dofs: errors.", "[initialDofs][err]" ) {
  /* DG order: 2
   * #sub-vertices: 8
   * #sub-cells: 7
   *
   *   ve2       ve0   ve3     ve1       ve4    ve5
   *    |---------|------|------|---------|------|
   *   -2    el1  0  el0 1  el2 2   el3   4  el4 5
   *
   *
   * solution (first sim):
   *
   * 2 -|      ***|      *    **|          |      *  <-- limited elements with a
   *    |         |    * |      |          |    * |      valid DG-solution
   *    |         |  *   |      |          |  *   |
   * 1 -|***      |*     |****  |*         |*     |
   *    |         |      |      |    *     |      |
   *    |         |      |      |        * |      |
   * 0 -|---***---|------|------|----------*------|
   *
   * Three quantities shifted by 0.0, 1.0, 2.0
   *
   * solution (second sim): mirror at x-axis (multiplied by -1).
   *
   *    ve |  x
   *     0 |  0
   *     1 |  2
   *     2 | -2
   *     3 |  1
   *     4 |  4
   *     5 |  5
   *
   *    el  | ve0 | ve1 | limited?
   *     0  |  0  |  3  | -
   *     1  |  0  |  2  | x
   *     2  |  3  |  1  | x
   *     3  |  4  |  2  | -
   *     4  |  4  |  5  | x
   */

  short l_elVe1[5][2] = { {0, 3},
                          {0, 2},
                          {3, 1},
                          {4, 1},
                          {4, 5} };

  short l_scSv1[3][2] = { { 1,2 },
                          { 0,1 },
                          { 2,3 } };

  typedef struct {
    double coords[1];
  } t_veChars1;

  t_veChars1 l_veChars1[6] = { {{ 0.0}},
                               {{ 2.0}},
                               {{-2.0}},
                               {{ 1.0}},
                               {{ 4.0}},
                               {{ 5.0}} };

  t_veChars1 l_svChars1[8] = { {{ 0.0                  }},
                               {{ 0.33333333333333333  }},
                               {{ 0.66666666666666666  }},
                               {{ 1.0                  }}
                             };

  typedef struct {
    int spType;
  } t_elChars1;

  t_elChars1 l_elChars1[5] = { { 1},
                               {17},
                               {17},
                               { 1},
                               {17} };

  // DG solution
  double l_dofsDg1[5][3][2][2];

  // sub-cell solution
  double l_dofsSc1[3][3][3][2];

  // init everything with invalid values
  for( unsigned short l_el = 0; l_el < 5; l_el++ )
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ )
      for( unsigned short l_md = 0; l_md < 2; l_md++ )
        for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
          l_dofsDg1[l_el][l_qt][l_md][l_cr] = std::numeric_limits< double >::max();

  for( unsigned short l_li = 0; l_li < 3; l_li++ )
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ )
      for( unsigned short l_sc = 0; l_sc < 3; l_sc++ )
        for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
          l_dofsSc1[l_li][l_qt][l_sc][l_cr] = std::numeric_limits< double >::max();

  // DG, el0, q0, cr0
  l_dofsDg1[0][0][0][0] = 1.5;
  l_dofsDg1[0][0][1][0] = 0.5;

  // DG, el3, q0, cr0
  l_dofsDg1[3][0][0][0] = 0.5;
  l_dofsDg1[3][0][1][0] = 0.5; // right to left: *-1

  // DG, el4, q0, cr0
  l_dofsDg1[4][0][0][0] = 1.5;
  l_dofsDg1[4][0][1][0] = 0.5;

  // set the remaining quantities
  for( unsigned short l_qt = 1; l_qt < 3; l_qt++ ) {
    l_dofsDg1[0][l_qt][0][0] = l_dofsDg1[0][0][0][0] + l_qt;
    l_dofsDg1[0][l_qt][1][0] = l_dofsDg1[0][0][1][0];

    l_dofsDg1[3][l_qt][0][0] = l_dofsDg1[3][0][0][0] + l_qt;
    l_dofsDg1[3][l_qt][1][0] = l_dofsDg1[3][0][1][0];

    l_dofsDg1[4][l_qt][0][0] = l_dofsDg1[4][0][0][0] + l_qt;
    l_dofsDg1[4][l_qt][1][0] = l_dofsDg1[4][0][1][0];
  }

  // set second run
  for( unsigned short l_el = 0; l_el < 5; l_el++ )
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ )
      for( unsigned short l_md = 0; l_md < 2; l_md++ )
        l_dofsDg1[l_el][l_qt][l_md][1] = -l_dofsDg1[l_el][l_qt][l_md][0];

  // SC, sc0, q0, cr0 (right-to-left)
  l_dofsSc1[0][0][0][0] = 0.0; // inner
  l_dofsSc1[0][0][1][0] = 2.0; // "left"
  l_dofsSc1[0][0][2][0] = 1.0; // "right"

  // SC, sc1, q0, cr0
  l_dofsSc1[1][0][0][0] = 1.0;
  l_dofsSc1[1][0][1][0] = 1.0;
  l_dofsSc1[1][0][2][0] = 2.0;

  // set remaining quantities
  for( unsigned short l_li = 0; l_li < 2; l_li++ )
    for( unsigned short l_qt = 1; l_qt < 3; l_qt++ )
      for( unsigned short l_sc = 0; l_sc < 3; l_sc++ )
        l_dofsSc1[l_li][l_qt][l_sc][0] = l_dofsSc1[l_li][0][l_sc][0] + l_qt;

  // set second fused run
  for( unsigned short l_li = 0; l_li < 2; l_li++ )
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ )
      for( unsigned short l_sc = 0; l_sc < 3; l_sc++ )
        l_dofsSc1[l_li][l_qt][l_sc][1] = -l_dofsSc1[l_li][l_qt][l_sc][0];

  // admissibility of the limited solution
  bool l_adm1[3][2] = { {false, false}, {false, false}, {true, true} };

  // reference solution which matches the initial numerical solution exactly
  std::string l_exprs1[2];
  // first fused simulation
  l_exprs1[0] = "// el1 \n\
                 if(      x < -2.0 + 2.0*(1.0/3.0) ) {\
                   q[0] := 1.0;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 else if( x < -2.0 + 2.0*(2.0/3.0) ) {\
                   q[0] := 0.0;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 else if( x <  0.0 ) {\
                   q[0] := 2.0;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 // el0 \n\
                 else if( x < 1.0 ) {\
                   q[0] := 1.0 + x;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 // el2 \n\
                 else if( x < 1.0 + 1.0*(2.0/3.0) ) {\
                   q[0] := 1.0;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 else if( x < 2.0 ) {\
                   q[0] := 2.0;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 // el3 \n\
                 else if( x < 4.0 ) {\
                   q[0] := 1.0 - (x-2.0)*0.5;\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 }\
                 // el4 \n\
                 else if( x < 5.0 ) {\
                   q[0] := 1.0 + (x-4.0);\
                   q[1] := q[0] + 1.0;\
                   q[2] := q[0] + 2.0;\
                 };\n";
  // second fused simulation
  l_exprs1[1] = l_exprs1[0] + "q[0] := -q[0];\
                               q[1] := -q[1];\
                               q[2] := -q[2];";

  // L1 error
  double l_l1[3][2]   = { {-1, -1},
                          {-1, -1},
                          {-1, -1} };
  // L2 error
  double l_l2p2[3][2]   = { {-1, -1},
                            {-1, -1},
                            {-1, -1} };
  // Linf error
  double l_lInf[3][2] = { {-1, -1},
                          {-1, -1},
                          {-1, -1} };

  // DG basis
  edge::dg::Basis l_basis1( LINE, 2);

  edge::setups::InitialDofs< LINE, 2, 3, 2>::err( (short) 0,
                                                  (short) 5,
                                                  (short) 0,
                                                  17,
                                                  l_exprs1,
                                                  l_basis1,
                                                  l_elVe1,
                                                  l_scSv1,
                                                  l_veChars1,
                                                  l_svChars1,
                                                  l_elChars1,
                                                  l_adm1,
                                                  l_dofsDg1,
                                                  l_dofsSc1,
                                                  l_l1,
                                                  l_l2p2,
                                                  l_lInf );

  // check for errors of matching reference solution
  for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
      REQUIRE( l_l1[l_qt][l_cr]   == Approx(0.0) );
      REQUIRE( l_l2p2[l_qt][l_cr] == Approx(0.0) );
      REQUIRE( l_lInf[l_qt][l_cr] == Approx(0.0) );
    }
  }

  // modify numerical solution
  l_dofsSc1[0][0][0][0] = -1.0; // qp1 off by -1
  l_dofsSc1[1][0][2][0] =  0.5; // qp2 off by 1.5
  l_dofsDg1[0][0][0][0] =  0.0; // alls qps off by -1.5

  l_dofsSc1[0][1][0][0] = -2.0; // qp1 off by -3
  l_dofsSc1[1][1][2][0] = -0.5; // qp2 off by -3.5
  l_dofsDg1[0][1][0][0] = -1.0; // all qps off by -3.5

  l_dofsSc1[0][2][0][1] =  3.0; // qp1 off by 5.0
  l_dofsSc1[1][2][2][1] =  2.5; // qp1 off by 6.5
  l_dofsDg1[0][2][0][1] =  2.0; // alls qps off by 5.5

  // manually compute errors through gauss legendre quadrature and max
  double l_l1Ut1[3][2] = { {   1.0 * (2.0 * 0.444444444444444)  // sc0
                             + 1.5 * (1.0 * 0.277777777777777)  // sc1
                             + 1.5 * (1.0 * 1.000000000000000), // dg0
                             0.0 },
                           {   3.0 * (2.0 * 0.444444444444444)  // sc0
                             + 3.5 * (1.0 * 0.277777777777777)  // sc1
                             + 3.5 * (1.0 * 1.000000000000000), // dg0
                             0.0 },
                           { 0.0,
                               5.0 * (2.0 * 0.444444444444444)  // sc0
                             + 6.5 * (1.0 * 0.277777777777777)  // sc1
                             + 5.5 * (1.0 * 1.000000000000000), // dg0
                           } };

  double l_l2p2Ut1[3][2] = { {   1.0 * 1.0 * (2.0 * 0.444444444444444)  // sc0
                               + 1.5 * 1.5 * (1.0 * 0.277777777777777)  // sc1
                               + 1.5 * 1.5 * (1.0 * 1.000000000000000), // dg0
                               0.0 },
                             {   3.0 * 3.0 * (2.0 * 0.444444444444444)  // sc0
                               + 3.5 * 3.5 * (1.0 * 0.277777777777777)  // sc1
                               + 3.5 * 3.5 * (1.0 * 1.000000000000000), // dg0
                               0.0 },
                             { 0.0,
                                 5.0 * 5.0 * (2.0 * 0.444444444444444)  // sc0
                               + 6.5 * 6.5 * (1.0 * 0.277777777777777)  // sc1
                               + 5.5 * 5.5 * (1.0 * 1.000000000000000), // dg0
                             } };

  double l_lInfUt1[3][2] = { { 1.5, 0.0 },
                             { 3.5, 0.0 },
                             { 0.0, 6.5 } };

  edge::setups::InitialDofs< LINE, 2, 3, 2>::err( (short) 0,
                                                  (short) 5,
                                                  (short) 0,
                                                  17,
                                                  l_exprs1,
                                                  l_basis1,
                                                  l_elVe1,
                                                  l_scSv1,
                                                  l_veChars1,
                                                  l_svChars1,
                                                  l_elChars1,
                                                  l_adm1,
                                                  l_dofsDg1,
                                                  l_dofsSc1,
                                                  l_l1,
                                                  l_l2p2,
                                                  l_lInf );

  // check errors
  for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
      REQUIRE( l_l1[l_qt][l_cr] == Approx(l_l1Ut1[l_qt][l_cr]) );
      REQUIRE( l_l2p2[l_qt][l_cr] == Approx(l_l2p2Ut1[l_qt][l_cr]) );
      REQUIRE( l_lInf[l_qt][l_cr] == Approx(l_lInfUt1[l_qt][l_cr]) );
    }
  }
}

#endif

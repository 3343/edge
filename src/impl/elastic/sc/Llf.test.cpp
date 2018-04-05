/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Unit tests of elastic local Lax-Friedrichs sub-cell solver.
 **/
#include <catch.hpp>
#define private public
#include "Llf.hpp"
#undef private

TEST_CASE( "Sub-cell elastic local Lax-Friedrichs: Initialization", "[elasticScLlf][init]" ) {
  /*
   * Our example problem:
   *
   * lp: limited elements + adjacent elements
   * el: dense elements
   * fa: face
   * elFaEl: face-adjacent elements
   *
   *   lp | el | vol
   *   0  | 4  | 5.0
   *   1  | 2  | 2.0
   *   2  | 3  | 3.0
   *   3  | 1  | 4.0
   *
   *   el | fa0 | fa1 | fa2
   *   0  | max | max | max
   *   1  | 0   | 1   | 0
   *   2  | 0   | 0   | 1
   *   3  | 1   | 1   | 0
   *   4  | 0   | 0   | 1
   *
   *   fa | out-x      | out-y     | out-z, t0: x,y,z, t1: x,y,z | face area
   *   0  |  1         |  0        | max                         | 0.25
   *   1  | -1         |  0        | max                         | 0.5
   *
   *   el | elFaEl0 | elFaEl1 | elFaEl2
   *   0  | max     | max     | max
   *   1  | 3       | 4       | 2
   *   2  | 4       | 1       | 3
   *   3  | 4       | 1       | 2
   *   4  | 1       | 3       | 2
   *
   *   el | lam     | mu      | rho
   *   0  | -1      | -1      | -1
   *   1  | 20.8E9  | 10.4E9  | 2600
   *   2  | 32.4E9  | 32.4E9  | 2700
   *   3  | 20.8E9  | 10.4E9  | 2600
   *   4  | 32.4E9  | 32.4E9  | 2700
   */

  const double l_dMax = std::numeric_limits< double >::max();
  struct {
    double outNormal[3];
    double tangent0[3];
    double tangent1[3];
    double area;
  } l_charsFa[2] = { { {    1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.25 },
                     { {   -1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.5 } };
  
  struct {
    double volume;
  } l_charsEl[5] = { {l_dMax}, {4.0}, {2.0}, {3.0}, {5.0} };

  struct {
    double lam;
    double mu;
    double rho;
  } l_matPars[5] = { {l_dMax, l_dMax, l_dMax},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0} };

  const int l_iMax = std::numeric_limits< int >::max();
  int l_lpEl[4]  = {4, 2, 3, 1};
  int l_elFa[5][3] = { {l_iMax, l_iMax, l_iMax},
                       {0, 1, 0},
                       {0, 0, 1},
                       {1, 1, 0},
                       {0, 0, 1} };
  int l_elFaEl[5][3] = { {l_iMax, l_iMax, l_iMax},
                         {3, 4, 2},
                         {4, 1, 3},
                         {4, 1, 2},
                         {1, 3, 2} };

  // dummy data: triangles don't have additional sub-faces
  int l_elVe[1][3] = { {0,0,0} };
  struct {
    double coords[3];
  } l_charsVe[1] = { { {0,0,0} } };

  // wrapper for the LLF solvers
  edge::elastic::sc::Llf< double,
                          TRIA3,
                          4,
                          8 > l_llf;

  // dynamic memory allocations
  edge::data::Dynamic l_dynMem;

  // allocate memory for solvers of 4 limited plus elements
  l_llf.alloc( 4,
               l_dynMem );

  // init the solvers
  l_llf.init( 0,
              4,
              l_lpEl,
              l_elVe,
              l_elFa,
              l_elFaEl,
              l_charsVe,
              l_charsFa,
              l_charsEl,
              l_matPars );

  /* Jacobian in x-direction:
   *   0      0 0       -(lam + 2.0 * mu) 0
   *   0      0 0              -lam       0
   *   0      0 0                0       -mu
   *   -1/rho 0 0                0        0
   *   0      0 -1/rho           0        0
   */
  // first jacobian
  double l_jac0[5][5] = { { 0,        0,  0,       -(20.8E9 + 2.0 * 10.4E9),  0      },
                          { 0,        0,  0,        -20.8E9,                  0      },
                          { 0,        0,  0,         0,                      -10.4E9 },
                          {-1.0/2600, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2600,  0,                       0      } };

  // second jacobian
  double l_jac1[5][5] = { { 0,        0,  0,       -(32.4E9 + 2.0 * 32.4E9),  0      },
                          { 0,        0,  0,        -32.4E9,                  0      },
                          { 0,        0,  0,         0,                      -32.4E9 },
                          {-1.0/2700, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2700,  0,                       0      } };

  // LLF contributions are averaged + left element's flux contribution
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_jac0[l_q1][l_q2] *= -0.5;
      l_jac1[l_q1][l_q2] *= -0.5;
    }
  }

  // viscous terms, +1: left element's viscosity contribution
  double l_vis0 = -0.5 * std::sqrt( (20.8E9 + 2 * 10.4E9)/2600 );
  double l_vis1 = -0.5 * std::sqrt( (32.4E9 + 2 * 32.4E9)/2700 );


  // scaling of sub-faces and sub-cells, first limited plus element, first face
  double l_sca = (0.25 / 7.0) * (49.0 / 5.0);

  // check viscosity and central flux contribution for first limited element, first face (normal points inside the element)
  REQUIRE( l_llf.m_vis[0][0  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[0][0+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[0][0  ][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[0][0+3][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, first limited plus element, second face
  l_sca = (0.25 / 7.0) * (49.0 / 5.0);

  // check viscosity and central flux contribution for first limited element, second face (normal points inside the element)
  REQUIRE( l_llf.m_vis[0][1  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[0][1+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[0][1  ][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[0][1+3][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, first limited plus element, third face
  l_sca = (0.5 / 7.0) * (49.0 / 5.0);

  // check viscosity and central flux contribution for first limited element, third face (normal points inside the element)
  REQUIRE( l_llf.m_vis[0][2  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[0][2+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[0][2  ][l_q1][l_q2] == Approx( l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[0][2+3][l_q1][l_q2] == Approx( l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }


  // scaling of sub-faces and sub-cells, second limited plus element, first face
  l_sca = (0.25 / 7.0) * (49.0 / 2.0);

  // check viscosity and central flux contribution for second limited element, first face (normal points out the element)
  REQUIRE( l_llf.m_vis[1][0  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[1][0+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[1][0  ][l_q1][l_q2] == Approx( l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[1][0+3][l_q1][l_q2] == Approx( l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, second limited plus element, second face
  l_sca = (0.25 / 7.0) * (49.0 / 2.0);

  // check viscosity and central flux contribution for second limited element, second face (normal points inside the element)
  REQUIRE( l_llf.m_vis[1][1  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[1][1+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[1][1  ][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[1][1+3][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, second limited plus element, third face
  l_sca = (0.5 / 7.0) * (49.0 / 2.0);

  // check viscosity and central flux contribution for second limited element, third face (normal points out the element)
  REQUIRE( l_llf.m_vis[1][2  ] == Approx( l_vis1*l_sca ) );
  REQUIRE( l_llf.m_vis[1][2+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[1][2  ][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[1][2+3][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }


  // scaling of sub-faces and sub-cells, third limited plus element, first face
  l_sca = (0.5 / 7.0) * (49.0 / 3.0);

  // check viscosity and central flux contribution for third limited element, first face (normal points out the element)
  REQUIRE( l_llf.m_vis[2][0  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[2][0+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[2][0  ][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[2][0+3][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, third limited plus element, second face
  l_sca = (0.5 / 7.0) * (49.0 / 3.0);

  // check viscosity and central flux contribution for third limited element, second face (normal points inside the element)
  REQUIRE( l_llf.m_vis[2][1  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[2][1+3] == Approx( l_vis0*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[2][1  ][l_q1][l_q2] == Approx( l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[2][1+3][l_q1][l_q2] == Approx( l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, third limited plus element, third face
  l_sca = (0.25 / 7.0) * (49.0 / 3.0);

  // check viscosity and central flux contribution for third limited element, third face (normal points inside the element)
  REQUIRE( l_llf.m_vis[2][2  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[2][2+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[2][2  ][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[2][2+3][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }


  // scaling of sub-faces and sub-cells, fourth limited plus element, first face
  l_sca = (0.25 / 7.0) * (49.0 / 4.0);

  // check viscosity and central flux contribution for fourth limited element, first face (normal points out the element)
  REQUIRE( l_llf.m_vis[3][0  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[3][0+3] == Approx( l_vis0*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[3][0  ][l_q1][l_q2] == Approx( l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[3][0+3][l_q1][l_q2] == Approx( l_jac0[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, fourth limited plus element, second face
  l_sca = (0.5 / 7.0) * (49.0 / 4.0);

  // check viscosity and central flux contribution for fourth limited element, second face (normal points out the element)
  REQUIRE( l_llf.m_vis[3][1  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[3][1+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[3][1  ][l_q1][l_q2] == Approx( -l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[3][1+3][l_q1][l_q2] == Approx( -l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }

  // scaling of sub-faces and sub-cells, fourth limited plus element, third face
  l_sca = (0.25 / 7.0) * (49.0 / 4.0);

  // check viscosity and central flux contribution for fourth limited element, third face (normal points out the element)
  REQUIRE( l_llf.m_vis[3][2  ] == Approx( l_vis0*l_sca ) );
  REQUIRE( l_llf.m_vis[3][2+3] == Approx( l_vis1*l_sca ) );

  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      REQUIRE( l_llf.m_cen[3][2  ][l_q1][l_q2] == Approx( l_jac0[l_q1][l_q2]*l_sca ) );
      REQUIRE( l_llf.m_cen[3][2+3][l_q1][l_q2] == Approx( l_jac1[l_q1][l_q2]*l_sca ) );
    }
  }
}

TEST_CASE( "Sub-cell elastic local Lax-Friedrichs: Is of the LLF solvers", "[elasticScLlf][llfIds]" ) {
  /* Tria3 elements, 2nd order:
   *   0,  1,  2: DG, left
   *   3,  4,  5: DG, right
   *   6,  7,  8: SC, left
   *   9, 10, 11: SC, right
   */

  unsigned short l_llfIds[2];
  bool l_left;
  
  //
  // heterogeneous
  //
  // type 0
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 0 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 3     );

  // type 1
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 1 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 4     );

  // type 2
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 2 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 5     );

  // type 3
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 3 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 3     );

  // type 4
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 4 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 4     );

  // type 5
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 5 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 5     );

  //
  // homogeneous
  //
  // type 0
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 6 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 0     );

  // type 1
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 7 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 1     );

  // type 2
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 8 , l_llfIds );
  
  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 2     );

  // type 3
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds( 9 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 0     );

  // type 4
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds(10 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 1     );

  // type 5
  l_left = edge::elastic::sc::Llf< double,
                                   TRIA3,
                                   2,
                                   3 >::llfIds(11 , l_llfIds );
  
  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 2     );

  /* Tet4 elements, 2nd order:
   *    0,  1,  2,  3        : DG-surf, left
   *    4,  5,  6,  7        : DG-surf, right
   *    8,  9, 10, 11, 12, 13: SC, left
   *   14, 15, 16, 17, 18, 19: SC, right
   */
  //
  // heterogeneous
  //
  // type 0
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 0, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 6     );

  // type 1
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 1, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 7     );

  // type 2
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 2, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 8     );

  // type 3
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 3, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 3     );
  REQUIRE( l_llfIds[1] == 9     );

  // type 4
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 4, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 6     );

  // type 5
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 5, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 7     );

  // type 6
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 6, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 8     );

  // type 7
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 7, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 3     );
  REQUIRE( l_llfIds[1] == 9     );

  //
  // homogeneous
  //
  // type 8
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 8, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 0     );

  // type 9
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 9, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 1     );

  // type 10
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 10, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 2     );

  // type 11
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 11, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 3     );
  REQUIRE( l_llfIds[1] == 3     );

  // type 12
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 12, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 4     );
  REQUIRE( l_llfIds[1] == 4     );

  // type 13
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 13, l_llfIds );

  REQUIRE( l_left      == true  );
  REQUIRE( l_llfIds[0] == 5     );
  REQUIRE( l_llfIds[1] == 5     );

  // type 14
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 14, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 0     );
  REQUIRE( l_llfIds[1] == 0     );

  // type 15
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 15, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 1     );
  REQUIRE( l_llfIds[1] == 1     );

  // type 16
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 16, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 2     );
  REQUIRE( l_llfIds[1] == 2     );

  // type 17
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 17, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 3     );
  REQUIRE( l_llfIds[1] == 3     );

  // type 18
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 18, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 4     );
  REQUIRE( l_llfIds[1] == 4     );

  // type 19
  l_left = edge::elastic::sc::Llf< double,
                                   TET4,
                                   2,
                                   3 >::llfIds( 19, l_llfIds );

  REQUIRE( l_left      == false );
  REQUIRE( l_llfIds[0] == 5     );
  REQUIRE( l_llfIds[1] == 5     );
}

TEST_CASE( "Sub-cell elastic local Lax-Friedrichs: Time Step", "[elasticScLlf][tsSc]" ) {
  /*
   * Setup: Analogue to initialization unit test with DOFs.
   *
   * Update of sub-cell DOFs for limited element 3.
   * Second order method: 3 sub-faces per DG-face, 9 sub-cell per DG-element.
   *
   * Sub-grid layout:
   *
   *     9
   *     | *
   *     |   *   14
   *  15 |     *
   *     |   7   *
   *     |         *
   *     7-----------8
   *     | *         | *
   *     |   *    2  |   *   13
   *  16 |     *     |     *
   *     |   8   *   |   6   *
   *     |         * |         *
   *     4-----------5-----------6
   *     | *         | *         | *
   *     |   *    0  |   *    1  |   *   12
   *  17 |     *     |     *     |     *
   *     |   3   *   |   4   *   |  5    *
   *     |         * |         * |         *
   *     0-----------1-----------2-----------3
   *           9           10          11
   *
   * Initial DOFs:
   *         xx              yy             xy               u        v           not default? | covered by unit test?
   *  sc |                |             |                |        |             |              |
   *   0 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              | yes
   *   1 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   2 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   3 |   2E9 7E9  2E9 | 6E9 2E9 6E9 | -100 -5E9 -100 |  6 2 6 | -30 -20 -30 | yes          |
   *   4 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   5 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   6 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              | yes
   *   7 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   8 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *   9 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  10 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  11 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  12 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  13 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  14 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  15 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  16 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   *  17 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20 |              |
   */

  // init sub-cell connectivity
  unsigned short l_scSfSc[18][3] = { {4,8,3},
                                     {5,6,4},
                                     {6,7,8},
                                     {9,0,17},
                                     {10,1,0},
                                     {11,12,1},
                                     {1,13,2},
                                     {2,14,15},
                                     {0,2,16},
                                     {3,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {4,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {5,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {5,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {6,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {7,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {7,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {8,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()},
                                     {3,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max()} };

  // init sub-face types
  unsigned short l_scTySf[9][3] = { {11,9,10},
                                    {11,9,10},
                                    {11,9,10},
                                    {0,7,2},
                                    {0,7,8},
                                    {0,1,8},
                                    {6,1,8},
                                    {6,1,2},
                                    {6,7,2} };

  // init sub-cell net-updates
  bool l_netUpSc[3] = { false, false, false };

  // init DOFs
  double l_scDofs[5][18][3];
  for( unsigned short l_sc = 0; l_sc < 18; l_sc++ ) {
    l_scDofs[0][l_sc][0] = 1E10;
    l_scDofs[0][l_sc][1] =  9E9;
    l_scDofs[0][l_sc][2] = 1E10;

    l_scDofs[1][l_sc][0] =  4E9;
    l_scDofs[1][l_sc][1] =  5E9;
    l_scDofs[1][l_sc][2] =  4E9;

    l_scDofs[2][l_sc][0] = -2E9;
    l_scDofs[2][l_sc][1] = -3E9;
    l_scDofs[2][l_sc][2] = -2E9;

    l_scDofs[3][l_sc][0] =    5;
    l_scDofs[3][l_sc][1] =    4;
    l_scDofs[3][l_sc][2] =    5;

    l_scDofs[4][l_sc][0] =  -20;
    l_scDofs[4][l_sc][1] =  -10;
    l_scDofs[4][l_sc][2] =  -20;
  }

  // modify sub-cell 3
  l_scDofs[0][3][0] =  2E9;
  l_scDofs[0][3][1] =  7E9;
  l_scDofs[0][3][2] =  2E9;

  l_scDofs[1][3][0] =  6E9;
  l_scDofs[1][3][1] =  2E9;
  l_scDofs[1][3][2] =  6E9;

  l_scDofs[2][3][0] = -100;
  l_scDofs[2][3][1] = -5E9;
  l_scDofs[2][3][2] = -100;

  l_scDofs[3][3][0] =    6;
  l_scDofs[3][3][1] =    2;
  l_scDofs[3][3][2] =    6;

  l_scDofs[4][3][0] =  -30;
  l_scDofs[4][3][1] =  -20;
  l_scDofs[4][3][2] =  -30;

  const double l_dMax = std::numeric_limits< double >::max();
  struct {
    double outNormal[3];
    double tangent0[3];
    double tangent1[3];
    double area;
  } l_charsFa[2] = { { {    1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.25 },
                     { {   -1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.5 } };
  
  struct {
    double volume;
  } l_charsEl[5] = { {l_dMax}, {4.0}, {2.0}, {3.0}, {5.0} };

  struct {
    double lam;
    double mu;
    double rho;
  } l_matPars[5] = { {l_dMax, l_dMax, l_dMax},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0} };

  const int l_iMax = std::numeric_limits< int >::max();
  int l_lpEl[4]  = {4, 2, 3, 1};
  int l_elFa[5][3] = { {l_iMax, l_iMax, l_iMax},
                       {0, 1, 0},
                       {0, 0, 1},
                       {1, 1, 0},
                       {0, 0, 1} };
  int l_elFaEl[5][3] = { {l_iMax, l_iMax, l_iMax},
                         {3, 4, 2},
                         {4, 1, 3},
                         {4, 1, 2},
                         {1, 3, 2} };

  // dummy data: triangles don't have additional sub-faces
  int l_elVe[1][3] = { {0,0,0} };
  struct {
    double coords[3];
  } l_charsVe[1] = { { {0,0,0} } };

  // wrapper for the LLF solvers
  edge::elastic::sc::Llf< double,
                          TRIA3,
                          2,
                          3 > l_llf;

  // dynamic memory allocations
  edge::data::Dynamic l_dynMem;

  // allocate memory for solvers of 4 limited plus elements
  l_llf.alloc( 4,
               l_dynMem );

  // init the solvers
  l_llf.init( 0,
              4,
              l_lpEl,
              l_elVe,
              l_elFa,
              l_elFaEl,
              l_charsVe,
              l_charsFa,
              l_charsEl,
              l_matPars );

  // resulting sub-cell dofs
  double l_scRes[5][9][3];

  l_llf.tsSc( 0.006,
              3,
              l_scSfSc,
              l_scTySf,
              l_netUpSc,
              l_scDofs,
              l_scRes );

  // check that fused sim 0 and 2 are identical
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ )
    for( unsigned short l_sc = 0; l_sc < 9; l_sc++ )
      REQUIRE( l_scRes[l_qt][l_sc][0] == l_scRes[l_qt][l_sc][2] );

  /*
   * check sub-cell 0:
   *   - homogeneous: 20.8E9  | 10.4E9  | 2600
   *   - DG-normals   (1,0), (-1,0), (1,0)
   *   - "DG-faces": 2, 0, 1
   *   - SC is on the right side of all sub-face: *(-1) for all updates
   *   - elements DOFs:
   *     - 0 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     - 4 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     - 8 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     - 3 |   2E9 7E9  2E9 | 6E9 2E9 6E9 | -100 -5E9 -100 |  6 2 6 | -30 -20 -30
   */

  // first jacobian
  double l_jac0[5][5] = { { 0,        0,  0,       -(20.8E9 + 2.0 * 10.4E9),  0      },
                          { 0,        0,  0,        -20.8E9,                  0      },
                          { 0,        0,  0,         0,                      -10.4E9 },
                          {-1.0/2600, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2600,  0,                       0      } };

  // second jacobian
  double l_jac1[5][5] = { { 0,        0,  0,       -(32.4E9 + 2.0 * 32.4E9),  0      },
                          { 0,        0,  0,        -32.4E9,                  0      },
                          { 0,        0,  0,         0,                      -32.4E9 },
                          {-1.0/2700, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2700,  0,                       0      } };

  // LLF contributions are averaged + left element's flux contribution
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_jac0[l_q1][l_q2] *= -0.5;
      l_jac1[l_q1][l_q2] *= -0.5;
    }
  }

  // viscous terms, +1: left element's viscosity cpntribution
  double l_vis0 = -0.5 * std::sqrt( (20.8E9 + 2 * 10.4E9)/2600 );
  double l_vis1 = -0.5 * std::sqrt( (32.4E9 + 2 * 32.4E9)/2700 );

  // reference solution
  double l_scUt1[5][3] = { { 1E10,  9E9, 1E10 },
                           {  4E9,  5E9,  4E9 },
                           { -2E9, -3E9, -2E9 },
                           {    5,    4,    5 },
                           {  -20,   -10, -20 } };

  // apply updates
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      l_scUt1[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt1[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][4][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

      l_scUt1[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt1[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][8][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

      l_scUt1[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][0][l_cr] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt1[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][3][l_cr] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
    }


    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
        l_scUt1[l_q1][l_cr] -=   l_scDofs[l_q2][0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt1[l_q1][l_cr] -=   l_scDofs[l_q2][4][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

        l_scUt1[l_q1][l_cr] -=   l_scDofs[l_q2][0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt1[l_q1][l_cr] -=   l_scDofs[l_q2][8][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

        l_scUt1[l_q1][l_cr] +=   l_scDofs[l_q2][0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt1[l_q1][l_cr] +=   l_scDofs[l_q2][3][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
      }
    }
  }

  // check the results
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_scUt1[l_qt][l_cr] == Approx( l_scRes[l_qt][0][l_cr] ) );
    }
  }


  /*
   * check sub-cell 6:
   *   - heterogeneous
   *     - own 20.8E9  | 10.4E9  | 2600
   *     - adj 32.4E9  | 32.4E9  | 2700
   *   - DG-normals   (1,0), (-1,0), (1,0)
   *   - "DG-faces": 0, 1, 2
   *   - SC is on the right side of all sub-face: *(-1) for all updates
   *   - elements DOFs:
   *     -  6 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     -  1 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     - 13 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *     -  2 |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   */

  // reference solution
  double l_scUt2[5][3] = { { 1E10,  9E9, 1E10 },
                           {  4E9,  5E9,  4E9 },
                           { -2E9, -3E9, -2E9 },
                           {    5,    4,    5 },
                           {  -20,   -10, -20 } };

  // apply updates
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      l_scUt2[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][ 0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt2[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][ 1][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

      l_scUt2[l_q1][l_cr] += l_vis1 * l_scDofs[l_q1][ 0][l_cr] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt2[l_q1][l_cr] -= l_vis1 * l_scDofs[l_q1][13][l_cr] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;

      l_scUt2[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][ 0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt2[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][ 2][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
    }


    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
        l_scUt2[l_q1][l_cr] +=   l_scDofs[l_q2][ 0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt2[l_q1][l_cr] +=   l_scDofs[l_q2][ 1][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

        l_scUt2[l_q1][l_cr] -=   l_scDofs[l_q2][ 0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt2[l_q1][l_cr] -=   l_scDofs[l_q2][13][l_cr]
                               * l_jac1[l_q1][l_q2] * (0.5  / 3.0) * (9.0 / 4.0) * 0.006;

        l_scUt2[l_q1][l_cr] +=   l_scDofs[l_q2][ 0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt2[l_q1][l_cr] +=   l_scDofs[l_q2][ 2][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      }
    }
  }

  // check the results
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_scUt2[l_qt][l_cr] == Approx( l_scRes[l_qt][6][l_cr] ) );
    }
  }



  // check the net-update computations for the first sub-cell (inner: same as before)
   l_netUpSc[0] = l_netUpSc[1] = l_netUpSc[2] = true;

  l_llf.tsSc( 0.006,
              3,
              l_scSfSc,
              l_scTySf,
              l_netUpSc,
              l_scDofs,
              l_scRes );

  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_scUt1[l_qt][l_cr] == Approx( l_scRes[l_qt][0][l_cr] ) );
    }
  }


  // check the net-update compuations for sub-cell 6
  // sub-cell 13 is at the DG-boundary -> DOFs are applied directly

  // reference solution
  double l_scUt3[5][3] = { { 1E10,  9E9, 1E10 },
                           {  4E9,  5E9,  4E9 },
                           { -2E9, -3E9, -2E9 },
                           {    5,    4,    5 },
                           {  -20,   -10, -20 } };

  // apply updates
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      l_scUt3[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][ 0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt3[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][ 1][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

      l_scUt3[l_q1][l_cr] += l_scDofs[l_q1][13][l_cr];

      l_scUt3[l_q1][l_cr] += l_vis0 * l_scDofs[l_q1][ 0][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      l_scUt3[l_q1][l_cr] -= l_vis0 * l_scDofs[l_q1][ 2][l_cr] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
    }


    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
        l_scUt3[l_q1][l_cr] +=   l_scDofs[l_q2][ 0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt3[l_q1][l_cr] +=   l_scDofs[l_q2][ 1][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;

        l_scUt3[l_q1][l_cr] +=   l_scDofs[l_q2][ 0][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
        l_scUt3[l_q1][l_cr] +=   l_scDofs[l_q2][ 2][l_cr]
                               * l_jac0[l_q1][l_q2] * (0.25 / 3.0) * (9.0 / 4.0) * 0.006;
      }
    }
  }

  // check the results
  for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_scUt3[l_qt][l_cr] == Approx( l_scRes[l_qt][6][l_cr] ) );
    }
  }

}

TEST_CASE( "Sub-cell elastic local Lax-Friedrichs: Net-updates at DG-face", "[elasticScLlf][nuFaSf]" ) {
  /*
   * Material parameters analogue to other settings.
   * We use the heterogeneous DG-face between elements 1 and 4 for this unit test.
   * 
   *   - lp: 3
   *   - fa: 1
   *   - face normal: (-1,0)
   *   - face area: 0.5
   *   - element volume: 4.0
   *
   *  el | sf |  xx              yy             xy               u        v
   *  0  | 0  |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *  0  | 1  |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *  0  | 2  |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *  1  | 0  |   2E9 7E9  2E9 | 6E9 2E9 6E9 | -100 -5E9 -100 |  6 2 6 | -30 -20 -30
   *  1  | 1  |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   *  1  | 2  |  1E10 9E9 1E10 | 4E9 5E9 4E9 | -2E9 -3E9 -2E9 |  5 4 5 | -20 -10 -20
   */
  // init DOFs
  double l_scDofs[2][5][3][3];

  for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
    for( unsigned short l_sf = 0; l_sf < 3; l_sf++ ) {
      l_scDofs[l_sd][0][l_sf][0] = 1E10;
      l_scDofs[l_sd][0][l_sf][1] =  9E9;
      l_scDofs[l_sd][0][l_sf][2] = 1E10;

      l_scDofs[l_sd][1][l_sf][0] =  4E9;
      l_scDofs[l_sd][1][l_sf][1] =  5E9;
      l_scDofs[l_sd][1][l_sf][2] =  4E9;

      l_scDofs[l_sd][2][l_sf][0] = -2E9;
      l_scDofs[l_sd][2][l_sf][1] = -3E9;
      l_scDofs[l_sd][2][l_sf][2] = -2E9;

      l_scDofs[l_sd][3][l_sf][0] =    5;
      l_scDofs[l_sd][3][l_sf][1] =    4;
      l_scDofs[l_sd][3][l_sf][2] =    5;

      l_scDofs[l_sd][4][l_sf][0] =  -20;
      l_scDofs[l_sd][4][l_sf][1] =  -10;
      l_scDofs[l_sd][4][l_sf][2] =  -20;
    }
  }

  // modify the right element, sub-face 0
  // modify sub-cell 3
  l_scDofs[1][0][0][0] =  2E9;
  l_scDofs[1][0][0][1] =  7E9;
  l_scDofs[1][0][0][2] =  2E9;

  l_scDofs[1][1][0][0] =  6E9;
  l_scDofs[1][1][0][1] =  2E9;
  l_scDofs[1][1][0][2] =  6E9;

  l_scDofs[1][2][0][0] = -100;
  l_scDofs[1][2][0][1] = -5E9;
  l_scDofs[1][2][0][2] = -100;

  l_scDofs[1][3][0][0] =    6;
  l_scDofs[1][3][0][1] =    2;
  l_scDofs[1][3][0][2] =    6;

  l_scDofs[1][4][0][0] =  -30;
  l_scDofs[1][4][0][1] =  -20;
  l_scDofs[1][4][0][2] =  -30;



  const double l_dMax = std::numeric_limits< double >::max();
  struct {
    double outNormal[3];
    double tangent0[3];
    double tangent1[3];
    double area;
  } l_charsFa[2] = { { {    1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.25 },
                     { {   -1.0,      0.0,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       { l_dMax,   l_dMax,  l_dMax},
                       0.5 } };
  
  struct {
    double volume;
  } l_charsEl[5] = { {l_dMax}, {4.0}, {2.0}, {3.0}, {5.0} };

  struct {
    double lam;
    double mu;
    double rho;
  } l_matPars[5] = { {l_dMax, l_dMax, l_dMax},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0},
                     {20.8E9, 10.4E9, 2600.0},
                     {32.4E9, 32.4E9, 2700.0} };

  const int l_iMax = std::numeric_limits< int >::max();
  int l_lpEl[4]  = {4, 2, 3, 1};
  int l_elFa[5][3] = { {l_iMax, l_iMax, l_iMax},
                       {0, 1, 0},
                       {0, 0, 1},
                       {1, 1, 0},
                       {0, 0, 1} };
  int l_elFaEl[5][3] = { {l_iMax, l_iMax, l_iMax},
                         {3, 4, 2},
                         {4, 1, 3},
                         {4, 1, 2},
                         {1, 3, 2} };

  // dummy data: triangles don't have additional sub-faces
  int l_elVe[1][3] = { {0,0,0} };
  struct {
    double coords[3];
  } l_charsVe[1] = { { {0,0,0} } };

  // wrapper for the LLF solvers
  edge::elastic::sc::Llf< double,
                          TRIA3,
                          2,
                          3 > l_llf;

  // dynamic memory allocations
  edge::data::Dynamic l_dynMem;

  // allocate memory for solvers of 4 limited plus elements
  l_llf.alloc( 4,
               l_dynMem );

  // init the solvers
  l_llf.init( 0,
              4,
              l_lpEl,
              l_elVe,
              l_elFa,
              l_elFaEl,
              l_charsVe,
              l_charsFa,
              l_charsEl,
              l_matPars );

  // derive reference solution
  double l_nuUt1[5][3][3];
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ )
    for( unsigned short l_sf = 0; l_sf < 3; l_sf++ )
      for( unsigned short l_cr = 0; l_cr < 3; l_cr++ )
        l_nuUt1[l_q1][l_sf][l_cr] = 0;


  // first jacobian
  double l_jac0[5][5] = { { 0,        0,  0,       -(20.8E9 + 2.0 * 10.4E9),  0      },
                          { 0,        0,  0,        -20.8E9,                  0      },
                          { 0,        0,  0,         0,                      -10.4E9 },
                          {-1.0/2600, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2600,  0,                       0      } };

  // second jacobian
  double l_jac1[5][5] = { { 0,        0,  0,       -(32.4E9 + 2.0 * 32.4E9),  0      },
                          { 0,        0,  0,        -32.4E9,                  0      },
                          { 0,        0,  0,         0,                      -32.4E9 },
                          {-1.0/2700, 0,  0,          0,                      0      },
                          { 0,        0, -1.0/2700,  0,                       0      } };

  // LLF contributions are averaged + left element's flux contribution
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      l_jac0[l_q1][l_q2] *= -0.5;
      l_jac1[l_q1][l_q2] *= -0.5;
    }
  }

  // viscous terms, +1: left element's viscosity cpntribution
  double l_vis1 = -0.5 * std::sqrt( (32.4E9 + 2 * 32.4E9)/2700 );

  // apply updates
  for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
    for( unsigned short l_sf = 0; l_sf < 3; l_sf++ ) {
      for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
        l_nuUt1[l_q1][l_sf][l_cr] += l_vis1 * l_scDofs[0][l_q1][l_sf][l_cr] * (0.5 / 3.0) * (9.0 / 4.0) * 0.006;
        l_nuUt1[l_q1][l_sf][l_cr] -= l_vis1 * l_scDofs[1][l_q1][l_sf][l_cr] * (0.5 / 3.0) * (9.0 / 4.0) * 0.006;
      }
    }


    for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
      for( unsigned short l_sf = 0; l_sf < 3; l_sf++ ) {
        for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
          l_nuUt1[l_q1][l_sf][l_cr] -=   l_scDofs[0][l_q2][l_sf][l_cr]
                                       * l_jac0[l_q1][l_q2] * (0.5 / 3.0) * (9.0 / 4.0) * 0.006;
          l_nuUt1[l_q1][l_sf][l_cr] -=   l_scDofs[1][l_q2][l_sf][l_cr]
                                       * l_jac1[l_q1][l_q2] * (0.5 / 3.0) * (9.0 / 4.0) * 0.006;
        }
      }
    }
  }

    // compute the result
    double l_res[5][3][3];

    l_llf.nuFaSf( 0.006,
                  3,
                  1,
                  l_scDofs,
                  l_res );

    // check the results
    for( unsigned short l_qt = 0; l_qt < 5; l_qt++ )
      for( unsigned short l_sf = 0; l_sf < 3; l_sf++ )
        for( unsigned short l_cr = 0; l_cr < 3; l_cr++ )
          REQUIRE( l_res[l_qt][l_sf][l_cr] == Approx( l_nuUt1[l_qt][l_sf][l_cr] ) );
}
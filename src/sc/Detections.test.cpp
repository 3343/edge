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
 * Unit test for the detection criteria.
 **/
#include <catch.hpp>
#define private public
#include "Detections.hpp"
#undef private


TEST_CASE( "Detection criteria: Admissibility through discrete maximum principle.", "[detections][dmp]" ) {
  /*
   * Our example:
   *
   * #quantities: 3
   * #runs: 2
   *
   * el: 4
   * elVeEl: 2, 1, 5, 9, 7
   *
   * current solution:
   *
   *                          min                                   max
   *  quantity     0           1           2             0           1           2
   *  run       0     1     0     1     0     1       0     1     0     1     0     1
   *  el      -0.1   0.3   1.0   2.0   2.0  -6.2     0.1   1.3   4.2   2.0   5.0  -6.2
   *  elVeEl0 -0.4   0.5   1.0   2.0  -5.1  -6.5     1.2   1.2   4.3   2.0   3.7   6.2
   *  elVeEl1 -0.1   0.2   3.1   2.0   2.2  -6.4     2.1   2.2   5.4   2.0   3.9   6.1
   *  elVeEl2 -0.6   1.2   1.0   2.0   3.3  -6.5     0.3   3.2   5.1   2.0   3.8   6.5
   *  elVeEl3 -0.2   0.2   2.0   2.0  -3.4  -6.3     0.1   4.2   4.2   2.0   5.0  -6.1
   *  elVeEl4  0.0   0.5   1.0   2.0   3.5  -6.2     0.3   2.2   4.3   2.0   4.0  -4.2
   *
   *  ext     -0.6   0.2   1.0   2.0  -5.1  -6.5     2.1   4.2   5.4   2.0   5.0   6.5
   *
   *  jumps
   *  2.7 4.0 4.4 0.0 10.1 13.0
   *
   *
   * candidate solution = current solution with a few modified values
   */
  double l_extCur1[10][2][3][2];
  double l_extCan1[10][2][3][2];

  // init
  for( unsigned short l_el = 0; l_el < 10; l_el++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCur1[l_el][0][l_qt][l_cr] = std::numeric_limits< double >::lowest();
        l_extCur1[l_el][1][l_qt][l_cr] = std::numeric_limits< double >::max();
        l_extCan1[l_el][0][l_qt][l_cr] = std::numeric_limits< double >::max();
        l_extCan1[l_el][1][l_qt][l_cr] = std::numeric_limits< double >::lowest();
      }
    }
  }

  // el
  l_extCur1[4][0][0][0] = -0.1;
  l_extCur1[4][0][0][1] =  0.3;

  l_extCur1[4][0][1][0] =  1.0;
  l_extCur1[4][0][1][1] =  2.0;

  l_extCur1[4][0][2][0] =  2.0;
  l_extCur1[4][0][2][1] = -6.2;

  l_extCur1[4][1][0][0] =  0.1;
  l_extCur1[4][1][0][1] =  1.3;

  l_extCur1[4][1][1][0] =  4.2;
  l_extCur1[4][1][1][1] =  2.0;

  l_extCur1[4][1][2][0] =  5.0;
  l_extCur1[4][1][2][1] = -6.2;

  // elVeEl0
  l_extCur1[2][0][0][0] = -0.4;
  l_extCur1[2][0][0][1] =  0.5;

  l_extCur1[2][0][1][0] =  1.0;
  l_extCur1[2][0][1][1] =  2.0;

  l_extCur1[2][0][2][0] = -5.1;
  l_extCur1[2][0][2][1] = -6.5;

  l_extCur1[2][1][0][0] =  1.2;
  l_extCur1[2][1][0][1] =  1.2;

  l_extCur1[2][1][1][0] =  4.3;
  l_extCur1[2][1][1][1] =  2.0;

  l_extCur1[2][1][2][0] =  3.7;
  l_extCur1[2][1][2][1] =  6.2;

  // elVeEl1
  l_extCur1[1][0][0][0] = -0.1;
  l_extCur1[1][0][0][1] =  0.2;

  l_extCur1[1][0][1][0] =  3.1;
  l_extCur1[1][0][1][1] =  2.0;

  l_extCur1[1][0][2][0] =  2.2;
  l_extCur1[1][0][2][1] = -6.4;

  l_extCur1[1][1][0][0] =  2.1;
  l_extCur1[1][1][0][1] =  2.2;

  l_extCur1[1][1][1][0] =  5.4;
  l_extCur1[1][1][1][1] =  2.0;

  l_extCur1[1][1][2][0] =  3.9;
  l_extCur1[1][1][2][1] =  6.1;

  // elVeEl2
  l_extCur1[5][0][0][0] = -0.6;
  l_extCur1[5][0][0][1] =  1.2;

  l_extCur1[5][0][1][0] =  1.0;
  l_extCur1[5][0][1][1] =  2.0;

  l_extCur1[5][0][2][0] =  3.3;
  l_extCur1[5][0][2][1] = -6.5;

  l_extCur1[5][1][0][0] =  0.3;
  l_extCur1[5][1][0][1] =  3.2;

  l_extCur1[5][1][1][0] =  5.1;
  l_extCur1[5][1][1][1] =  2.0;

  l_extCur1[5][1][2][0] =  3.8;
  l_extCur1[5][1][2][1] =  6.5;

  // elVeEl3
  l_extCur1[9][0][0][0] = -0.2;
  l_extCur1[9][0][0][1] =  0.2;

  l_extCur1[9][0][1][0] =  2.0;
  l_extCur1[9][0][1][1] =  2.0;

  l_extCur1[9][0][2][0] = -3.4;
  l_extCur1[9][0][2][1] = -6.3;

  l_extCur1[9][1][0][0] =  0.1;
  l_extCur1[9][1][0][1] =  4.2;

  l_extCur1[9][1][1][0] =  4.2;
  l_extCur1[9][1][1][1] =  2.0;

  l_extCur1[9][1][2][0] =  5.0;
  l_extCur1[9][1][2][1] = -6.1;

  // elVeEl4
  l_extCur1[7][0][0][0] =  0.0;
  l_extCur1[7][0][0][1] =  0.5;

  l_extCur1[7][0][1][0] =  1.0;
  l_extCur1[7][0][1][1] =  2.0;

  l_extCur1[7][0][2][0] =  3.5;
  l_extCur1[7][0][2][1] = -6.2;

  l_extCur1[7][1][0][0] =  0.3;
  l_extCur1[7][1][0][1] =  2.2;

  l_extCur1[7][1][1][0] =  4.3;
  l_extCur1[7][1][1][1] =  2.0;

  l_extCur1[7][1][2][0] =  4.0;
  l_extCur1[7][1][2][1] = -4.2;

  int  l_elId1 = 4;
  int  l_elVeEl1[5];
  l_elVeEl1[0] = 2;
  l_elVeEl1[1] = 1;
  l_elVeEl1[2] = 5;
  l_elVeEl1[3] = 9;
  l_elVeEl1[4] = 7;

  int  l_nElVeEl = 5;
  bool l_admiss1[2];

  /**
   * Test #1, current == candidate
   **/
  // copy over to candidate solution
  for( unsigned short l_ex = 0; l_ex < 1; l_ex++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCan1[4][l_ex][l_qt][l_cr] = l_extCur1[4][l_ex][l_qt][l_cr];
        l_extCan1[2][l_ex][l_qt][l_cr] = l_extCur1[2][l_ex][l_qt][l_cr];
        l_extCan1[1][l_ex][l_qt][l_cr] = l_extCur1[1][l_ex][l_qt][l_cr];
        l_extCan1[5][l_ex][l_qt][l_cr] = l_extCur1[5][l_ex][l_qt][l_cr];
        l_extCan1[9][l_ex][l_qt][l_cr] = l_extCur1[9][l_ex][l_qt][l_cr];
        l_extCan1[7][l_ex][l_qt][l_cr] = l_extCur1[7][l_ex][l_qt][l_cr];
      }
    }
  }

  edge::sc::Detections< 3, 2 >::dmp( l_extCur1[l_elId1],
                                     l_extCur1,
                                     l_extCan1[l_elId1],
                                     l_nElVeEl,
                                     l_elVeEl1,
                                     l_admiss1 );

  // check admissibility
  REQUIRE( l_admiss1[0] == true );
  REQUIRE( l_admiss1[1] == true );


  /**
   * Test #2
   **/
  // copy over to candidate solution
  for( unsigned short l_ex = 0; l_ex < 1; l_ex++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCan1[4][l_ex][l_qt][l_cr] = l_extCur1[4][l_ex][l_qt][l_cr];
      }
    }
  }

  // modify values
  l_extCan1[4][0][0][1] = 0.1961; // valid
  l_extCan1[4][1][2][0] = 5.0102; // invalid

  edge::sc::Detections< 3, 2 >::dmp( l_extCur1[l_elId1],
                                     l_extCur1,
                                     l_extCan1[l_elId1],
                                     l_nElVeEl,
                                     l_elVeEl1,
                                     l_admiss1 );

  // check admissibility
  REQUIRE( l_admiss1[0] == false );
  REQUIRE( l_admiss1[1] == true );

  /**
   * Test #3
   **/
  // copy over to candidate solution
  for( unsigned short l_ex = 0; l_ex < 1; l_ex++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCan1[4][l_ex][l_qt][l_cr] = l_extCur1[4][l_ex][l_qt][l_cr];
      }
    }
  }

  // modify values
  l_extCan1[4][0][0][1] = 0.1959; // invalid
  l_extCan1[4][1][2][0] = 5.010;  // valid

  edge::sc::Detections< 3, 2 >::dmp( l_extCur1[l_elId1],
                                     l_extCur1,
                                     l_extCan1[l_elId1],
                                     l_nElVeEl,
                                     l_elVeEl1,
                                     l_admiss1 );

  // check admissibility
  REQUIRE( l_admiss1[0] == true );
  REQUIRE( l_admiss1[1] == false );

  /**
   * Test #4
   **/
  // copy over to candidate solution
  for( unsigned short l_ex = 0; l_ex < 1; l_ex++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCan1[4][l_ex][l_qt][l_cr] = l_extCur1[4][l_ex][l_qt][l_cr];
      }
    }
  }

  // modify values
  l_extCan1[4][0][0][1] = 0.1959; // invalid
  l_extCan1[4][1][2][0] = 5.0102; // invalid

  edge::sc::Detections< 3, 2 >::dmp( l_extCur1[l_elId1],
                                     l_extCur1,
                                     l_extCan1[l_elId1],
                                     l_nElVeEl,
                                     l_elVeEl1,
                                     l_admiss1 );

  // check admissibility
  REQUIRE( l_admiss1[0] == false );
  REQUIRE( l_admiss1[1] == false );

  /**
   * Test #5
   **/
  // copy over to candidate solution
  for( unsigned short l_ex = 0; l_ex < 1; l_ex++ ) {
    for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        l_extCan1[4][l_ex][l_qt][l_cr] = l_extCur1[4][l_ex][l_qt][l_cr];
      }
    }
  }

  // modify values
  l_extCan1[4][0][0][1] = 0.1961; // valid
  l_extCan1[4][1][2][0] = 5.010;  // valid

  edge::sc::Detections< 3, 2 >::dmp( l_extCur1[l_elId1],
                                     l_extCur1,
                                     l_extCan1[l_elId1],
                                     l_nElVeEl,
                                     l_elVeEl1,
                                     l_admiss1 );

  // check admissibility
  REQUIRE( l_admiss1[0] == true );
  REQUIRE( l_admiss1[1] == true );
}

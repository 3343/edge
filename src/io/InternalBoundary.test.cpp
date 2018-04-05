/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Unit tests for internal boundary plotter.
 **/

#include <catch.hpp>

#define private public
#include "InternalBoundary.hpp"
#undef private

TEST_CASE( "Tests the assembly of the sub-face mesh for triangular elements.", "[io/InternalBoundary][meshTria3]" ) {
  /*
   * Our example setup:
   *
   *  1 *         *         *
   *    * *       0 *       * 2
   *    *   *     0   *     *   2
   *    *  4  *   0  1  *   *  2  2
   *    *       * 0       * *       2
   *  0 *33333333**11111111************
   *   -1         0        1          2 
   * 
   * #bf: 4
   * 
   * bf | be 0 | be 1 | fId 0 | fId 1
   * 0  | 0    | inf  | 2     | inf
   * 1  | 0    | inf  | 0     | inf
   * 2  | 1    | inf  | 1     | inf
   * 3  | 2    | inf  | 0     | inf
   *
   * be | de
   * 0  | 1
   * 1  | 2
   * 2  | 4
   *
   * el |   veCrds 0 | veCrds 1 | veCrds 2
   * 0  | 0:  inf inf  | 0: inf inf  | 4: inf inf 
   * 1  | 2:  0.0 0.0  | 1: 1.0 0.0  | 8: 0.0 1.0
   * 2  | 1:  1.0 0.0  | 3: 2.0 0.0  | 7: 1.0 1.0
   * 3  | 0:  inf inf  | 0: inf inf  | 0: inf inf
   * 4  | 5: -1.0 0.0  | 2: 0.0 0.0  | 6: -1.0 1.0
   */

  // dynamic memory
  edge::data::Dynamic l_dynMem;

  int l_bfBe[4][2] = { { 0, std::numeric_limits< int >::max() },
                       { 0, std::numeric_limits< int >::max() },
                       { 1, std::numeric_limits< int >::max() },
                       { 2, std::numeric_limits< int >::max() }
                     };

  int l_beEl[3] = { 1, 2, 4 };

  struct {
    double coords[3];
  } l_charsVe[9] = { { std::numeric_limits< double >::max(), std::numeric_limits< double >::max() },
                     {  1.0, 0.0 },
                     {  0.0, 0.0 },
                     {  2.0, 0.0 },
                     { std::numeric_limits< double >::max(), std::numeric_limits< double >::max() },
                     { -1.0, 0.0 },
                     { -1.0, 1.0 },
                     {  1.0, 1.0 },
                     {  0.0, 1.0 } };

  int l_elVe[5][3] = { {0, 0, 4},
                       {2, 1, 8},
                       {1, 3, 7},
                       {0, 0, 0},
                       {5, 2, 6} };

  struct  {
    double coords[2];
  } l_charsSv[10] = { {  0.0,                0.0                },
                      {  0.3333333333333333, 0.0                },
                      {  0.6666666666666666, 0.0                },
                      {  1.0,                0.0                },
                      {  0.0,                0.3333333333333333 },
                      {  0.3333333333333333, 0.3333333333333333 },
                      {  0.6666666666666666, 0.3333333333333333 },
                      {  0.0,                0.6666666666666666 },
                      {  0.3333333333333333, 0.6666666666666666 },
                      {  0.0,                1.0                } };

  struct {
    unsigned short fIdBfEl[2];
  } l_charsBf[4] = {
    {2, std::numeric_limits< unsigned short >::max() },
    {0, std::numeric_limits< unsigned short >::max() },
    {1, std::numeric_limits< unsigned short >::max() },
    {0, std::numeric_limits< unsigned short >::max() }
  };

  unsigned short l_scSv[18][3] = {
    {  1,5,4,  },
    {  2,6,5,  },
    {  5,8,7,  },
    {  0,1,4,  },
    {  1,2,5,  },
    {  2,3,6,  },
    {  5,6,8,  },
    {  7,8,9,  },
    {  4,5,7,  },
    {  0,1,std::numeric_limits< unsigned short >::max(),  },
    {  1,2,std::numeric_limits< unsigned short >::max(),  },
    {  2,3,std::numeric_limits< unsigned short >::max(),  },
    {  3,6,std::numeric_limits< unsigned short >::max(),  },
    {  6,8,std::numeric_limits< unsigned short >::max(),  },
    {  8,9,std::numeric_limits< unsigned short >::max(),  },
    {  9,7,std::numeric_limits< unsigned short >::max(),  },
    {  7,4,std::numeric_limits< unsigned short >::max(),  },
    {  4,0,std::numeric_limits< unsigned short >::max(),  }
  };

  // internal boundary writer
  std::string l_path = "/tmp/output/wpppp";
  edge::io::InternalBoundary< int, TRIA3, 2 > l_iBndWriter( l_path );

  // allocate memory
  l_iBndWriter.alloc( 4,
                      33,
                      l_dynMem );

  // initialize
  l_iBndWriter.init( 4,
                     l_scSv,
                     l_bfBe,
                     l_beEl,
                     l_elVe,
                     l_charsSv,
                     l_charsVe,
                     l_charsBf );

  // check internal data structures
  REQUIRE( l_iBndWriter.m_svCrds[0][0][0] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][1][0] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][2][0] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][3][0] == Approx(0.0) );

  REQUIRE( l_iBndWriter.m_svCrds[0][0][1] == Approx(1.0    ) );
  REQUIRE( l_iBndWriter.m_svCrds[0][1][1] == Approx(2.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][2][1] == Approx(1.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][3][1] == Approx(0.0    ) );

  REQUIRE( l_iBndWriter.m_svCrds[0][0][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][1][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][2][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[0][3][2] == Approx(0.0) );


  REQUIRE( l_iBndWriter.m_svCrds[1][0][0] == Approx(0.0)     );
  REQUIRE( l_iBndWriter.m_svCrds[1][1][0] == Approx(1.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][2][0] == Approx(2.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][3][0] == Approx(1.0)     );

  REQUIRE( l_iBndWriter.m_svCrds[1][0][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][1][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][2][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][3][1] == Approx(0.0) );

  REQUIRE( l_iBndWriter.m_svCrds[1][0][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][1][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][2][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[1][3][2] == Approx(0.0) );


  REQUIRE( l_iBndWriter.m_svCrds[2][0][0] == Approx(1.0 + 1.0)     );
  REQUIRE( l_iBndWriter.m_svCrds[2][1][0] == Approx(1.0 + 2.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][2][0] == Approx(1.0 + 1.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][3][0] == Approx(1.0 + 0.0)     );

  REQUIRE( l_iBndWriter.m_svCrds[2][0][1] == Approx(0.0)     );
  REQUIRE( l_iBndWriter.m_svCrds[2][1][1] == Approx(1.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][2][1] == Approx(2.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][3][1] == Approx(1.0)     );

  REQUIRE( l_iBndWriter.m_svCrds[2][0][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][1][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][2][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[2][3][2] == Approx(0.0) );


  REQUIRE( l_iBndWriter.m_svCrds[3][0][0] == Approx(-1.0 + 0.0)     );
  REQUIRE( l_iBndWriter.m_svCrds[3][1][0] == Approx(-1.0 + 1.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][2][0] == Approx(-1.0 + 2.0/3.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][3][0] == Approx(-1.0 + 1.0)     );

  REQUIRE( l_iBndWriter.m_svCrds[3][0][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][1][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][2][1] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][3][1] == Approx(0.0) );

  REQUIRE( l_iBndWriter.m_svCrds[3][0][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][1][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][2][2] == Approx(0.0) );
  REQUIRE( l_iBndWriter.m_svCrds[3][3][2] == Approx(0.0) );

  REQUIRE( l_iBndWriter.m_bfSfSv[0][0][0] == 0 );
  REQUIRE( l_iBndWriter.m_bfSfSv[0][0][1] == 1 );

  REQUIRE( l_iBndWriter.m_bfSfSv[0][1][0] == 1 );
  REQUIRE( l_iBndWriter.m_bfSfSv[0][1][1] == 2 );

  REQUIRE( l_iBndWriter.m_bfSfSv[0][2][0] == 2 );
  REQUIRE( l_iBndWriter.m_bfSfSv[0][2][1] == 3 );


  REQUIRE( l_iBndWriter.m_bfSfSv[1][0][0] == 4 );
  REQUIRE( l_iBndWriter.m_bfSfSv[1][0][1] == 5 );

  REQUIRE( l_iBndWriter.m_bfSfSv[1][1][0] == 5 );
  REQUIRE( l_iBndWriter.m_bfSfSv[1][1][1] == 6 );

  REQUIRE( l_iBndWriter.m_bfSfSv[1][2][0] == 6 );
  REQUIRE( l_iBndWriter.m_bfSfSv[1][2][1] == 7 );


  REQUIRE( l_iBndWriter.m_bfSfSv[2][0][0] == 8 );
  REQUIRE( l_iBndWriter.m_bfSfSv[2][0][1] == 9 );

  REQUIRE( l_iBndWriter.m_bfSfSv[2][1][0] == 9 );
  REQUIRE( l_iBndWriter.m_bfSfSv[2][1][1] == 10 );

  REQUIRE( l_iBndWriter.m_bfSfSv[2][2][0] == 10 );
  REQUIRE( l_iBndWriter.m_bfSfSv[2][2][1] == 11 );


  REQUIRE( l_iBndWriter.m_bfSfSv[3][0][0] == 12 );
  REQUIRE( l_iBndWriter.m_bfSfSv[3][0][1] == 13 );

  REQUIRE( l_iBndWriter.m_bfSfSv[3][1][0] == 13 );
  REQUIRE( l_iBndWriter.m_bfSfSv[3][1][1] == 14 );

  REQUIRE( l_iBndWriter.m_bfSfSv[3][2][0] == 14 );
  REQUIRE( l_iBndWriter.m_bfSfSv[3][2][1] == 15 );

  // check the buffer copy
  double l_data[4][3][7];
  for( unsigned short l_bf = 0; l_bf < 4; l_bf++ )
    for( unsigned short l_sf = 0; l_sf < 3; l_sf++ )
      for( unsigned short l_qt = 0; l_qt < 7; l_qt++ )
        l_data[l_bf][l_sf][l_qt] = l_bf*3*7 + l_sf*7 + l_qt;

  l_iBndWriter.copy( 1,              // 1 is first bf
                     3,              // 3 total bfs
                     4,              // plot 4 qts
                     7,              // stride is 7
                     l_data[0][0] ); // start at 0,0

  REQUIRE( l_iBndWriter.m_buffer[ 0] == Approx(1*3*7 + 0*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 1] == Approx(1*3*7 + 1*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 2] == Approx(1*3*7 + 2*7.0 + 0) );

  REQUIRE( l_iBndWriter.m_buffer[ 3] == Approx(2*3*7 + 0*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 4] == Approx(2*3*7 + 1*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 5] == Approx(2*3*7 + 2*7.0 + 0) );

  REQUIRE( l_iBndWriter.m_buffer[ 6] == Approx(3*3*7 + 0*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 7] == Approx(3*3*7 + 1*7.0 + 0) );
  REQUIRE( l_iBndWriter.m_buffer[ 8] == Approx(3*3*7 + 2*7.0 + 0) );


  REQUIRE( l_iBndWriter.m_buffer[ 9] == Approx(1*3*7 + 0*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[10] == Approx(1*3*7 + 1*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[11] == Approx(1*3*7 + 2*7.0 + 1) );

  REQUIRE( l_iBndWriter.m_buffer[12] == Approx(2*3*7 + 0*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[13] == Approx(2*3*7 + 1*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[14] == Approx(2*3*7 + 2*7.0 + 1) );

  REQUIRE( l_iBndWriter.m_buffer[15] == Approx(3*3*7 + 0*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[16] == Approx(3*3*7 + 1*7.0 + 1) );
  REQUIRE( l_iBndWriter.m_buffer[17] == Approx(3*3*7 + 2*7.0 + 1) );


  REQUIRE( l_iBndWriter.m_buffer[18] == Approx(1*3*7 + 0*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[19] == Approx(1*3*7 + 1*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[20] == Approx(1*3*7 + 2*7.0 + 2) );

  REQUIRE( l_iBndWriter.m_buffer[21] == Approx(2*3*7 + 0*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[22] == Approx(2*3*7 + 1*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[23] == Approx(2*3*7 + 2*7.0 + 2) );

  REQUIRE( l_iBndWriter.m_buffer[24] == Approx(3*3*7 + 0*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[25] == Approx(3*3*7 + 1*7.0 + 2) );
  REQUIRE( l_iBndWriter.m_buffer[26] == Approx(3*3*7 + 2*7.0 + 2) );


  REQUIRE( l_iBndWriter.m_buffer[27] == Approx(1*3*7 + 0*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[28] == Approx(1*3*7 + 1*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[29] == Approx(1*3*7 + 2*7.0 + 3) );

  REQUIRE( l_iBndWriter.m_buffer[30] == Approx(2*3*7 + 0*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[31] == Approx(2*3*7 + 1*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[32] == Approx(2*3*7 + 2*7.0 + 3) );

  REQUIRE( l_iBndWriter.m_buffer[33] == Approx(3*3*7 + 0*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[34] == Approx(3*3*7 + 1*7.0 + 3) );
  REQUIRE( l_iBndWriter.m_buffer[35] == Approx(3*3*7 + 2*7.0 + 3) );

  std::string m_varNames[4];
  char const * m_varNamesC[4];

  for( unsigned short l_qt = 0; l_qt < 4; l_qt++ ) {
    m_varNames[l_qt] =  std::to_string( l_qt );
    m_varNamesC[l_qt] = m_varNames[l_qt].c_str();
  }

  l_iBndWriter.write( 0,           // 1 is first bf
                      4,           // 3 total bfs
                      4,           // plot 4 qts
                      7,           // stride is 7
                      m_varNamesC,
                      l_data[0][0] ); // start at 0,0

  l_iBndWriter.write( 0,           // 1 is first bf
                      4,           // 3 total bfs
                      4,           // plot 4 qts
                      7,           // stride is 7
                      m_varNamesC,
                      l_data[0][0] ); // start at 0,0

}
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
 * Unit test for the data layout.
 **/
#include <catch.hpp>
#define private public
#include "DataLayout.hpp"
#undef private


TEST_CASE( "Data Layout: Sparse adjacency.", "[dataLayout][spAd]" ) {
  /*
   * Our example
   *
   *   en | tg | inner | send tg | recv tg | rank | enEn
   *   0  | 0  | x     |         |         |      | 1,x,5
   *   1  | 0  | x     |         |         |      | 2,0,3
   *   2  | 0  | x     |         |         |      | x,0,6
   *   3  | 0  |       | 0       |         | 0    | 2,x,8
   *   4  | 0  |       | 0       |         | 0    | 7,x,5
   *   5  | 0  |       | 1       |         | 0    | 8,4,2
   *   6  | 0  |       | 0       |         | 1    | 0,x,x
   *   7  | 0  |       |         | 0       | 0    | x,3,4
   *   8  | 0  |       |         | 0       | 0    | 5,4,x
   *   9  | 0  |       |         | 0       | 1    | 5,4,6
   */
  t_enLayout l_tg1;
  l_tg1.nEnts = 10;
  l_tg1.timeGroups.resize(1);
  l_tg1.timeGroups[0].nEntsOwn = 7;
  l_tg1.timeGroups[0].nEntsNotOwn = 3;
  l_tg1.timeGroups[0].inner.first = 0;
  l_tg1.timeGroups[0].inner.size  = 3;
  l_tg1.timeGroups[0].send.resize(3);
  l_tg1.timeGroups[0].receive.resize(3);
  l_tg1.timeGroups[0].neRanks.resize(3);
  l_tg1.timeGroups[0].neTgs.resize(3);

  l_tg1.timeGroups[0].send[0].first = 3;
  l_tg1.timeGroups[0].send[0].size  = 2;
  l_tg1.timeGroups[0].send[1].first = 5;
  l_tg1.timeGroups[0].send[1].size  = 1;
  l_tg1.timeGroups[0].send[2].first = 6;
  l_tg1.timeGroups[0].send[2].size  = 1;

  l_tg1.timeGroups[0].receive[0].first = 7;
  l_tg1.timeGroups[0].receive[0].size  = 2;
  l_tg1.timeGroups[0].receive[1].first = 9;
  l_tg1.timeGroups[0].receive[1].size  = 0;
  l_tg1.timeGroups[0].receive[2].first = 9;
  l_tg1.timeGroups[0].receive[2].size  = 1;

  double  l_dataRaw1[3*10][87];
  double (*l_dataPtrs1[10][3])[87];
  int l_lpFaLp1[10][3] = {
    {1, std::numeric_limits< int >::max(), 5},
    {2,0,3},
    {std::numeric_limits< int >::max(),0,6},
    {2,std::numeric_limits< int >::max(),8},
    {8,std::numeric_limits< int >::max(),5},
    {8,4,2},
    {0,std::numeric_limits< int >::max(),std::numeric_limits< int >::max()},
    {std::numeric_limits< int >::max(),3,4},
    {5,4,std::numeric_limits< int >::max()},
    {5,4,6}
  };

  std::vector< std::vector< unsigned char * > > l_send1;
  std::vector< std::vector< unsigned char * > > l_recv1;

  edge::data::DataLayout::adj( 3,
                               l_lpFaLp1[0],
                               l_tg1,
                               l_dataRaw1,
                               l_dataPtrs1,
                               l_send1,
                               l_recv1 );

  // check the pointers
  REQUIRE( l_dataPtrs1[0][0] == l_dataRaw1 +  0 );
  REQUIRE( l_dataPtrs1[0][1] == nullptr         );
  REQUIRE( l_dataPtrs1[0][2] == l_dataRaw1 +  1 );

  REQUIRE( l_dataPtrs1[1][0] == l_dataRaw1 +  2 );
  REQUIRE( l_dataPtrs1[1][1] == l_dataRaw1 +  3 );
  REQUIRE( l_dataPtrs1[1][2] == l_dataRaw1 +  4 );

  REQUIRE( l_dataPtrs1[2][0] == nullptr         );
  REQUIRE( l_dataPtrs1[2][1] == l_dataRaw1 +  5 );
  REQUIRE( l_dataPtrs1[2][2] == l_dataRaw1 +  6 );

  // send
  REQUIRE( l_dataPtrs1[3][0] == l_dataRaw1 +  7 );
  REQUIRE( l_dataPtrs1[3][1] == nullptr         );
  REQUIRE( l_dataPtrs1[3][2] == l_dataRaw1 + 12 );

  REQUIRE( l_dataPtrs1[4][0] == l_dataRaw1 + 13 );
  REQUIRE( l_dataPtrs1[4][1] == nullptr         );
  REQUIRE( l_dataPtrs1[4][2] == l_dataRaw1 +  8 );

  REQUIRE( l_dataPtrs1[5][0] == nullptr         );
  REQUIRE( l_dataPtrs1[5][1] == l_dataRaw1 +  9 );
  REQUIRE( l_dataPtrs1[5][2] == l_dataRaw1 + 10 );

  REQUIRE( l_dataPtrs1[6][0] == l_dataRaw1 + 11 );
  REQUIRE( l_dataPtrs1[6][1] == nullptr         );
  REQUIRE( l_dataPtrs1[6][2] == nullptr         );

  REQUIRE( l_dataPtrs1[7][0] == nullptr         );
  REQUIRE( l_dataPtrs1[7][1] == l_dataRaw1 + 14 );
  REQUIRE( l_dataPtrs1[7][2] == l_dataRaw1 + 15 );

  REQUIRE( l_dataPtrs1[8][0] == nullptr         );
  REQUIRE( l_dataPtrs1[8][1] == l_dataRaw1 + 16 );
  REQUIRE( l_dataPtrs1[8][2] == nullptr         );

  REQUIRE( l_dataPtrs1[9][0] == nullptr         );
  REQUIRE( l_dataPtrs1[9][1] == nullptr         );
  REQUIRE( l_dataPtrs1[9][2] == l_dataRaw1 + 17 );

  REQUIRE( l_send1.size() == 1 );
  REQUIRE( l_recv1.size() == 1 );
  REQUIRE( l_send1[0].size() == 3+1 );
  REQUIRE( l_recv1[0].size() == 3+1 );

  REQUIRE( l_send1[0][0] == (unsigned char*) (l_dataRaw1 + 12) );
  REQUIRE( l_send1[0][1] == (unsigned char*) (l_dataRaw1 + 14) );
  REQUIRE( l_send1[0][2] == (unsigned char*) (l_dataRaw1 + 14) );
  REQUIRE( l_send1[0][3] == (unsigned char*) (l_dataRaw1 + 14) );

  REQUIRE( l_recv1[0][0] == (unsigned char*) (l_dataRaw1 + 14) );
  REQUIRE( l_recv1[0][1] == (unsigned char*) (l_dataRaw1 + 17) );
  REQUIRE( l_recv1[0][2] == (unsigned char*) (l_dataRaw1 + 17) );
  REQUIRE( l_recv1[0][3] == (unsigned char*) (l_dataRaw1 + 18) );
}

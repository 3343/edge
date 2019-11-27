/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Unit tests for sparse entities.
 **/

#include <catch.hpp>
#include "SparseEntities.hpp"

TEST_CASE( "SparseEnts: Number of sparse entities.", "[nSp][SparseEnts]" ) {
  typedef struct { unsigned short spType; } t_enChars;

  t_enChars l_chars1[10];
  l_chars1[0].spType = 0;
  l_chars1[1].spType = 3; // 1
  l_chars1[2].spType = 2;
  l_chars1[3].spType = 4;
  l_chars1[4].spType = 3; // 2
  l_chars1[5].spType = 0;
  l_chars1[6].spType = 8;
  l_chars1[7].spType = 2;
  l_chars1[8].spType = 2;
  l_chars1[9].spType = 3; // 3

  REQUIRE( edge::data::SparseEntities::nSp( 10, 3, l_chars1 ) == 3 );
}

TEST_CASE( "SparseEnts: Extract a sparse from a dense entity layout.", "[denseToSparse][SparseEnts]" ) {
  // sparse type
  int_spType l_spType;
  // face char
  std::vector< t_faceChars > l_faChars;
  // dense and sparse layouts
  t_enLayout l_deLayout;
  t_enLayout l_spLayout;

  /*
   * Test 1
   *
   * entity "type": faces
   * #faces: 38
   * #time regions: 2
   *
   * Dense Layout
   *                 5/1                     8/2                               7/2
   * inner: [ 0, 1, 2, 1234, 3 || 4, 1234, 5, 1234, 6, 7, 8, 9 || 1234, 10, 11, 12, 13, 14, 1234 ] 
   *        <------ tg0 ------><------------ tg1 --------------><------------ tg2 --------------->
   *
   *                               9/4                             4/2                  5/1
   * send: [ 15, 16, 1234, 1234, 17, 1234, 18, 19, 1234 || 20, 1234, 21, 1234 || 22, 1234, 23, 24, 25 ]
   *       <------------------- tg0 --------------------><------- tg1 --------><--------- tg2 -------->
   */
  l_spType = 1234;

  l_faChars.resize(38);
  l_faChars[ 0].spType = 0;
  l_faChars[ 1].spType = 1;
  l_faChars[ 2].spType = 2;
  l_faChars[ 3].spType = 1234;
  l_faChars[ 4].spType = 3;
  l_faChars[ 5].spType = 4;
  l_faChars[ 6].spType = 1234;
  l_faChars[ 7].spType = 5;
  l_faChars[ 8].spType = 1234;
  l_faChars[ 9].spType = 6;
  l_faChars[10].spType = 7;
  l_faChars[11].spType = 8;
  l_faChars[12].spType = 9;
  l_faChars[13].spType = 1234;
  l_faChars[14].spType = 10;
  l_faChars[15].spType = 11;
  l_faChars[16].spType = 12;
  l_faChars[17].spType = 13;
  l_faChars[18].spType = 14;
  l_faChars[19].spType = 1234;
  l_faChars[20].spType = 15;
  l_faChars[21].spType = 16;
  l_faChars[22].spType = 1234;
  l_faChars[23].spType = 1234;
  l_faChars[24].spType = 17;
  l_faChars[25].spType = 1234;
  l_faChars[26].spType = 18;
  l_faChars[27].spType = 19;
  l_faChars[28].spType = 1234;
  l_faChars[29].spType = 20;
  l_faChars[30].spType = 1234;
  l_faChars[31].spType = 21;
  l_faChars[32].spType = 1234;
  l_faChars[33].spType = 22;
  l_faChars[34].spType = 1234;
  l_faChars[35].spType = 23;
  l_faChars[36].spType = 24;
  l_faChars[37].spType = 25;

  std::size_t l_nElsInDe[3] = {5, 8, 7};
  std::size_t l_nElsSeDe[3] = {9, 4, 5};
  std::size_t l_nElsInSp[3] = {0, 0, 0};
  std::size_t l_nElsSeSp[3] = {0, 0, 0};

  // derive sparse layout
  edge::data::SparseEntities::denseToSparse(  3,
                                              l_spType,
                                             &l_faChars[0],
                                              l_nElsInDe,
                                              l_nElsSeDe,
                                              l_nElsInSp,
                                              l_nElsSeSp );

  // check the results
  REQUIRE( l_nElsInSp[0] == 1 );
  REQUIRE( l_nElsInSp[1] == 2 );
  REQUIRE( l_nElsInSp[2] == 2 );

  REQUIRE( l_nElsSeSp[0] == 4 );
  REQUIRE( l_nElsSeSp[1] == 2 );
  REQUIRE( l_nElsSeSp[2] == 1 );
}

// TODO: Adjust unit test to MPI-parallel implementation

#if 0
TEST_CASE( "SparseEnts: Extract a sparse from a dense entity layout based on adjacent entities.", "[denseToSparseAdj][SparseEnts]" ) {
  // sparse type
  int_spType l_spType;
  // face char
  std::vector< t_faceChars > l_faChars;
  // dense and sparse layouts
  t_enLayout l_deLayout;
  t_enLayout l_spLayout;

  /*
   * Test 1
   *
   * #entities: 38
   * #time regions: 2
   * each dense entities is adjacent to two other entieis
   *
   *
   * Dense Layout
   *                                        13-0          13-1           2-2             13-0          13-1      2-2      mpi-tg
   * tg 0: [ 0-1, 0-3, 1-3, 3-4, 0-1 || 1-0, 2-4, 0-1 | 4-0, 3-1 | 0-0, 1-1, 0-1 || 2-4, 0-0, 0-1 | 1-0, 0-0 | 3-0, 4-0 ] adjacent entities
   *                         x                x          x                           x                               x    sparse type (1234 of adjacent) matchin?
   *       <------ inner ------------><------------- send -----------------------><-------------- recv ----------------->
   *
   *                                                           4-3     5-1   5-2       4-3       5-1     5-2
   * tg 1: [ 0-3, 3-0, 1-4, 4-1, 0-0, 4-3, 1-0, 3-0, 2-4 || 0-1, 0-4 | 1-1 | 2-1 || 0-3, 4-0 | 0-0 3-0 | 1-3 ]
   *                    x    x                        x           x           x           x
   *       <------------------- inner ------------------><------ send -----------><----------- recv --------->
   *
   *
   * Dense Layout of the adjacent entities:
   *
   * [ 1111, 1111,  1234, 1111, 1234 ]
   *    0     1       2     3     4
   */
  l_spType = 1234;

  // setup the dense layout
  l_deLayout.nEnts = 38;
  l_deLayout.timeGroups.resize( 2 );

  l_deLayout.timeGroups[0].nEntsOwn    = 13;
  l_deLayout.timeGroups[0].nEntsNotOwn =  7;
  l_deLayout.timeGroups[0].inner.first =  0;
  l_deLayout.timeGroups[0].inner.size  =  5;

  l_deLayout.timeGroups[0].neRanks.resize( 3 );
  l_deLayout.timeGroups[0].neRanks[0] = 13;
  l_deLayout.timeGroups[0].neRanks[1] = 13;
  l_deLayout.timeGroups[0].neRanks[2] = 2;

  l_deLayout.timeGroups[0].neTgs.resize( 3 );
  l_deLayout.timeGroups[0].neTgs[0] = 0;
  l_deLayout.timeGroups[0].neTgs[1] = 1;
  l_deLayout.timeGroups[0].neTgs[2] = 2;

  l_deLayout.timeGroups[0].send.resize(    3 );
  l_deLayout.timeGroups[0].receive.resize( 3 );

  l_deLayout.timeGroups[0].send[0].first    = 5;
  l_deLayout.timeGroups[0].send[0].size     = 3;

  l_deLayout.timeGroups[0].send[1].first    = 8;
  l_deLayout.timeGroups[0].send[1].size     = 2;

  l_deLayout.timeGroups[0].send[2].first    = 10;
  l_deLayout.timeGroups[0].send[2].size     = 3;

  l_deLayout.timeGroups[0].receive[0].first = 13;
  l_deLayout.timeGroups[0].receive[0].size  = 3;

  l_deLayout.timeGroups[0].receive[1].first = 16;
  l_deLayout.timeGroups[0].receive[1].size  = 2;

  l_deLayout.timeGroups[0].receive[2].first = 18;
  l_deLayout.timeGroups[0].receive[2].size  = 2;

  l_deLayout.timeGroups[1].nEntsOwn    = 13;
  l_deLayout.timeGroups[1].nEntsNotOwn = 5;
  l_deLayout.timeGroups[1].inner.first = 20;
  l_deLayout.timeGroups[1].inner.size  = 9;

  l_deLayout.timeGroups[1].neRanks.resize( 3 );
  l_deLayout.timeGroups[1].neRanks[0] = 4;
  l_deLayout.timeGroups[1].neRanks[1] = 5;
  l_deLayout.timeGroups[1].neRanks[2] = 5;

  l_deLayout.timeGroups[1].neTgs.resize( 3 );
  l_deLayout.timeGroups[1].neTgs[0] = 3;
  l_deLayout.timeGroups[1].neTgs[1] = 1;
  l_deLayout.timeGroups[1].neTgs[2] = 2;

  l_deLayout.timeGroups[1].send.resize(    3 );
  l_deLayout.timeGroups[1].receive.resize( 3 );

  l_deLayout.timeGroups[1].send[0].first    = 29;
  l_deLayout.timeGroups[1].send[0].size     = 2;

  l_deLayout.timeGroups[1].send[1].first    = 31;
  l_deLayout.timeGroups[1].send[1].size     = 1;

  l_deLayout.timeGroups[1].send[2].first    = 32;
  l_deLayout.timeGroups[1].send[2].size     = 1;

  l_deLayout.timeGroups[1].receive[0].first = 33;
  l_deLayout.timeGroups[1].receive[0].size  = 2;

  l_deLayout.timeGroups[1].receive[1].first = 35;
  l_deLayout.timeGroups[1].receive[1].size  = 2;

  l_deLayout.timeGroups[1].receive[2].first = 37;
  l_deLayout.timeGroups[1].receive[2].size  = 1;

  /*
   * Test 1
   *
   * #entities: 38
   * #time regions: 2
   * each dense entities is adjacent to two other entieis
   *
   *
   * Dense Layout
   *                                        13-0          13-1           2-2             13-0          13-1      2-2      mpi-tg
   * tg 0: [ 0-1, 0-3, 1-3, 3-4, 0-1 || 1-0, 2-4, 0-1 | 4-0, 3-1 | 0-0, 1-1, 0-1 || 2-4, 0-0, 0-1 | 1-0, 0-0 | 3-0, 4-0 ] adjacent entities
   *                         x                x          x                           x                               x    sparse type (1234 of adjacent) matchin?
   *       <------ inner ------------><------------- send -----------------------><-------------- recv ----------------->
   *
   *                                                           4-3     5-1   5-2       4-3       5-1     5-2
   * tg 1: [ 0-3, 3-0, 1-4, 4-1, 0-0, 3-3, 1-0, 3-0, 2-4 || 0-1, 0-4 | 1-1 | 2-1 || 0-3, 4-0 | 0-0 3-0 | 1-3 ]
   *                    x    x                        x           x           x           x
   *       <------------------- inner ------------------><------ send -----------><----------- recv --------->
   *
   *
   * Dense Layout of the adjacent entities:
   *
   * [ 1111, 1111,  1234, 1111, 1234 ]
   *    0     1       2     3     4
   */

  // setup the characteristics of the adjacent elements
  l_faChars.resize(5);
  l_faChars[0].spType = 1111;
  l_faChars[1].spType = 1111;
  l_faChars[2].spType = 1234;
  l_faChars[3].spType = 1111;
  l_faChars[4].spType = 1234;

  // setup the adjacency information
  int l_enToEn[38][2] = { {0, 1}, {0, 3}, {1, 3}, {3, 4}, {0, 1},
                          {1, 0}, {2, 4}, {0, 1}, {4, 0}, {3, 1}, {0, 0 }, {1, 1}, {0, 1},
                          {2, 4}, {0, 0}, {0, 1}, {1, 0}, {0, 0}, {3, 0 }, {4, 0},

                          {0, 3}, {3, 0}, {1, 4}, {4, 1}, {0, 0}, {4, 3}, {1, 0}, {3, 0}, {2, 4},
                          {0, 1}, {0, 4}, {1, 1}, {2, 1},
                          {0, 3}, {4, 0}, {0, 0}, {3, 0}, {1, 3} };

  // derive sparse layout
  edge::data::SparseEntities::denseToSparseAdj(  2,
                                                 l_enToEn[0],
                                                 l_spType,
                                                &l_faChars[0],
                                                 l_deLayout,
                                                 l_spLayout );

  // check the results
  REQUIRE( l_spLayout.timeGroups[0].nEntsOwn    == 3 );
  REQUIRE( l_spLayout.timeGroups[0].nEntsNotOwn == 2 );

  REQUIRE( l_spLayout.timeGroups[0].inner.first == 0 );
  REQUIRE( l_spLayout.timeGroups[0].inner.size  == 1 );

  REQUIRE( l_spLayout.timeGroups[0].neRanks[0] == 13 );
  REQUIRE( l_spLayout.timeGroups[0].neRanks[1] == 13 );
  REQUIRE( l_spLayout.timeGroups[0].neRanks[2] == 2  );

  REQUIRE( l_spLayout.timeGroups[0].neTgs[0] == 0 );
  REQUIRE( l_spLayout.timeGroups[0].neTgs[1] == 1 );
  REQUIRE( l_spLayout.timeGroups[0].neTgs[2] == 2 );

  REQUIRE( l_spLayout.timeGroups[0].send[0].first == 1 );
  REQUIRE( l_spLayout.timeGroups[0].send[0].size == 1 );

  REQUIRE( l_spLayout.timeGroups[0].send[1].first == 2 );
  REQUIRE( l_spLayout.timeGroups[0].send[1].size == 1 );

  REQUIRE( l_spLayout.timeGroups[0].send[2].first == 3 );
  REQUIRE( l_spLayout.timeGroups[0].send[2].size == 0 );

  REQUIRE( l_spLayout.timeGroups[0].receive[0].first == 3 );
  REQUIRE( l_spLayout.timeGroups[0].receive[0].size == 1 );

  REQUIRE( l_spLayout.timeGroups[0].receive[1].first == 4 );
  REQUIRE( l_spLayout.timeGroups[0].receive[1].size == 0 );

  REQUIRE( l_spLayout.timeGroups[0].receive[2].first == 4 );
  REQUIRE( l_spLayout.timeGroups[0].receive[2].size == 1 );


  REQUIRE( l_spLayout.timeGroups[1].nEntsOwn    == 6 );
  REQUIRE( l_spLayout.timeGroups[1].nEntsNotOwn == 1 );

  REQUIRE( l_spLayout.timeGroups[1].inner.first == 5 );
  REQUIRE( l_spLayout.timeGroups[1].inner.size  == 4 );

  REQUIRE( l_spLayout.timeGroups[1].neRanks[0] == 4 );
  REQUIRE( l_spLayout.timeGroups[1].neRanks[1] == 5 );
  REQUIRE( l_spLayout.timeGroups[1].neRanks[2] == 5 );

  REQUIRE( l_spLayout.timeGroups[1].neTgs[0] == 3 );
  REQUIRE( l_spLayout.timeGroups[1].neTgs[1] == 1 );
  REQUIRE( l_spLayout.timeGroups[1].neTgs[2] == 2 );

  REQUIRE( l_spLayout.timeGroups[1].send[0].first == 9 );
  REQUIRE( l_spLayout.timeGroups[1].send[0].size == 1 );

  REQUIRE( l_spLayout.timeGroups[1].send[1].first == 10 );
  REQUIRE( l_spLayout.timeGroups[1].send[1].size == 0 );

  REQUIRE( l_spLayout.timeGroups[1].send[2].first == 10 );
  REQUIRE( l_spLayout.timeGroups[1].send[2].size == 1 );

  REQUIRE( l_spLayout.timeGroups[1].receive[0].first == 11 );
  REQUIRE( l_spLayout.timeGroups[1].receive[0].size == 1 );

  REQUIRE( l_spLayout.timeGroups[1].receive[1].first == 12 );
  REQUIRE( l_spLayout.timeGroups[1].receive[1].size == 0 );

  REQUIRE( l_spLayout.timeGroups[1].receive[2].first == 12 );
  REQUIRE( l_spLayout.timeGroups[1].receive[2].size == 0 );
}
#endif

TEST_CASE( "SparseEnts: First sparse ids and sizes of subregions based on their entities.", "[subRgnsSpIdEn][SparseEnts]" ) {
  /**
   * Our setup:
   *
   *   dense ids:  [  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19 ]
   *   sparse ids: [        0     1  2        3          4       5   6           7       8 ]
   *
   *   region ex1:                     [3          |4      |5      |6              ]
   *                                        sr0        sr1    sr2           sr3
   *
   *   region ex2:        [0                 |3                        |7              ]
   *                               sr0                   sr1                  sr2
   *
   * Input example 1:
   *   region first: 6
   *   sub-region #0 size: 4
   *   sub-region #1 size: 2
   *   sub-region #2 size: 2
   *   sub-region #3 size: 4
   *
   * Output example 1:
   *   first: 3 4 5 6
   *   sizes: 1 1 1 2
   *
   * Input example 2:
   *   region first: 2
   *   sub-region #0 size: 6
   *   sub-region #1 size: 7
   *   sub-region #2 size: 4
   *
   * Output example 2:
   *   first: 0 3 7
   *   sizes: 3 4 1
   **/

  typedef struct { int spType; } Chars;
  Chars l_deChars[20];
  l_deChars[ 0].spType = 15;
  l_deChars[ 1].spType = 25;
  l_deChars[ 2].spType = 99;
  l_deChars[ 3].spType = 55;
  l_deChars[ 4].spType = 99;
  l_deChars[ 5].spType = 99;
  l_deChars[ 6].spType = 15;
  l_deChars[ 7].spType = 15;
  l_deChars[ 8].spType = 99;
  l_deChars[ 9].spType = 15;
  l_deChars[10].spType = 15;
  l_deChars[11].spType = 99;
  l_deChars[12].spType = 10;
  l_deChars[13].spType = 99;
  l_deChars[14].spType = 99;
  l_deChars[15].spType = 12;
  l_deChars[16].spType = 25;
  l_deChars[17].spType = 99;
  l_deChars[18].spType = 55;
  l_deChars[19].spType = 99;

  // example 1
  int l_subRgnSizes1[4] = { 4, 2, 2, 4 };
  int l_first1[4];
  int l_sizes1[4];

  // call our derivation
  edge::data::SparseEntities::subRgnsSpId( 6,
                                           4,
                                           l_subRgnSizes1,
                                           99,
                                           l_deChars,
                                           l_first1,
                                           l_sizes1 );

  REQUIRE( l_first1[0] == 3 );
  REQUIRE( l_first1[1] == 4 );
  REQUIRE( l_first1[2] == 5 );
  REQUIRE( l_first1[3] == 6 );

  REQUIRE( l_sizes1[0] == 1 );
  REQUIRE( l_sizes1[1] == 1 );
  REQUIRE( l_sizes1[2] == 1 );
  REQUIRE( l_sizes1[3] == 2 );


  // example 2
  int l_subRgnSizes2[4] = { 6, 7, 4 };
  int l_first2[3];
  int l_sizes2[3];

  // call our derivation
  edge::data::SparseEntities::subRgnsSpId( 2,
                                           3,
                                           l_subRgnSizes2,
                                           99,
                                           l_deChars,
                                           l_first2,
                                           l_sizes2 );

  REQUIRE( l_first2[0] == 0 );
  REQUIRE( l_first2[1] == 3 );
  REQUIRE( l_first2[2] == 7 );

  REQUIRE( l_sizes2[0] == 3 );
  REQUIRE( l_sizes2[1] == 4 );
  REQUIRE( l_sizes2[2] == 1 );
}

TEST_CASE( "SparseEnts: Links sparse entities to dense entities.", "[linkSpDe][SparseEnts]" ) {
  /*
   * Setup:
   *
   *    SpId | DeId         SpId | DeId
   *         | 0               0 | 2
   *         | 1               1 | 4
   *       0 | 2               2 | 5
   *         | 3               3 | 7
   *       1 | 4
   *       2 | 5
   *         | 6
   *       3 | 7
   *         | 8
   */
  typedef struct{ unsigned short spType; } Chars;

  Chars l_chars1[9] = { {8}, {4}, {33}, {0}, {33}, {33}, {17}, {33}, {18} };
  int l_link1[4];

  edge::data::SparseEntities::linkSpDe( 9,
                                        33,
                                        l_chars1,
                                        l_link1 );

  // check the results
  REQUIRE( l_link1[0] == 2 );
  REQUIRE( l_link1[1] == 4 );
  REQUIRE( l_link1[2] == 5 );
  REQUIRE( l_link1[3] == 7 );
}

TEST_CASE( "SparseEnts: Links sparse entities to sparse entities.", "[linkSpSp][SparseEnts]" ) {
  /*
   * Setup:
     *
     *   Input (implicit):    Output
     *
     *    SpId0 | SpId1 | DeId         SpId0 | SpId1
     *          | 0     | 0                0 | 1
     *          |       | 1                1 | -
     *        0 | 1     | 2                2 | -
     *          |       | 3                3 | 2
     *        1 |       | 4
     *        2 |       | 5
     *          |       | 6
     *        3 | 2     | 7
     *          |       | 8
   */
  typedef struct{ unsigned short spType; } Chars;

  unsigned short const l_spType0 = 16;
  unsigned short const l_spType1 = 32;
  unsigned short const l_spTypeB = 16 + 32;

  Chars l_chars1[9] = { {8+l_spType1}, {4}, {4+l_spTypeB}, {0}, {l_spType0}, {l_spType0}, {12}, {l_spTypeB}, {13} };
  int l_link1[4];

  edge::data::SparseEntities::linkSpSp( 9,
                                        l_spType0,
                                        l_spType1,
                                        l_chars1,
                                        l_link1 );

  // check the results
  REQUIRE( l_link1[0] == 1 );
  REQUIRE( l_link1[1] == std::numeric_limits< int >::max() );
  REQUIRE( l_link1[2] == std::numeric_limits< int >::max() );
  REQUIRE( l_link1[3] == 2 );
}

TEST_CASE( "SparseEnts: Links sparse entities based on adjacency information (single sparse type).", "[linkSpAdjSst][SparseEnts]" ) {
  /**
   * Setup:
   *
   *   Dense adjacency      faces with        dense illustration  compressed sparse info
   *   information of       sparse ids        of derived          as given as output
   *   elements connected   (implicit info):  adjacency info:     of this function:
   *   to faces:
   *
   *   el | fa              fa | spId         el | spId           elSpId | faSpId
   *   0  | 1-4             0  | -            0  | x-x            0      | 2-1
   *   1  | 5-3             1  | -            1  | 2-1            1      | 3-x
   *   2  | 9-7             2  | 0            2  | 3-x            2      | x-1
   *   3  | 8-3             3  | 1            3  | x-1            3      | 0-x
   *   4  | 2-1             4  | -            4  | 0-x            4      | 1-x
   *   5  | 4-6             5  | 2            5  | x-x            5      | 3-0
   *   6  | 3-8             6  | -            6  | 1-x
   *   7  | 0-4             7  | -            7  | x-x
   *   8  | 9-2             8  | -            8  | 3-0
   *   9  | 4-1             9  | 3            9  | x-x
   *   10 | 4-1                               10 | x-x
   **/

  // set input info info
  int l_nEn = 11;
  int l_nAdjPerEn = 2;

  int l_enToEn[11][2] = { {1,4}, {5,3}, {9,7}, {8,3}, {2,1},
                          {4,6}, {3,8}, {0,4}, {9,2}, {4,1},
                          {4,1} };

  unsigned short l_spType = 87;

  typedef struct{ unsigned short spType; } CharsAdj;
  CharsAdj l_charsAdj[10];
  l_charsAdj[0].spType = 5;
  l_charsAdj[1].spType = 8;
  l_charsAdj[2].spType = 87;
  l_charsAdj[3].spType = 87;
  l_charsAdj[4].spType = 5;
  l_charsAdj[5].spType = 87;
  l_charsAdj[6].spType = 0;
  l_charsAdj[7].spType = 5;
  l_charsAdj[8].spType = 3;
  l_charsAdj[9].spType = 87;

  int l_spLink[6][2];

  // link the sparse entities
  edge::data::SparseEntities::linkSpAdjSst( l_nEn,
                                            l_nAdjPerEn,
                                            l_enToEn[0],
                                            l_spType,
                                            l_charsAdj,
                                            l_spLink[0] );
  // check the result
  REQUIRE( l_spLink[0][0] == 2 );
  REQUIRE( l_spLink[0][1] == 1 );

  REQUIRE( l_spLink[1][0] == 3 );
  REQUIRE( l_spLink[1][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink[2][0] == std::numeric_limits< int >::max() );
  REQUIRE( l_spLink[2][1] == 1 );

  REQUIRE( l_spLink[3][0] == 0 );
  REQUIRE( l_spLink[3][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink[4][0] == 1 );
  REQUIRE( l_spLink[4][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink[5][0] == 3 );
  REQUIRE( l_spLink[5][1] == 0 );
}

TEST_CASE( "SparseEnts: Links sparse entities based on adjacency information (double sparse type).", "[linkSpAdjDst][SparseEnts]" ) {
  /*
   * Test case for the flat version:
   *
   *   fa | el  | spId           el | spId         fa | spId           faSpId | elSpId
   *   0  | 1-4 |                0  | -            0  | x-x            0      | 2-1
   *   1  | 5-3 |  0             1  | -            1  | 2-1            1      | 3-x
   *   2  | 9-7 |  1             2  | 0            2  | 3-x            2      | x-1
   *   3  | 8-3 |  2             3  | 1            3  | x-1            3      | 0-x
   *   4  | 2-1 |  3             4  | -            4  | 0-x            4      | 1-x
   *   5  | 2-6 | ***            5  | 2            5  | x-x            5      | 3-0
   *   6  | 3-8 |  4             6  | -            6  | 1-x
   *   7  | 0-4 |                7  | -            7  | x-x
   *   8  | 9-2 |  5             8  | -            8  | 3-0
   *   9  | 4-1 |                9  | 3            9  | x-x
   *   10 | 4-9 | ***                              10 | x-x
   */

  // set input info
  int l_nEn1 = 11;
  int l_nAdjPerEn1 = 2;

  int l_enEn1[11][2] = { {1,4}, {5,3}, {9,7}, {8,3}, {2,1},
                         {2,6}, {3,8}, {0,4}, {9,2}, {4,1},
                         {4,9} };

  unsigned short l_spTypeFrom1 = 87;
  unsigned short l_spTypeTo1   = 43;

  typedef struct{ unsigned short spType; } CharsAdj;
  CharsAdj l_charsFrom1[11];
  l_charsFrom1[0].spType  = 5;
  l_charsFrom1[1].spType  = l_spTypeFrom1;
  l_charsFrom1[2].spType   = l_spTypeFrom1;
  l_charsFrom1[3].spType  = l_spTypeFrom1;
  l_charsFrom1[4].spType  = l_spTypeFrom1;
  l_charsFrom1[5].spType  = 9;
  l_charsFrom1[6].spType  = l_spTypeFrom1;
  l_charsFrom1[7].spType  = 5;
  l_charsFrom1[8].spType  = l_spTypeFrom1;
  l_charsFrom1[9].spType  = 4;
  l_charsFrom1[10].spType = 3;

  CharsAdj l_charsTo1[10];
  l_charsTo1[0].spType = 5;
  l_charsTo1[1].spType = 8;
  l_charsTo1[2].spType = l_spTypeTo1;
  l_charsTo1[3].spType = l_spTypeTo1;
  l_charsTo1[4].spType = 5;
  l_charsTo1[5].spType = l_spTypeTo1;
  l_charsTo1[6].spType = 0;
  l_charsTo1[7].spType = 5;
  l_charsTo1[8].spType = 3;
  l_charsTo1[9].spType = l_spTypeTo1;

  int l_spLink1[6][2];

  // link the sparse entities
  edge::data::SparseEntities::linkSpAdjDst( l_nEn1,
                                            l_nAdjPerEn1,
                                            l_enEn1[0],
                                            l_spTypeFrom1,
                                            l_spTypeTo1,
                                            l_charsFrom1,
                                            l_charsTo1,
                                            l_spLink1[0] );

  // check the result
  REQUIRE( l_spLink1[0][0] == 2 );
  REQUIRE( l_spLink1[0][1] == 1 );

  REQUIRE( l_spLink1[1][0] == 3 );
  REQUIRE( l_spLink1[1][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink1[2][0] == std::numeric_limits< int >::max() );
  REQUIRE( l_spLink1[2][1] == 1 );

  REQUIRE( l_spLink1[3][0] == 0 );
  REQUIRE( l_spLink1[3][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink1[4][0] == 1 );
  REQUIRE( l_spLink1[4][1] == std::numeric_limits< int >::max() );

  REQUIRE( l_spLink1[5][0] == 3 );
  REQUIRE( l_spLink1[5][1] == 0 );

  /*
   * Test case for the dynamic version.
   * In comparison to the flat version, the number of adjacent entities is dynamic and no undefined links are returned.
   *
   *   ve | el        | spId           el | spId         ve | spId           veSpId | elSpId
   *   0  | 1-4       |  -             0  | -            0  | x-x            0      | 2-1
   *   1  | 8-5-0-3   |  0             1  | -            1  | x-2-x-1        1      | -
   *   2  | 8-7-4     |  1             2  | 0            2  | x-x-x          2      | 1
   *   3  | 8-3-6-7-8 |  2             3  | 1            3  | x-1-x-x-x      3      | 0
   *   4  | 2-1       |  3             4  | -            4  | 0-x            4      | 1
   *   5  | 2-6-5-3   |  -             5  | 2            5  | x-x-x-x        5      | 3-0-2-1
   *   6  | 3-8       |  4             6  | -            6  | 1-x
   *   7  | 0-4       |  -             7  | -            7  | x-x
   *   8  | 9-2-5-3-7 |  5             8  | -            8  | 3-0-2-1-x
   *   9  | 4-1-9-8   |  -             9  | 3            9  | x-x-x-x
   *   10 | 4-9       |  -                               10 | x-x
   */

  // set input info
  int l_nEn2 = 11;

  int l_enEnRaw2[35] = { 1, 4,          //  0
                         8, 5, 0, 3,    //  2
                         8, 7, 4,       //  6
                         8, 3, 6, 7, 8, //  9
                         2, 1,          // 14
                         2, 6, 5, 3,    // 16
                         3, 8,          // 20
                         0, 4,          // 22
                         9, 2, 5, 3, 7, // 24
                         4, 1, 9, 8,    // 29
                         4, 9 };        // 33

  int *l_enEnPtr2[12] = { l_enEnRaw2+ 0,
                          l_enEnRaw2+ 2,
                          l_enEnRaw2+ 6,
                          l_enEnRaw2+ 9,
                          l_enEnRaw2+14,
                          l_enEnRaw2+16,
                          l_enEnRaw2+20,
                          l_enEnRaw2+22,
                          l_enEnRaw2+24,
                          l_enEnRaw2+29,
                          l_enEnRaw2+33,
                          l_enEnRaw2+35 }; // ghost

  unsigned short l_spTypeFrom2 = 129;
  unsigned short l_spTypeTo2   =  97;

  CharsAdj l_charsFrom2[11];
  l_charsFrom2[0].spType  = 5;
  l_charsFrom2[1].spType  = l_spTypeFrom2;
  l_charsFrom2[2].spType  = l_spTypeFrom2;
  l_charsFrom2[3].spType  = l_spTypeFrom2;
  l_charsFrom2[4].spType  = l_spTypeFrom2;
  l_charsFrom2[5].spType  = 9;
  l_charsFrom2[6].spType  = l_spTypeFrom2;
  l_charsFrom2[7].spType  = 5;
  l_charsFrom2[8].spType  = l_spTypeFrom2;
  l_charsFrom2[9].spType  = 4;
  l_charsFrom2[10].spType = 3;

  CharsAdj l_charsTo2[10];
  l_charsTo2[0].spType = 5;
  l_charsTo2[1].spType = 8;
  l_charsTo2[2].spType = l_spTypeTo2;
  l_charsTo2[3].spType = l_spTypeTo2;
  l_charsTo2[4].spType = 5;
  l_charsTo2[5].spType = l_spTypeTo2;
  l_charsTo2[6].spType = 0;
  l_charsTo2[7].spType = 5;
  l_charsTo2[8].spType = 3;
  l_charsTo2[9].spType = l_spTypeTo2;

  // derive sizes of link
  int l_nRaw2 = -1;
  int l_nPtr2 = -1;

  // test derivation of number of entries
  edge::data::SparseEntities::nLinkSpAdjDst( l_nEn2,
                                             l_enEnPtr2,
                                             l_spTypeFrom2,
                                             l_spTypeTo2,
                                             l_charsFrom2,
                                             l_charsTo2,
                                             l_nRaw2,
                                             l_nPtr2 );

  // check the result
  REQUIRE( l_nRaw2 == 9 );
  REQUIRE( l_nPtr2 == 6 );

  // derive the link itself
  int l_raw2[9];
  int *l_ptr2[6+1];

  // test derivation of sparse link
  edge::data::SparseEntities::linkSpAdjDst( l_nEn2,
                                            l_enEnPtr2,
                                            l_spTypeFrom2,
                                            l_spTypeTo2,
                                            l_charsFrom2,
                                            l_charsTo2,
                                            l_raw2,
                                            l_ptr2 );


  /*
   * Test case for the dynamic version.
   * In comparison to the flat version, the number of adjacent entities is dynamic and no undefined links are returned.
   *
   *   ve | el        | spId           el | spId         ve | spId           veSpId | elSpId
   *   0  | 1-4       |  -             0  | -            0  | x-x            0      | 2-1
   *   1  | 8-5-0-3   |  0             1  | -            1  | x-2-x-1        1      | -
   *   2  | 8-7-4     |  1             2  | 0            2  | x-x-x          2      | 1
   *   3  | 8-3-6-7-8 |  2             3  | 1            3  | x-1-x-x-x      3      | 0
   *   4  | 2-1       |  3             4  | -            4  | 0-x            4      | 1
   *   5  | 2-6-5-3   |  -             5  | 2            5  | x-x-x-x        5      | 3-0-2-1
   *   6  | 3-8       |  4             6  | -            6  | 1-x
   *   7  | 0-4       |  -             7  | -            7  | x-x
   *   8  | 9-2-5-3-7 |  5             8  | -            8  | 3-0-2-1-x
   *   9  | 4-1-9-8   |  -             9  | 3            9  | x-x-x-x
   *   10 | 4-9       |  -                               10 | x-x
   */

  // check the results
  REQUIRE( l_raw2[0] == 2 );
  REQUIRE( l_raw2[1] == 1 );

  REQUIRE( l_raw2[2] == 1 );

  REQUIRE( l_raw2[3] == 0 );

  REQUIRE( l_raw2[4] == 1 );

  REQUIRE( l_raw2[5] == 3 );
  REQUIRE( l_raw2[6] == 0 );
  REQUIRE( l_raw2[7] == 2 );
  REQUIRE( l_raw2[8] == 1 );

  REQUIRE( l_ptr2[0] == l_raw2+0 );
  REQUIRE( l_ptr2[1] == l_raw2+2 );
  REQUIRE( l_ptr2[2] == l_raw2+2 );
  REQUIRE( l_ptr2[3] == l_raw2+3 );
  REQUIRE( l_ptr2[4] == l_raw2+4 );
  REQUIRE( l_ptr2[5] == l_raw2+5 );
  REQUIRE( l_ptr2[6] == l_raw2+9 );
}

TEST_CASE( "SparseEnts: Inherit sparse information within entities.", "[inherit][SparseEnts]" ) {
  /*
   * Out setup:
   *   "parent" sparse type: 01
   *   "child"  sparse type: 10
   *
   *   Input:      Output: 
   *   en | bits   en | bits
   *   0  | 01     0  | 11
   *   1  | 11     1  | 11
   *   2  | 01     2  | 11
   *   3  | 10     3  | 10
   *   4  | 01     4  | 11
   *   5  | 00     5  | 00
   */

  // setup chars
  typedef struct { unsigned short spType; } t_enChars;

  t_enChars l_enChars1[6] = { {1}, {3}, {1}, {2}, {1}, {0} };

  // call inheritance
  edge::data::SparseEntities::inherit( 6,
                                       1,
                                       2,
                                       l_enChars1 );

  // check the result
  REQUIRE( l_enChars1[0].spType == 3 );
  REQUIRE( l_enChars1[1].spType == 3 );
  REQUIRE( l_enChars1[2].spType == 3 );
  REQUIRE( l_enChars1[3].spType == 2 );
  REQUIRE( l_enChars1[4].spType == 3 );
  REQUIRE( l_enChars1[5].spType == 0 );
}

TEST_CASE( "SparseEnts: Propagate sparse information to adjacent entities.", "[propAdj][SparseEnts]" ) {
  /*
   * Our setup for the flat array version
   *
   *   Input bits  Input bits  Adjacency of    Result (critical bit is
   *   of en0:     of en1:     entity 0 to 1:  the second (x): xy:
   *
   *   en0 | bits  en0 | bits  en0 | en1       en1
   *   0   | 00    0   | 00    0   | 0-2       0 | 00
   *   1   | 11    1   | 01    1   | 1-2       1 | 11
   *   2   | 00    2   | 00    2   | 3-1       2 | 10
   *   3   | 10    3   | 01    3   | 4-2       3 | 01
   *   4   | 11    4   | 01    4   | 1-4       4 | 11
   *   5   | 00                5   | 2-3
   */

  typedef struct { unsigned short spType; } t_enChars;

  // set up en0-chars
  t_enChars l_en0Chars1[6];
  l_en0Chars1[0].spType = 0;
  l_en0Chars1[1].spType = 3;
  l_en0Chars1[2].spType = 0;
  l_en0Chars1[3].spType = 2;
  l_en0Chars1[4].spType = 3;
  l_en0Chars1[5].spType = 0;

  // set up en1-chars
  t_enChars l_en1Chars1[5];
  l_en1Chars1[0].spType = 0;
  l_en1Chars1[1].spType = 1;
  l_en1Chars1[2].spType = 0;
  l_en1Chars1[3].spType = 1;
  l_en1Chars1[4].spType = 1;

  // set up adjacency info
  int l_en0ToEn11[6][2] = { {0,2}, {1,2}, {3,1}, {4,2}, {1,4}, {2,3} };

  // call the propagation
  edge::data::SparseEntities::propAdj( 6,
                                       2,
                                       l_en0ToEn11[0],
                                       2,
                                       2,
                                       l_en0Chars1,
                                       l_en1Chars1 );

  // check the result
  REQUIRE( l_en1Chars1[0].spType == 0 );
  REQUIRE( l_en1Chars1[1].spType == 3 );
  REQUIRE( l_en1Chars1[2].spType == 2 );
  REQUIRE( l_en1Chars1[3].spType == 1 );
  REQUIRE( l_en1Chars1[4].spType == 3 );


  /*
   * Our setup for the array of pointer version
   *
   *   Input bits  Input bits  Adjacency of    Result (critical bit is
   *   of en0:     of en1:     entity 0 to 1:  the second (x): xy:
   *
   *   en0 | bits  en0 | bits  en0 | en1       en1
   *   0   | 00    0   | 00    0   | 0         0 | 10
   *   1   | 11    1   | 01    1   | 1-2       1 | 11
   *   2   | 00    2   | 00    2   | 3-1-2     2 | 10
   *   3   | 10    3   | 01    3   | 4-2       3 | 01
   *   4   | 11    4   | 01    4   | 1-4-0     4 | 11
   *   5   | 00                5   | 2-3-0-1
   */

  // set up en0-chars
  t_enChars l_en0Chars2[6];
  l_en0Chars2[0].spType = 0;
  l_en0Chars2[1].spType = 3;
  l_en0Chars2[2].spType = 0;
  l_en0Chars2[3].spType = 2;
  l_en0Chars2[4].spType = 3;
  l_en0Chars2[5].spType = 0;

  // set up en1-chars
  t_enChars l_en1Chars2[5];
  l_en1Chars2[0].spType = 0;
  l_en1Chars2[1].spType = 1;
  l_en1Chars2[2].spType = 0;
  l_en1Chars2[3].spType = 1;
  l_en1Chars2[4].spType = 1;

  // set up adjacency info   0   1     3       6     8       11      15
  int l_en0ToEn12[15] =    { 0,  1,2,  3,1,2,  4,2,  1,4,0,  2,3,0,1 };
  int *l_ptrs[7] = { l_en0ToEn12+0, l_en0ToEn12+1, l_en0ToEn12+3, l_en0ToEn12+6, l_en0ToEn12+8, l_en0ToEn12+11, l_en0ToEn12+15 };

  // call the propagation
  edge::data::SparseEntities::propAdj( 6,
                                       l_ptrs,
                                       2,
                                       2,
                                       l_en0Chars2,
                                       l_en1Chars2 );

  // check the result
  REQUIRE( l_en1Chars2[0].spType == 2 );
  REQUIRE( l_en1Chars2[1].spType == 3 );
  REQUIRE( l_en1Chars2[2].spType == 2 );
  REQUIRE( l_en1Chars2[3].spType == 1 );
  REQUIRE( l_en1Chars2[4].spType == 3 );
}

TEST_CASE( "Sparse Ents: Derivation of dense ids from Cartesian points", "[ptToMesh][SparseEnts]") {
  t_enLayout l_elLayout;

  l_elLayout.timeGroups.resize( 2 );
  l_elLayout.timeGroups[0].nEntsOwn    = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[1].nEntsOwn    = 3;
  l_elLayout.timeGroups[1].nEntsNotOwn = 4;

  // setup receiver coordinates
  real_mesh l_recvCrds[6][3] = { { -1.00, -1.00, -2.00 }, // not part of the "domain"
                                 {  0.15,  0.15,  0.15 }, // inside a tet
                                 {  0.00,  0.05,  0.00 }, // hits an edge
                                 {  0.00,  0.00,  0.00 }, // vertex
                                 {  2.00,  0.00,  0.00 }, // outside
                                 {  0.20,  0.10,  0.05 }  // inside
                               };

  t_vertexChars l_veChars[8] = { {{0.0, 0.0, 0.0}, 0}, // reference tet as dummy tet
                                 {{1.0, 0.0, 0.0}, 0},
                                 {{0.0, 1.0, 0.0}, 0},
                                 {{0.0, 0.0, 1.0}, 0}, 
                                 {{0.0, 0.0, 8.0}, 0}, // shifted ref-tet
                                 {{1.0, 0.0, 8.0}, 0},
                                 {{0.0, 1.0, 8.0}, 0},
                                 {{0.0, 1.0, 9.0}, 0} };

  // setup tets
  int l_enVe[10][4] = { {4,5,6,7}, // tg 1
                        {4,5,6,7},
                        {0,0,0,0},
                        {4,5,6,7}, // tg 2
                        {0,1,2,3},
                        {0,1,2,3},
                        {0,1,2,3}, // ghost tg 2
                        {4,5,6,7},
                        {4,5,6,7},
                        {0,1,2,3} };

  // For plotting in Mathematica, you can use:
  //  Graphics3D[ {{Opacity[.3], Tetrahedron[]}, {PointSize[0.1], Point[{0.15, 0.15, 0.15}] }}, Axes -> True ]

  // dense ids
  int l_de[6];

  edge::data::SparseEntities::ptToEn( TET4,
                                      6,
                                      l_recvCrds,
                                      l_elLayout,
                                      l_enVe[0],
                                      l_veChars,
                                      l_de );

  REQUIRE( l_de[0] == 4 );
}
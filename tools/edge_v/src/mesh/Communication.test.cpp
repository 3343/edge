/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Tests the derivation of communication structures.
 **/
#include <catch.hpp>
#define private public
#include "Communication.h"
#undef private

TEST_CASE( "Tests the derivation of the partition-time group pairs.", "[communication][getPaTgPairs]" ) {
  /*
   * Our test case:
   *
   *  id | elFaEl   | partition | time group | pairs
   *   0 |          | 0         | 2          |
   *   1 |          | 0         | 2          |
   *   2 |          | 1         | 0          |
   *  ==============|===========|============|
   *   3 |  6  4  2 | 2         | 1          | (1, 0)
   *   4 |  1  5  9 | 2         | 1          | (0, 2), (4, 2)
   *   5 |  7  2 12 | 2         | 1          | (1, 0), (5, 0)
   *   6 | 10  5  4 | 2         | 1          | (4, 1)
   *   7 |  3 11  0 | 2         | 1          | (5, 2), (0, 2)
   *  ==============|===========|============|
   *   8 |          | 3         | 1          |
   *   9 |          | 4         | 2          |
   *  10 |          | 4         | 1          |
   *  11 |          | 5         | 2          |
   *  12 |          | 5         | 0
   *
   * Sorted, unique pairs: (0, 2), (1, 0), (4, 1), (4, 2), (5, 0), (5, 2)
   */

  // assemble test structures
  std::size_t l_elFaEl[13][3];
  for( unsigned short l_el = 0; l_el < 13; l_el++ )
    for( unsigned short l_fa = 0; l_fa < 3; l_fa++ )
      l_elFaEl[l_el][l_fa] = std::numeric_limits< std::size_t >::max();

  l_elFaEl[3][0] =  6; l_elFaEl[3][1] =  4; l_elFaEl[3][2] =  2;
  l_elFaEl[4][0] =  1; l_elFaEl[4][1] =  5; l_elFaEl[4][2] =  9;
  l_elFaEl[5][0] =  7; l_elFaEl[5][1] =  2; l_elFaEl[5][2] = 12;
  l_elFaEl[6][0] = 10; l_elFaEl[6][1] =  5; l_elFaEl[6][2] =  4;
  l_elFaEl[7][0] =  3; l_elFaEl[7][1] = 11; l_elFaEl[7][2] =  0;

  //                            0  1  2  3  4  5  6  7  8  9 10 11 12
  std::size_t l_elPa[13]    = { 0, 0, 1, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5 };
  unsigned short l_elTg[13] = { 2, 2, 0, 1, 1, 1, 1, 1, 1, 2, 1, 2, 0 };

  // get the unique and sort pairs
  std::set< std::pair< std::size_t, unsigned short > > l_pairs;
  edge_v::mesh::Communication::getPaTgPairs( edge_v::TRIA3,
                                             3,
                                             5,
                                             l_elFaEl[0],
                                             l_elTg,
                                             l_elPa,
                                             l_pairs );


  // check the results
  auto l_it = l_pairs.begin();

  REQUIRE( std::get<0>(*l_it) == 0 );
  REQUIRE( std::get<1>(*l_it) == 2 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 1 );
  REQUIRE( std::get<1>(*l_it) == 0 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 4 );
  REQUIRE( std::get<1>(*l_it) == 1 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 4 );
  REQUIRE( std::get<1>(*l_it) == 2 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 5 );
  REQUIRE( std::get<1>(*l_it) == 0 );

  l_it++;
  REQUIRE( std::get<0>(*l_it) == 5 );
  REQUIRE( std::get<1>(*l_it) == 2 );
}

TEST_CASE( "Tests the derivation of communicating elements.", "[communication][getPaElComm]" ) {
  /*
   * Test case:
   *
   * id | partition | send? | elFaEl
   *  0 | 0         | 0     | 0, 2, max
   *  1 | 0         | 0     | max, 2, 0
   *  2 | 0         | 1     | 14, 7, 1
   *  3 | 1         | 0     | 4, 5, max
   *  4 | 1         | 0     | 3, 3, 3
   *  5 | 1         | 1     | 6, 7, 8
   *  6 | 1         | 1     | 7, 3, 0
   *  7 | 1         | 1     | 15, 8, 2
   *  8 | 2         | 1     | 7, 6, 9
   *  9 | 2         | 1     | 10, 8, 2
   * 10 | 3         | 0     | 11, 12, 14
   * 11 | 3         | 0     | 15, 10, 13
   * 12 | 3         | 0     | 11, 13, 15
   * 13 | 3         | 1     | 12, 7, 6
   * 14 | 3         | 1     | 6, 8, 2
   * 15 | 3         | 1     | 14, 8, 9
   */
  //                         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  std::size_t l_elPa[16] = { 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3 };
  std::size_t l_elFaEl[16][3] = { {0, 2, std::numeric_limits< std::size_t >::max() },
                                  {std::numeric_limits< std::size_t >::max(), 2, 0},
                                  {14, 7, 1},
                                  {4, 5, std::numeric_limits< std::size_t >::max()},
                                  {3, 3, 3},
                                  {6, 7, 8},
                                  {7, 3, 0},
                                  {15, 8, 2},
                                  {7, 6, 9},
                                  {10, 8, 2},
                                  {11, 12, 14},
                                  {15, 10, 13},
                                  {11, 13, 15},
                                  {12, 7, 6},
                                  {6, 8, 2},
                                  {14, 8, 9} };
  std::size_t l_first[4] = {0};
  std::size_t l_size[4] = {0};
  edge_v::mesh::Communication::getPaElComm( edge_v::TRIA3,
                                            16,
                                            l_elFaEl[0],
                                            l_elPa,
                                            l_first,
                                            l_size );

  // check the results
  REQUIRE( l_first[0] ==  2 );
  REQUIRE( l_size[0]  ==  1 );

  REQUIRE( l_first[1] ==  5 );
  REQUIRE( l_size[1]  ==  3 );

  REQUIRE( l_first[2] ==  8 );
  REQUIRE( l_size[2]  ==  2 );

  REQUIRE( l_first[3] == 13 );
  REQUIRE( l_size[3]  ==  3 );
}

TEST_CASE( "Tests the derivation of send messages for a single communication region.", "[communication][getMsgsSend]" ) {
  /*
   * Our test case:
   *
   *  id | elFaEl   | partition | time group | pairs
   *   0 |     7    | 0         | 2          |
   *   1 |  4       | 0         | 2          |
   *   2 |     5 3  | 1         | 0          |
   *  ==============|===========|============|
   *   3 |  6  4  2 | 2         | 1          | (1, 0)
   *   4 |  1  5  9 | 2         | 1          | (0, 2), (4, 2)
   *   5 |  7  2 12 | 2         | 1          | (1, 0), (5, 0)
   *   6 | 10  5  4 | 2         | 1          | (4, 1)
   *   7 |  3 11  0 | 2         | 1          | (5, 2), (0, 2)
   *  ==============|===========|============|
   *   8 |          | 3         | 1          |
   *   9 |  4       | 4         | 2          |
   *  10 |        6 | 4         | 1          |
   *  11 |        7 | 5         | 2          |
   *  12 |     5    | 5         | 0
   *
   * Sorted, unique pairs: (0, 2), (1, 0), (4, 1), (4, 2), (5, 0), (5, 2)
   */
  // assemble test structures
  std::size_t l_elFaEl[13][3];
  for( unsigned short l_el = 0; l_el < 13; l_el++ )
    for( unsigned short l_fa = 0; l_fa < 3; l_fa++ )
      l_elFaEl[l_el][l_fa] = std::numeric_limits< std::size_t >::max();

  l_elFaEl[3][0] =  6; l_elFaEl[3][1] =  4; l_elFaEl[3][2] =  2;
  l_elFaEl[4][0] =  1; l_elFaEl[4][1] =  5; l_elFaEl[4][2] =  9;
  l_elFaEl[5][0] =  7; l_elFaEl[5][1] =  2; l_elFaEl[5][2] = 12;
  l_elFaEl[6][0] = 10; l_elFaEl[6][1] =  5; l_elFaEl[6][2] =  4;
  l_elFaEl[7][0] =  3; l_elFaEl[7][1] = 11; l_elFaEl[7][2] =  0;

  l_elFaEl[ 0][1] = 7;
  l_elFaEl[ 1][0] = 4;
  l_elFaEl[ 2][1] = 5;
  l_elFaEl[ 2][2] = 3;
  l_elFaEl[ 9][0] = 4;
  l_elFaEl[10][2] = 6;
  l_elFaEl[11][2] = 7;
  l_elFaEl[12][1] = 5;

  //                            0  1  2  3  4  5  6  7  8  9 10 11 12
  std::size_t l_elPa[13]    = { 0, 0, 1, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5 };
  unsigned short l_elTg[13] = { 2, 2, 0, 1, 1, 1, 1, 1, 1, 2, 1, 2, 0 };

  // get the send messages
  std::vector< edge_v::mesh::Communication::Message > l_send;
  edge_v::mesh::Communication::getMsgsSend( edge_v::TRIA3,
                                            3,
                                            5,
                                            l_elFaEl[0],
                                            l_elPa,
                                            l_elTg,
                                            l_send );

  // check the results
  REQUIRE( l_send.size() == 6 );

  REQUIRE( l_send[0].pa == 0 );
  REQUIRE( l_send[0].tg == 2 );
  REQUIRE( l_send[0].el.size() == 2 );
  REQUIRE( l_send[0].fa.size() == 2 );
  REQUIRE( l_send[0].el[0] == 4 );
  REQUIRE( l_send[0].fa[0] == 0 );
  REQUIRE( l_send[0].el[1] == 7 );
  REQUIRE( l_send[0].fa[1] == 2 );

  REQUIRE( l_send[1].pa == 1 );
  REQUIRE( l_send[1].tg == 0 );
  REQUIRE( l_send[1].el.size() == 2 );
  REQUIRE( l_send[1].fa.size() == 2 );
  REQUIRE( l_send[1].el[0] == 3 );
  REQUIRE( l_send[1].fa[0] == 2 );
  REQUIRE( l_send[1].el[1] == 5 );
  REQUIRE( l_send[1].fa[1] == 1 );

  REQUIRE( l_send[2].pa == 4 );
  REQUIRE( l_send[2].tg == 1 );
  REQUIRE( l_send[2].el.size() == 1 );
  REQUIRE( l_send[2].fa.size() == 1 );
  REQUIRE( l_send[2].el[0] == 6 );
  REQUIRE( l_send[2].fa[0] == 0 );

  REQUIRE( l_send[3].pa == 4 );
  REQUIRE( l_send[3].tg == 2 );
  REQUIRE( l_send[3].el.size() == 1 );
  REQUIRE( l_send[3].fa.size() == 1 );
  REQUIRE( l_send[3].el[0] == 4 );
  REQUIRE( l_send[3].fa[0] == 2 );

  REQUIRE( l_send[4].pa == 5 );
  REQUIRE( l_send[4].tg == 0 );
  REQUIRE( l_send[4].el.size() == 1 );
  REQUIRE( l_send[4].fa.size() == 1 );
  REQUIRE( l_send[4].el[0] == 5 );
  REQUIRE( l_send[4].fa[0] == 2 );

  REQUIRE( l_send[5].pa == 5 );
  REQUIRE( l_send[5].tg == 2 );
  REQUIRE( l_send[5].el.size() == 1 );
  REQUIRE( l_send[5].fa.size() == 1 );
  REQUIRE( l_send[5].el[0] == 7 );
  REQUIRE( l_send[5].fa[0] == 1 );
}

TEST_CASE( "Tests the derivation of the global communication structure.", "[communication][getGlobal]" ) {
  /*
   * Our test case:
   *
   *  id | elFaEl   | partition | time group | comm partners
   *   0 |  7  2  1 | 0         | 2          | (2, 3), (1, 0)
   *   1 |  0  2  3 | 0         | 2          | (1, 0), (2, 1)
   *   2 |  0  1  6 | 1         | 0          | (0, 2), (2, 2)
   *   3 |  8  1  4 | 2         | 1          | (3, 1), (0, 2)
   *   4 |  11 5  7 | 2         | 2          | (5, 2)
   *   5 |  4  11 6 | 2         | 2          | (5, 2)
   *   6 | 12  2  7 | 2         | 2          | (5, 4), (1, 0)
   *   7 |  3 11  0 | 2         | 3          | (5, 2), (0, 2)
   *   8 |  3  8  8 | 3         | 1          | (2, 1)
   *   9 | 10  9  9 | 4         | 0          |
   *  10 |  9 10 12 | 4         | 1          | (5, 0)
   *  11 |  5  4  7 | 5         | 2          | (2, 2), (2, 2), (2, 3)
   *  12 |  6 11 10 | 5         | 4          | (2, 2), (4, 1)
   *
   */
  // assemble test structures
  std::size_t l_elFaEl[13][3] = { {7,  2,  1}, {  0,  2,  3}, { 0,  1,  6}, {8, 1, 4}, {11, 5, 7},
                                  {4, 11,  6}, { 12,  2,  7}, { 3, 11,  0}, {3, 8, 8}, {10, 9, 9},
                                  {9, 10, 12}, {  5,  4,  7}, { 6, 11, 10} };

  //                            0  1  2  3  4  5  6  7  8  9 10 11 12
  std::size_t l_elPa[13]    = { 0, 0, 1, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5 };
  unsigned short l_elTg[13] = { 2, 2, 0, 1, 2, 2, 2, 3, 1, 0, 1, 2, 4 };


  std::vector< edge_v::mesh::Communication::Partition > l_struct;
  edge_v::mesh::Communication::getStruct( edge_v::TRIA3,
                                          13,
                                          l_elFaEl[0],
                                          l_elPa,
                                          l_elTg,
                                          l_struct );

  // check the results
  REQUIRE( l_struct[0].tr.size() == 1 );
  REQUIRE( l_struct[0].tr[0].tg == 2 );
  REQUIRE( l_struct[0].tr[0].send.size() == 3 );
  REQUIRE( l_struct[0].tr[0].recv.size() == 3 );

  REQUIRE( l_struct[0].tr[0].send[0].pa == 1 );
  REQUIRE( l_struct[0].tr[0].send[0].tg == 0 );
  REQUIRE( l_struct[0].tr[0].send[0].el.size() == 2 );
  REQUIRE( l_struct[0].tr[0].send[0].fa.size() == 2 );
  REQUIRE( l_struct[0].tr[0].send[0].el[0] == 0 );
  REQUIRE( l_struct[0].tr[0].send[0].fa[0] == 1 );
  REQUIRE( l_struct[0].tr[0].send[0].elAd[0] == 2 );
  REQUIRE( l_struct[0].tr[0].send[0].faAd[0] == 0 );
  REQUIRE( l_struct[0].tr[0].send[0].el[1] == 1 );
  REQUIRE( l_struct[0].tr[0].send[0].fa[1] == 1 );
  REQUIRE( l_struct[0].tr[0].send[0].elAd[1] == 2 );
  REQUIRE( l_struct[0].tr[0].send[0].faAd[1] == 1 );
  REQUIRE( l_struct[0].tr[0].recv[0].pa == 1 );
  REQUIRE( l_struct[0].tr[0].recv[0].tg == 0 );
  REQUIRE( l_struct[0].tr[0].recv[0].el.size() == 2 );
  REQUIRE( l_struct[0].tr[0].recv[0].fa.size() == 2 );
  REQUIRE( l_struct[0].tr[0].recv[0].elAd[0] == 2 );
  REQUIRE( l_struct[0].tr[0].recv[0].faAd[0] == 0 );
  REQUIRE( l_struct[0].tr[0].recv[0].elAd[1] == 2 );
  REQUIRE( l_struct[0].tr[0].recv[0].faAd[1] == 1 );

  REQUIRE( l_struct[0].tr[0].send[1].pa == 2 );
  REQUIRE( l_struct[0].tr[0].send[1].tg == 1 );
  REQUIRE( l_struct[0].tr[0].send[1].el.size() == 1 );
  REQUIRE( l_struct[0].tr[0].send[1].fa.size() == 1 );
  REQUIRE( l_struct[0].tr[0].send[1].el[0] == 1 );
  REQUIRE( l_struct[0].tr[0].send[1].fa[0] == 2 );

  REQUIRE( l_struct[0].tr[0].send[2].pa == 2 );
  REQUIRE( l_struct[0].tr[0].send[2].tg == 3 );
  REQUIRE( l_struct[0].tr[0].send[2].el.size() == 1 );
  REQUIRE( l_struct[0].tr[0].send[2].fa.size() == 1 );
  REQUIRE( l_struct[0].tr[0].send[2].el[0] == 0 );
  REQUIRE( l_struct[0].tr[0].send[2].fa[0] == 0 );


  REQUIRE( l_struct[1].tr.size() == 1 );
  REQUIRE( l_struct[1].tr[0].tg == 0 );
  REQUIRE( l_struct[1].tr[0].send.size() == 2 );
  REQUIRE( l_struct[1].tr[0].recv.size() == 2 );

  REQUIRE( l_struct[1].tr[0].send[0].pa == 0 );
  REQUIRE( l_struct[1].tr[0].send[0].tg == 2 );
  REQUIRE( l_struct[1].tr[0].send[0].el.size() == 2 );
  REQUIRE( l_struct[1].tr[0].send[0].fa.size() == 2 );
  REQUIRE( l_struct[1].tr[0].send[0].el[0] == 2 );
  REQUIRE( l_struct[1].tr[0].send[0].fa[0] == 0 );
  REQUIRE( l_struct[1].tr[0].send[0].el[1] == 2 );
  REQUIRE( l_struct[1].tr[0].send[0].fa[1] == 1 );

  REQUIRE( l_struct[2].tr.size() == 3 );
  REQUIRE( l_struct[2].tr[0].tg == 1 );
  REQUIRE( l_struct[2].tr[0].send.size() == 2 );
  REQUIRE( l_struct[2].tr[0].recv.size() == 2 );

  REQUIRE( l_struct[2].tr[0].send[0].pa == 0 );
  REQUIRE( l_struct[2].tr[0].send[0].tg == 2 );
  REQUIRE( l_struct[2].tr[0].send[0].el.size() == 1 );
  REQUIRE( l_struct[2].tr[0].send[0].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[0].send[0].el[0] == 3 );
  REQUIRE( l_struct[2].tr[0].send[0].fa[0] == 1 );

  REQUIRE( l_struct[2].tr[0].send[1].pa == 3 );
  REQUIRE( l_struct[2].tr[0].send[1].tg == 1 );
  REQUIRE( l_struct[2].tr[0].send[1].el.size() == 1 );
  REQUIRE( l_struct[2].tr[0].send[1].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[0].send[1].el[0] == 3 );
  REQUIRE( l_struct[2].tr[0].send[1].fa[0] == 0 );


  REQUIRE( l_struct[2].tr[1].tg == 2 );
  REQUIRE( l_struct[2].tr[1].send.size() == 3 );
  REQUIRE( l_struct[2].tr[1].recv.size() == 3 );

  REQUIRE( l_struct[2].tr[1].send[0].pa == 1 );
  REQUIRE( l_struct[2].tr[1].send[0].tg == 0 );
  REQUIRE( l_struct[2].tr[1].send[0].el.size() == 1 );
  REQUIRE( l_struct[2].tr[1].send[0].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[1].send[0].el[0] == 6 );
  REQUIRE( l_struct[2].tr[1].send[0].fa[0] == 1 );

  REQUIRE( l_struct[2].tr[1].send[1].pa == 5 );
  REQUIRE( l_struct[2].tr[1].send[1].tg == 2 );
  REQUIRE( l_struct[2].tr[1].send[1].el.size() == 2 );
  REQUIRE( l_struct[2].tr[1].send[1].fa.size() == 2 );
  REQUIRE( l_struct[2].tr[1].send[1].el[0] == 4 );
  REQUIRE( l_struct[2].tr[1].send[1].fa[0] == 0 );
  REQUIRE( l_struct[2].tr[1].send[1].el[1] == 5 );
  REQUIRE( l_struct[2].tr[1].send[1].fa[1] == 1 );

  REQUIRE( l_struct[2].tr[1].send[2].pa == 5 );
  REQUIRE( l_struct[2].tr[1].send[2].tg == 4 );
  REQUIRE( l_struct[2].tr[1].send[2].el.size() == 1 );
  REQUIRE( l_struct[2].tr[1].send[2].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[1].send[2].el[0] == 6 );
  REQUIRE( l_struct[2].tr[1].send[2].fa[0] == 0 );

  REQUIRE( l_struct[2].tr[2].tg == 3 );
  REQUIRE( l_struct[2].tr[2].send.size() == 2 );
  REQUIRE( l_struct[2].tr[2].recv.size() == 2 );

  REQUIRE( l_struct[2].tr[2].send[0].pa == 0 );
  REQUIRE( l_struct[2].tr[2].send[0].tg == 2 );
  REQUIRE( l_struct[2].tr[2].send[0].el.size() == 1 );
  REQUIRE( l_struct[2].tr[2].send[0].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[2].send[0].el[0] == 7 );
  REQUIRE( l_struct[2].tr[2].send[0].fa[0] == 2 );

  REQUIRE( l_struct[2].tr[2].send[1].pa == 5 );
  REQUIRE( l_struct[2].tr[2].send[1].tg == 2 );
  REQUIRE( l_struct[2].tr[2].send[1].el.size() == 1 );
  REQUIRE( l_struct[2].tr[2].send[1].fa.size() == 1 );
  REQUIRE( l_struct[2].tr[2].send[1].el[0] == 7 );
  REQUIRE( l_struct[2].tr[2].send[1].fa[0] == 1 );


  REQUIRE( l_struct[3].tr.size() == 1 );
  REQUIRE( l_struct[3].tr[0].tg == 1 );
  REQUIRE( l_struct[3].tr[0].send.size() == 1 );
  REQUIRE( l_struct[3].tr[0].recv.size() == 1 );

  REQUIRE( l_struct[3].tr[0].send[0].pa == 2 );
  REQUIRE( l_struct[3].tr[0].send[0].tg == 1 );
  REQUIRE( l_struct[3].tr[0].send[0].el.size() == 1 );
  REQUIRE( l_struct[3].tr[0].send[0].fa.size() == 1 );
  REQUIRE( l_struct[3].tr[0].send[0].el[0] == 8 );
  REQUIRE( l_struct[3].tr[0].send[0].fa[0] == 0 );


  REQUIRE( l_struct[4].tr.size() == 1 );
  REQUIRE( l_struct[4].tr[0].tg == 1 );
  REQUIRE( l_struct[4].tr[0].send.size() == 1 );
  REQUIRE( l_struct[4].tr[0].recv.size() == 1 );

  REQUIRE( l_struct[4].tr[0].send[0].pa == 5 );
  REQUIRE( l_struct[4].tr[0].send[0].tg == 4 );
  REQUIRE( l_struct[4].tr[0].send[0].el.size() == 1 );
  REQUIRE( l_struct[4].tr[0].send[0].fa.size() == 1 );
  REQUIRE( l_struct[4].tr[0].send[0].el[0] == 10 );
  REQUIRE( l_struct[4].tr[0].send[0].fa[0] ==  2 );


  REQUIRE( l_struct[5].tr.size() == 2 );
  REQUIRE( l_struct[5].tr[0].tg == 2 );
  REQUIRE( l_struct[5].tr[0].send.size() == 2 );
  REQUIRE( l_struct[5].tr[0].recv.size() == 2 );

  REQUIRE( l_struct[5].tr[0].send[0].pa == 2 );
  REQUIRE( l_struct[5].tr[0].send[0].tg == 2 );
  REQUIRE( l_struct[5].tr[0].send[0].el.size() == 2 );
  REQUIRE( l_struct[5].tr[0].send[0].fa.size() == 2 );
  REQUIRE( l_struct[5].tr[0].send[0].el[0] == 11 );
  REQUIRE( l_struct[5].tr[0].send[0].fa[0] ==  0 );
  REQUIRE( l_struct[5].tr[0].send[0].el[1] == 11 );
  REQUIRE( l_struct[5].tr[0].send[0].fa[1] ==  1 );

  REQUIRE( l_struct[5].tr[1].tg == 4 );
  REQUIRE( l_struct[5].tr[1].send.size() == 2 );
  REQUIRE( l_struct[5].tr[1].recv.size() == 2 );
  REQUIRE( l_struct[5].tr[1].send[0].el.size() == 1 );
  REQUIRE( l_struct[5].tr[1].send[0].fa.size() == 1 );
  REQUIRE( l_struct[5].tr[1].send[0].pa == 2 );
  REQUIRE( l_struct[5].tr[1].send[0].tg == 2 );
  REQUIRE( l_struct[5].tr[1].send[0].el[0] == 12 );
  REQUIRE( l_struct[5].tr[1].send[0].fa[0] ==  0 );
  REQUIRE( l_struct[5].tr[1].send[1].pa == 4 );
  REQUIRE( l_struct[5].tr[1].send[1].tg == 1 );
  REQUIRE( l_struct[5].tr[1].send[1].el.size() == 1 );
  REQUIRE( l_struct[5].tr[1].send[1].fa.size() == 1 );
  REQUIRE( l_struct[5].tr[1].send[1].el[0] == 12 );
  REQUIRE( l_struct[5].tr[1].send[1].fa[0] ==  2 );

  REQUIRE( l_struct[5].tr[1].recv[0].pa == 2 );
  REQUIRE( l_struct[5].tr[1].recv[0].tg == 2 );
  REQUIRE( l_struct[5].tr[1].recv[0].el.size() == 1 );
  REQUIRE( l_struct[5].tr[1].recv[0].fa.size() == 1 );
  REQUIRE( l_struct[5].tr[1].recv[0].elAd[0] == 6 );
  REQUIRE( l_struct[5].tr[1].recv[0].faAd[0] == 0 );

  REQUIRE( l_struct[5].tr[1].recv[1].pa == 4 );
  REQUIRE( l_struct[5].tr[1].recv[1].tg == 1 );
  REQUIRE( l_struct[5].tr[1].recv[1].el.size() == 1 );
  REQUIRE( l_struct[5].tr[1].recv[1].fa.size() == 1 );
  REQUIRE( l_struct[5].tr[1].recv[1].elAd[0] == 10 );
  REQUIRE( l_struct[5].tr[1].recv[1].faAd[0] ==  2 );

  // go through the constructor
  std::size_t l_nPaEls[6] = {2, 1, 5, 1, 2, 2};

  edge_v::mesh::Communication l_comm( 5,
                                      edge_v::TRIA3,
                                      13,
                                      l_elFaEl[0],
                                      6,
                                      l_nPaEls,
                                      l_elTg );

  // check simplified comm structure
  REQUIRE( l_comm.getStruct(0)[4*0+0+0] == 3 );

  REQUIRE( l_comm.getStruct(0)[4*0+1+0] == 2 );
  REQUIRE( l_comm.getStruct(0)[4*0+1+1] == 1 );
  REQUIRE( l_comm.getStruct(0)[4*0+1+2] == 0 );
  REQUIRE( l_comm.getStruct(0)[4*0+1+3] == 2 );

  REQUIRE( l_comm.getStruct(0)[4*1+1+0] == 2 );
  REQUIRE( l_comm.getStruct(0)[4*1+1+1] == 2 );
  REQUIRE( l_comm.getStruct(0)[4*1+1+2] == 1 );
  REQUIRE( l_comm.getStruct(0)[4*1+1+3] == 1 );

  REQUIRE( l_comm.getStruct(0)[4*2+1+0] == 2 );
  REQUIRE( l_comm.getStruct(0)[4*2+1+1] == 2 );
  REQUIRE( l_comm.getStruct(0)[4*2+1+2] == 3 );
  REQUIRE( l_comm.getStruct(0)[4*2+1+3] == 1 );

  REQUIRE( l_comm.getStruct(1)[4*0+0+0] == 2 );

  REQUIRE( l_comm.getStruct(1)[4*0+1+0] == 0 );
  REQUIRE( l_comm.getStruct(1)[4*0+1+1] == 0 );
  REQUIRE( l_comm.getStruct(1)[4*0+1+2] == 2 );
  REQUIRE( l_comm.getStruct(1)[4*0+1+3] == 2 );

  REQUIRE( l_comm.getStruct(1)[4*1+1+0] == 0 );
  REQUIRE( l_comm.getStruct(1)[4*1+1+1] == 2 );
  REQUIRE( l_comm.getStruct(1)[4*1+1+2] == 2 );
  REQUIRE( l_comm.getStruct(1)[4*1+1+3] == 1 );

  REQUIRE( l_comm.getStruct(2)[4*0+0+0] == 7 );

  REQUIRE( l_comm.getStruct(2)[4*0+1+0] == 1 );
  REQUIRE( l_comm.getStruct(2)[4*0+1+1] == 0 );
  REQUIRE( l_comm.getStruct(2)[4*0+1+2] == 2 );
  REQUIRE( l_comm.getStruct(2)[4*0+1+3] == 1 );

  REQUIRE( l_comm.getStruct(2)[4*1+1+0] == 1 );
  REQUIRE( l_comm.getStruct(2)[4*1+1+1] == 3 );
  REQUIRE( l_comm.getStruct(2)[4*1+1+2] == 1 );
  REQUIRE( l_comm.getStruct(2)[4*1+1+3] == 1 );

  REQUIRE( l_comm.getStruct(2)[4*2+1+0] == 2 );
  REQUIRE( l_comm.getStruct(2)[4*2+1+1] == 1 );
  REQUIRE( l_comm.getStruct(2)[4*2+1+2] == 0 );
  REQUIRE( l_comm.getStruct(2)[4*2+1+3] == 1 );

  REQUIRE( l_comm.getStruct(2)[4*3+1+0] == 2 );
  REQUIRE( l_comm.getStruct(2)[4*3+1+1] == 5 );
  REQUIRE( l_comm.getStruct(2)[4*3+1+2] == 2 );
  REQUIRE( l_comm.getStruct(2)[4*3+1+3] == 2 );

  // [...]

  REQUIRE( l_comm.m_sendRecvOff[0] == 0 );
  REQUIRE( l_comm.m_sendRecvOff[1] == 4 );
  REQUIRE( l_comm.m_sendRecvOff[2] == 7 );
  REQUIRE( l_comm.m_sendRecvOff[3] == 15 );
  REQUIRE( l_comm.m_sendRecvOff[4] == 16 );
  REQUIRE( l_comm.m_sendRecvOff[5] == 17 );
  REQUIRE( l_comm.m_sendRecvOff[6] == 22 );

  REQUIRE( l_comm.getSendEl(0)[0] == 0 - 0);
  REQUIRE( l_comm.getSendFa(0)[0] == 1    );

  REQUIRE( l_comm.getSendEl(0)[1] == 1 - 0);
  REQUIRE( l_comm.getSendFa(0)[1] == 1    );

  REQUIRE( l_comm.getSendEl(0)[2] == 1 - 0);
  REQUIRE( l_comm.getSendFa(0)[2] == 2    );

  REQUIRE( l_comm.getSendEl(0)[3] == 0 - 0);
  REQUIRE( l_comm.getSendFa(0)[3] == 0    );


  REQUIRE( l_comm.getSendEl(1)[0] == 2 - 2);
  REQUIRE( l_comm.getSendFa(1)[0] == 0    );

  REQUIRE( l_comm.getSendEl(1)[1] == 2 - 2);
  REQUIRE( l_comm.getSendFa(1)[1] == 1    );

  REQUIRE( l_comm.getSendEl(1)[2] == 2 - 2);
  REQUIRE( l_comm.getSendFa(1)[2] == 2    );


  REQUIRE( l_comm.getSendEl(2)[0] == 3 - 3);
  REQUIRE( l_comm.getSendFa(2)[0] == 1    );

  REQUIRE( l_comm.getSendEl(2)[1] == 3 - 3);
  REQUIRE( l_comm.getSendFa(2)[1] == 0    );

  REQUIRE( l_comm.getSendEl(2)[2] == 6 - 3);
  REQUIRE( l_comm.getSendFa(2)[2] == 1    );

  REQUIRE( l_comm.getSendEl(2)[3] == 4 - 3);
  REQUIRE( l_comm.getSendFa(2)[3] == 0    );

  REQUIRE( l_comm.getSendEl(2)[4] == 5 - 3);
  REQUIRE( l_comm.getSendFa(2)[4] == 1    );

  REQUIRE( l_comm.getSendEl(2)[5] == 6 - 3);
  REQUIRE( l_comm.getSendFa(2)[5] == 0    );

  REQUIRE( l_comm.getSendEl(2)[6] == 7 - 3);
  REQUIRE( l_comm.getSendFa(2)[6] == 2    );

  REQUIRE( l_comm.getSendEl(3)[0] == 8 - 8);
  REQUIRE( l_comm.getSendFa(3)[0] == 0    );

  // [...]

  REQUIRE( l_comm.getRecvEl(0)[0] == 0 - 0 );
  REQUIRE( l_comm.getRecvFa(0)[0] == 1     );

  REQUIRE( l_comm.getRecvEl(0)[1] == 1 - 0 );
  REQUIRE( l_comm.getRecvFa(0)[1] == 1     );

  REQUIRE( l_comm.getRecvEl(0)[2] == 1 - 0 );
  REQUIRE( l_comm.getRecvFa(0)[2] == 2     );

  REQUIRE( l_comm.getRecvEl(0)[3] == 0 - 0 );
  REQUIRE( l_comm.getRecvFa(0)[3] == 0     );


  REQUIRE( l_comm.getRecvEl(1)[0] == 2 - 2 );
  REQUIRE( l_comm.getRecvFa(1)[0] == 0     );

  REQUIRE( l_comm.getRecvEl(1)[1] == 2 - 2 );
  REQUIRE( l_comm.getRecvFa(1)[1] == 1     );

  REQUIRE( l_comm.getRecvEl(1)[2] == 2 - 2 );
  REQUIRE( l_comm.getRecvFa(1)[2] == 2     );


  REQUIRE( l_comm.getRecvEl(2)[0] == 3 - 3);
  REQUIRE( l_comm.getRecvFa(2)[0] == 1 );

  REQUIRE( l_comm.getRecvEl(2)[1] == 3 - 3);
  REQUIRE( l_comm.getRecvFa(2)[1] == 0 );

  REQUIRE( l_comm.getRecvEl(2)[2] == 6 - 3);
  REQUIRE( l_comm.getRecvFa(2)[2] == 1 );

  REQUIRE( l_comm.getRecvEl(2)[3] == 5 - 3);
  REQUIRE( l_comm.getRecvFa(2)[3] == 1 );

  REQUIRE( l_comm.getRecvEl(2)[4] == 4 - 3);
  REQUIRE( l_comm.getRecvFa(2)[4] == 0 );

  REQUIRE( l_comm.getRecvEl(2)[5] == 6 - 3);
  REQUIRE( l_comm.getRecvFa(2)[5] == 0 );

  // [...]
}
/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Unit tests of the distributed memory implementation.
 **/
#include <catch.hpp>

#define private public
#include "Mpi.h"
#undef private

TEST_CASE( "Initialization of the communication layout", "[initLayout]" ) {
  edge::parallel::Mpi l_mpi;
#ifdef USE_MPI
  /**
   * Setup example entity layout:
   *
   *        inner    |       send                                | receive
   * tg0 [  0 - 10 ] | 0:[ 11 - 21 ] 3:[ 22 - 41 ] 4:[ 42 - 91 ] | 0:[ 92 - 103] 3:[104 - 109] 4:[110-112]
   * tg1 [113 - 117] | 1:[118 - 119] 3:[120 - 158]               | 1:[159 - 171] 3:[172 - 188]
   * tg2 [189 - 920]
   *
   **/
  t_enLayout l_enLa;

  // resize vectors
  l_enLa.nEnts = 921;
  l_enLa.timeGroups.resize( 3 );

  l_enLa.timeGroups[0].send.resize(    3 );
  l_enLa.timeGroups[0].receive.resize( 3 );
  l_enLa.timeGroups[0].neRanks.resize( 3 );

  l_enLa.timeGroups[1].send.resize(    2 );
  l_enLa.timeGroups[1].receive.resize( 2 );
  l_enLa.timeGroups[1].neRanks.resize( 2 );

  // assign regions
  l_enLa.timeGroups[0].nEntsOwn         =  91;
  l_enLa.timeGroups[0].nEntsNotOwn      =  21;
  l_enLa.timeGroups[0].inner.first      =   0;
  l_enLa.timeGroups[0].inner.size       =  11;

  l_enLa.timeGroups[0].send[0].first    =  11;
  l_enLa.timeGroups[0].send[0].size     =  11;
  l_enLa.timeGroups[0].send[1].first    =  22;
  l_enLa.timeGroups[0].send[1].size     =  20;
  l_enLa.timeGroups[0].send[2].first    =  42;
  l_enLa.timeGroups[0].send[2].size     =  50;

  l_enLa.timeGroups[0].receive[0].first =  92;
  l_enLa.timeGroups[0].receive[0].size  =  12;
  l_enLa.timeGroups[0].receive[1].first = 104;
  l_enLa.timeGroups[0].receive[1].size  =   6;
  l_enLa.timeGroups[0].receive[2].first = 110;
  l_enLa.timeGroups[0].receive[2].size  =   3;

  l_enLa.timeGroups[0].neRanks[0] = 0;
  l_enLa.timeGroups[0].neRanks[1] = 3;
  l_enLa.timeGroups[0].neRanks[2] = 4;

  l_enLa.timeGroups[1].nEntsOwn         =  46;
  l_enLa.timeGroups[1].nEntsNotOwn      =  30;
  l_enLa.timeGroups[1].inner.first      = 113;
  l_enLa.timeGroups[1].inner.size       =   5;

  l_enLa.timeGroups[1].send[0].first    = 118;
  l_enLa.timeGroups[1].send[0].size     =   2;
  l_enLa.timeGroups[1].send[1].first    = 120;
  l_enLa.timeGroups[1].send[1].size     =  39;

  l_enLa.timeGroups[1].receive[0].first = 159;
  l_enLa.timeGroups[1].receive[0].size  =  13;
  l_enLa.timeGroups[1].receive[1].first = 172;
  l_enLa.timeGroups[1].receive[1].size  =  17;

  l_enLa.timeGroups[1].neRanks[0] = 1;
  l_enLa.timeGroups[1].neRanks[1] = 3;

  l_enLa.timeGroups[2].nEntsOwn         = 732;
  l_enLa.timeGroups[2].nEntsNotOwn      =   0;
  l_enLa.timeGroups[2].inner.first      = 189;
  l_enLa.timeGroups[2].inner.size       = 732;


  // create a dummy pointer
  double l_data;
  double *l_dPtr = &l_data;

  // create dummy size per entry
  std::size_t l_bytes = 20 * 9 * 8 * sizeof(double);

  edge::parallel::g_rank = 2;
  edge::parallel::g_nRanks = 21;

  l_mpi.initLayout( l_enLa, l_dPtr, l_bytes, 2, 5 );

  // check for sizes and ne ranks
  REQUIRE( l_mpi.m_send.size() == 3 );
  REQUIRE( l_mpi.m_recv.size() == 3 );

  for( int_tg l_tg = 0; l_tg < 3; l_tg++ ) {
    REQUIRE( l_mpi.m_send[l_tg].size() == l_enLa.timeGroups[l_tg].neRanks.size() );
    REQUIRE( l_mpi.m_recv[l_tg].size() == l_enLa.timeGroups[l_tg].neRanks.size() );

    for( std::size_t l_ne = 0 ; l_ne < l_mpi.m_send[l_tg].size(); l_ne++ ) {
      REQUIRE( l_mpi.m_send[l_tg][l_ne].rank == l_enLa.timeGroups[l_tg].neRanks[l_ne] );
      REQUIRE( l_mpi.m_recv[l_tg][l_ne].rank == l_enLa.timeGroups[l_tg].neRanks[l_ne] );
    }
  }

  // check tags
// TODO: Refactor the tests for the tags according to structure negleting MPI-info.
  REQUIRE( l_mpi.m_send[0][0].tag ==   2     * (5 * 21)    // own rank
                                     + (2+0) *  21         // time group
                                     + 0                ); // neighboring rank

  REQUIRE( l_mpi.m_send[0][1].tag ==   2     * (5 * 21)    // own rank
                                     + (2+0) *  21         // time group
                                     + 3                ); // neighboring rank

  REQUIRE( l_mpi.m_send[0][2].tag ==   2     * (5 * 21)    // own rank
                                     + (2+0) *  21         // time group
                                     + 4                ); // neighboring rank

  REQUIRE( l_mpi.m_recv[0][0].tag ==   0     * (5 * 21)    // neighboring rank
                                     + (2+0) *  21         // time group
                                     + 2                ); // own rank

  REQUIRE( l_mpi.m_recv[0][1].tag ==   3     * (5 * 21)    // neighboring rank
                                     + (2+0) *  21         // time group
                                     + 2                ); // own rank

  REQUIRE( l_mpi.m_recv[0][2].tag ==   4     * (5 * 21)    // neighboring rank
                                     + (2+0) *  21         // time group
                                     + 2                ); // own rank


  REQUIRE( l_mpi.m_send[1][0].tag ==   2     * (5 * 21)    // own rank
                                     + (2+1) *  21         // time group
                                     + 1                ); // neighboring rank

  REQUIRE( l_mpi.m_send[1][1].tag ==   2     * (5 * 21)    // own rank
                                     + (2+1) *  21         // time group
                                     + 3                ); // neighboring rank

  REQUIRE( l_mpi.m_recv[1][0].tag ==   1     * (5 * 21)    // neighboring rank
                                     + (2+1) *  21         // time group
                                     + 2                ); // own rank

  REQUIRE( l_mpi.m_recv[1][1].tag ==   3     * (5 * 21)    // neighboring rank
                                     + (2+1) *  21         // time group
                                     + 2                ); // own rank

  // check the ptrs and message size
  REQUIRE( l_mpi.m_send[0][0].ptr  == l_dPtr +  11 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_send[0][0].size ==           11 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_send[0][1].ptr  == l_dPtr +  22 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_send[0][1].size ==           20 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_send[0][2].ptr  == l_dPtr +  42 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_send[0][2].size ==           50 * (20 * 9 * 8) * 8 );


  REQUIRE( l_mpi.m_recv[0][0].ptr  == l_dPtr +  92 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_recv[0][0].size ==           12 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_recv[0][1].ptr  == l_dPtr + 104 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_recv[0][1].size ==            6 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_recv[0][2].ptr  == l_dPtr + 110 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_recv[0][2].size ==            3 * (20 * 9 * 8) * 8 );


  REQUIRE( l_mpi.m_send[1][0].ptr  == l_dPtr + 118 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_send[1][0].size ==            2 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_send[1][1].ptr  == l_dPtr + 120 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_send[1][1].size ==           39 * (20 * 9 * 8) * 8 );


  REQUIRE( l_mpi.m_recv[1][0].ptr  == l_dPtr + 159 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_recv[1][0].size ==           13 * (20 * 9 * 8) * 8 );

  REQUIRE( l_mpi.m_recv[1][1].ptr  == l_dPtr + 172 * (20 * 9 * 8)     );
  REQUIRE( l_mpi.m_recv[1][1].size ==           17 * (20 * 9 * 8) * 8 );

  // check number of messages
  REQUIRE( l_mpi.m_nMsgs == 10 );
#endif
}

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
 * Tests regular discretization into tets.
 **/

#include <catch.hpp>

#include "io/logging.h"

#define private public
#include "Tet.h"
#undef private

// test-fixture which sets MPI-info manually
class RegTet{
  private:
    int m_rank;
    int m_nRanks;
  public:
    RegTet(){
      // store current info
      m_rank   = edge::parallel::g_rank;
      m_nRanks = edge::parallel::g_nRanks;

      edge::parallel::g_rank = 0;
      edge::parallel::g_nRanks = 4096;
    };
    ~RegTet() {
      // restore MPI-info
      edge::parallel::g_rank   = m_rank;
      edge::parallel::g_nRanks = m_nRanks;
    }
};

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Number of entities", "[regTet][nEn]" ) {
  edge::mesh::regular::Tet l_tet;
  double l_dummy[3];

  // number of elements
  int_el l_nEl;

  unsigned int l_nHex[3] = {13,9,3};
  int l_mpiNe[3][2] = { {0,0}, {0,0}, {0,0} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  l_nEl = l_tet.getNEl();

  REQUIRE( l_nEl == 13*9*3*5 );
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Tets in hexes", "[regTet][inHex]" ) {
  edge::mesh::regular::Tet l_tet;
  double l_dummy[3];

  // detection of receive elements
  bool l_recv;

  unsigned int l_nHex[3] = {15,17,21};
  int l_mpiNe[3][2] = { {1,2}, {3,4}, {5,6} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_tet.m_periodic = true;

  int_el l_id[3] = {13, 2, 14};

  l_recv = l_tet.isRecvEl( 0, l_id );
  REQUIRE( l_recv == false );

  l_recv = l_tet.isRecvEl( 1, l_id );
  REQUIRE( l_recv == false );

  l_id[0] = -1; l_id[1] = 2; l_id[2] = 14;

  l_recv = l_tet.isRecvEl( 0, l_id );
  REQUIRE( l_recv == false );

  l_recv = l_tet.isRecvEl( 1, l_id );
  REQUIRE( l_recv == true );

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_recv = l_tet.isRecvEl( 0, l_id );
  REQUIRE( l_recv == true );

  l_id[0] = -1; l_id[1] = -1; l_id[2] = 14;

  l_recv = l_tet.isRecvEl( 1, l_id );
  REQUIRE( l_recv == false );

  l_id[0] = 0; l_id[1] = 0; l_id[2] = -1;

  l_recv = l_tet.isRecvEl( 0, l_id );
  REQUIRE( l_recv == false );

  l_recv = l_tet.isRecvEl( 3, l_id );
  REQUIRE( l_recv == true );

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_recv = l_tet.isRecvEl( 3, l_id );
  REQUIRE( l_recv == true );

  // detection of send elements
  bool l_send[3];

  l_id[0] = 4; l_id[1] = 5; l_id[2] = 10;

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == false );

  l_id[0] = 0; l_id[1] = 5; l_id[2] = 10;

  l_tet.isSendEl( 0, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == false );

  l_id[0] = 0; l_id[1] = 0; l_id[2] = 10;

  l_tet.isSendEl( 1, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == false );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == false );

  l_id[0] = 0; l_id[1] = 0; l_id[2] = 20;

  l_tet.isSendEl( 0, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == false );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == true );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == true );

  l_tet.isSendEl( 2, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == false );

  l_tet.isSendEl( 0, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == false );

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == true );

  l_tet.isSendEl( 4, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == true );

  l_mpiNe[0][0] = -1; l_mpiNe[0][1] = 2;
  l_mpiNe[1][0] = -1; l_mpiNe[1][1] = 4;
  l_mpiNe[2][0] = -1; l_mpiNe[2][1] = 6;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == false );
  REQUIRE( l_send[2] == true );

  l_mpiNe[0][0] = -1; l_mpiNe[0][1] =  2;
  l_mpiNe[1][0] =  3; l_mpiNe[1][1] =  4;
  l_mpiNe[2][0] = -1; l_mpiNe[2][1] = -1;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == false );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == false );

  l_mpiNe[0][0] =  1; l_mpiNe[0][1] = -1;
  l_mpiNe[1][0] =  3; l_mpiNe[1][1] = -1;
  l_mpiNe[2][0] = -1; l_mpiNe[2][1] = -1;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_tet.isSendEl( 3, l_id, l_send );
  REQUIRE( l_send[0] == true );
  REQUIRE( l_send[1] == true );
  REQUIRE( l_send[2] == false );

  // detection of inner elements
  bool l_inner;

  l_mpiNe[0][0] =  1; l_mpiNe[0][1] = -1;
  l_mpiNe[1][0] =  3; l_mpiNe[1][1] = -1;
  l_mpiNe[2][0] = -1; l_mpiNe[2][1] = -1;

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  l_id[0] = 10; l_id[1] = 14; l_id[2] = 9;

  l_inner = l_tet.isInnerEl( 3, l_id );
  REQUIRE( l_inner == true );

  l_mpiNe[0][0] = -1; l_mpiNe[0][1] = 2;
  l_mpiNe[1][0] = -1; l_mpiNe[1][1] = 4;
  l_mpiNe[2][0] = -1; l_mpiNe[2][1] = 6;

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  l_id[0] = 0; l_id[1] = 0; l_id[2] = 0;

  l_inner = l_tet.isInnerEl( 3, l_id );
  REQUIRE( l_inner == true );

  l_id[0] = 0; l_id[1] = 0; l_id[2] = 0;

  l_mpiNe[0][0] =  1; l_mpiNe[0][1] = 2;
  l_mpiNe[1][0] =  3; l_mpiNe[1][1] = 4;
  l_mpiNe[2][0] =  5; l_mpiNe[2][1] = 6;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  l_inner = l_tet.isInnerEl( 4, l_id );
  REQUIRE( l_inner == true );

  l_inner = l_tet.isInnerEl( 2, l_id );
  REQUIRE( l_inner == true );

  l_inner = l_tet.isInnerEl( 3, l_id );
  REQUIRE( l_inner == false );

  l_inner = l_tet.isInnerEl( 1, l_id );
  REQUIRE( l_inner == false );
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Mesh Structure", "[regTet][struct]" ) {
  edge::mesh::regular::Tet l_tet;
  double l_dummy[3];

  unsigned int l_nHex[3] = {3,5,2};
  int l_mpiNe[3][2] = { {0,0}, {0,0}, {0,0} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size     == 150 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send.size()    == 0 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive.size() == 0 );

  l_mpiNe[0][0] = 3; l_mpiNe[0][1] = 0;
  l_mpiNe[1][0] = 0; l_mpiNe[1][1] = 0;
  l_mpiNe[2][0] = 0; l_mpiNe[2][1] = 0;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size      == 150 - (5*2*2) );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size    == 5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size == 5*2*2 );

  l_mpiNe[0][0] = 3; l_mpiNe[0][1] = 3;
  l_mpiNe[1][0] = 0; l_mpiNe[1][1] = 0;
  l_mpiNe[2][0] = 0; l_mpiNe[2][1] = 0;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size      == 150 - (5*2*2*2) );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size    == 5*2*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send.size()     == 1       );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size == 5*2*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive.size()  == 1 );

  l_mpiNe[0][0] = 0; l_mpiNe[0][1] = 0;
  l_mpiNe[1][0] = 0; l_mpiNe[1][1] = 0;
  l_mpiNe[2][0] = 1; l_mpiNe[2][1] = 1;

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size      == 150 - (3*5*2*2) );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size    == 3*5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size == 3*5*2*2 );

  l_mpiNe[0][0] = 0; l_mpiNe[0][1] = 0;
  l_mpiNe[1][0] = 0; l_mpiNe[1][1] = 0;
  l_mpiNe[2][0] = 1; l_mpiNe[2][1] = 2;

  l_tet.init( true,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size      == 150 - (3*5*2) - (3*5*2) );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size    == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[1].size    == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[1].size == 3*5*2 );

  l_mpiNe[0][0] = 1; l_mpiNe[0][1] = 4;
  l_mpiNe[1][0] = 2; l_mpiNe[1][1] = 5;
  l_mpiNe[2][0] = 3; l_mpiNe[2][1] = 6;

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  int_el l_nIn = l_tet.getNElInner();
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size == l_nIn );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size    == 5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[1].size    == 5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[2].size    == 3*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[3].size    == 3*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[4].size    == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[5].size    == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size == 5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[1].size == 5*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[2].size == 3*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[3].size == 3*2*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[4].size == 3*5*2 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[5].size == 3*5*2 );
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Adjacencies", "[regTet][adj]" ) {
  edge::mesh::regular::Tet l_tet;
  double l_dummy[3];

  unsigned int l_nHex[3] = {17,7,23};
  int l_mpiNe[3][2] = { {0,0}, {0,0}, {0,0} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );
  l_tet.m_periodic=false;

  /*
   * hex: x=15, y=5, z=13 has id 14*19*9 + 6*19 + 16 = 2524.
   * vertex id 0 of this hex is 14*20*10 + 6*20*16 + 16 = 2936
   * tet 0  in the hex is given by: (13*7*17 + 5*17 + 15)*5 = 8235 
   * the hex type of 0 is B)
   */
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size == 17*7*23*5 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks.size() == 0 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send.size() == 0 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive.size() == 0 );

  REQUIRE( l_tet.m_hexes[2524].el[0][0] == 8235       );
  REQUIRE( l_tet.m_hexes[2524].ve[0]    == 2936       );
  REQUIRE( l_tet.m_hexes[2524].ve[1]    == 2937       );
  REQUIRE( l_tet.m_hexes[2524].ve[3]    == 2936+20    );
  REQUIRE( l_tet.m_hexes[2524].ve[4]    == 2936+20*10 );

  // get the real data
  int_el l_nEl = l_tet.getNEl();
  int_el (*l_elVe)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];

  l_tet.getElVe( &l_elVe[0] );

  // check for the same data in the derived data structure
  REQUIRE( l_elVe[8235][0] == 2936       );
  REQUIRE( l_elVe[8235][1] == 2937       );
  REQUIRE( l_elVe[8235][2] == 2936+20    );
  REQUIRE( l_elVe[8235][3] == 2936+20*10 );

  delete[] l_elVe;

  l_mpiNe[0][0] = 1; l_mpiNe[0][1] = 4;
  l_mpiNe[1][0] = 2; l_mpiNe[1][1] = 5;
  l_mpiNe[2][0] = 3; l_mpiNe[2][1] = 6;

  // check faces in an MPI-setup
  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  REQUIRE( l_tet.m_faLayout.nEnts == 17*7*23*10 +   // default
                                     7 *23*2    +   // left bnd
                                     17*23*2    +   // front bd
                                     17* 7*2     ); // bottom bnd
  REQUIRE( l_tet.m_faLayout.timeGroups[0].nEntsOwn         == 17*7*23*10 );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[0].size     == 0          );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[1].size     == 7*23*2     );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[2].size     == 0          );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[3].size     == 17*23*2    );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[4].size     == 0          );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].send[5].size     == 17*7*2     );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[0].size  == 7*23*2     );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[1].size  == 0          );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[2].size  == 17*23*2    );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[3].size  == 0          );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[4].size  == 17*7*2     );
  REQUIRE( l_tet.m_faLayout.timeGroups[0].receive[5].size  == 0          );

  // hex: x=15, y=5, z=13 has id 14*19*9 + 6*19 + 16 = 2524.
  // all faces of the hex are inner faces, the respective ID are relative to inner faces only
  // 13*7*17+5*17+15 = 1,647 owned hexes appear before this hex; every hex adds 10 faces,
  // where the bnd-faces are all MPI and thus appear later.
  // -> the first face of our hex would be 1,647*10=16,470 if it weren't for the right and back MPI-bnd.
  // -> right MPI: 13*7*2  +  5*2 = 192
  // -> back  MPI: 13*17*2        = 442
  // -> 16,470 - 192 - 442 = 15836
  REQUIRE( l_tet.m_hexes[2524].fa[2]  == 15836 ); // 1st newly created face
  REQUIRE( l_tet.m_hexes[2524].fa[3]  == 15837 ); // 2nd
  REQUIRE( l_tet.m_hexes[2524].fa[6]  == 15838 ); // 3rd
  REQUIRE( l_tet.m_hexes[2524].fa[7]  == 15839 ); // 4th
  REQUIRE( l_tet.m_hexes[2524].fa[10] == 15840 ); // 5th
  REQUIRE( l_tet.m_hexes[2524].fa[11] == 15841 ); // 6th
  REQUIRE( l_tet.m_hexes[2524].fa[12] == 15842 ); // 7th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[13] == 15843 ); // 8th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[14] == 15844 ); // 9th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[15] == 15845 ); // 10th (inner)

  // left hex id  : 2524-1    = 2523 
  // front hex id:  2524-19   = 2505
  // bottom hex id: 2524-19*9 = 2353
  REQUIRE( l_tet.m_hexes[2524].fa[0]  == l_tet.m_hexes[2523].fa[2] );  // copied from left element
  REQUIRE( l_tet.m_hexes[2524].fa[1]  == l_tet.m_hexes[2523].fa[3] );  // copied from left element
  REQUIRE( l_tet.m_hexes[2524].fa[4]  == l_tet.m_hexes[2505].fa[6] );  // copied from front element
  REQUIRE( l_tet.m_hexes[2524].fa[5]  == l_tet.m_hexes[2505].fa[7] );  // copied from front element
  REQUIRE( l_tet.m_hexes[2524].fa[8]  == l_tet.m_hexes[2353].fa[10] ); // copied from bottom element
  REQUIRE( l_tet.m_hexes[2524].fa[9]  == l_tet.m_hexes[2353].fa[11] ); // copied from bottom element

  // check through face's vertices
  int_el l_nFa = l_tet.getNFa();
  int_el (*l_faVe)[3] = ( int_el(*)[3] ) new int_el[3 * l_nFa];

  l_tet.getFaVe( l_faVe );

  // face type is B)
  REQUIRE( l_faVe[15836][0] == l_tet.m_hexes[2524].ve[1] );
  REQUIRE( l_faVe[15836][1] == l_tet.m_hexes[2524].ve[2] );
  REQUIRE( l_faVe[15836][2] == l_tet.m_hexes[2524].ve[6] );

  REQUIRE( l_faVe[15844][0] == l_tet.m_hexes[2524].ve[1] );
  REQUIRE( l_faVe[15844][1] == l_tet.m_hexes[2524].ve[4] );
  REQUIRE( l_faVe[15844][2] == l_tet.m_hexes[2524].ve[6] );

  REQUIRE( l_faVe[15845][0] == l_tet.m_hexes[2524].ve[3] );
  REQUIRE( l_faVe[15845][1] == l_tet.m_hexes[2524].ve[4] );
  REQUIRE( l_faVe[15845][2] == l_tet.m_hexes[2524].ve[6] );

  delete[] l_faVe;
 
  /*
   * check elements adjacent to faces
   *
   * now that send-elements have higher ids than inner-elements, we have to remove them
   * from the previous counting wich reached 8235 elements until hex x=15, y=5, z=13.
   *
   *
   * #inner-elements in bottom layer:
   *
   * (16x2*2 + 5x2*2 + 15x5x3 + 2) = 311
   * 
   *    type A) corner
   *      |
   *     \|/
   *     *********************...*************************
   *     *1    *             16x2                        *
   *     *     *                                         *
   *     *********************...*************************
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     .  5  .                                   .  5  .
   *     .  x  .            15x5x3                 .  x  .
   *     .  2  .                                   .  2  .
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     * 
   *     *********************...*************************
   *     *                                         *     *
   *     *                   16x2                  *    1* <- type A) corner
   *     *********************...*************************
   *
   *
   * #innner-elements in layers 1-13:
   *
   * (2*4 + 15x3*2 + 5x3*2 + 15x5x5)*12 = 6036
   *
   *     *********************...*************************
   *     *2    *             15x3                  *    2*
   *     *     *                                   *     *
   *     *********************...*************************
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     .  5  .                                   .  5  .
   *     .  x  .            15x5x5                 .  x  .
   *     .  3  .                                   .  3  .
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *     *                                   *     * 
   *     *********************...*************************
   *     *     *                                   *     *
   *     *2    *             15x3                  *    2*
   *     *********************...*************************
   *
   * #inner-elements in the current layer:
   *
   * 2*2 + 15x3 + 5x3 + 4x3 + 15x4x5 + 14x5  = 446
   *
   *     *********************...*************************
   *     *     *                                   *     *
   *     *     *                                   *     *
   *     *********************...*************************
   *     *     *       14x5                  *     *     *
   *     *     *                             *     *     *
   *     *     * *************...*************************
   *     *     *                               /|\ *     *
   *     .  5  .                                |  .  4  .
   *     .  x  .                          hex x=15 .  x  .
   *     .  3  .                              y= 5 .  3  .
   *     *     *                              z=13 *     *
   *     *     *                                   *     *
   *     *     *          15x4x5                   *     *
   *     *     *                                   *     * 
   *     *********************...*************************
   *     *     *                                   *     *
   *     *2    *             15x3                  *    2*
   *     *********************...*************************
   *
   * -> tet0 of hex 2524 has id 6793.
   */

  // first inner element in hex 0, 0, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*1+ 1].el[2][0] == 0    );

  // first inner element in hex 1, 0, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*1+ 2].el[2][0] == 2    );
  // second inner element in hex 1, 0, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*1+ 2].el[4][0] == 3    );

  // first inner element in hex 16, 0, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*1+17].el[2][0] == 16*2   );

  // first inner element in 0, 1, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*2+ 1].el[2][0] == 16*2+1 );

  // first inner element in hex 15, 1, 0
  REQUIRE( l_tet.m_hexes[(19*9)*1+19*2+16].el[2][0] == 16*2+1 + 2 + 14*3 );
  // first inner element in hex 0,0,1
  REQUIRE( l_tet.m_hexes[(19*9)*2+19*1+ 1].el[1][0] == 311 );


  // first inner element of hex 15, 5, 13
  REQUIRE( l_tet.m_hexes[2524].el[0][0] == 6793 );

  l_nFa = l_tet.getNFa();
  int_el (*l_faEl)[2] = ( int_el(*)[2] ) new int_el[2 * l_nFa];

  l_tet.getFaEl( l_faEl );

  // 15, 5, 13 is a type B) element; 16, 5, 13 is an MPI hex
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 0] ][0] == 6793-5   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 0] ][1] == 6793     );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 1] ][0] == 6793-1   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 1] ][1] == 6793+4   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 2] ][0] == 6793+1   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 2] ][1] == 6793+5   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 3] ][0] == 6793+3   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 3] ][1] == 6793+5+2 );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 4] ][1] == 6793     );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 5] ][1] == 6793+3   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 6] ][0] == 6793+1   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 7] ][0] == 6793+4   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 8] ][1] == 6793+0   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[ 9] ][1] == 6793+1   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[10] ][0] == 6793+3   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[11] ][0] == 6793+4   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[12] ][0] == 6793+0   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[12] ][1] == 6793+2   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[13] ][0] == 6793+1   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[13] ][1] == 6793+2   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[14] ][0] == 6793+2   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[14] ][1] == 6793+3   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[15] ][0] == 6793+2   );
  REQUIRE( l_faEl[ l_tet.m_hexes[2524].fa[15] ][1] == 6793+4   );

  delete[] l_faEl;

  // check faces adjacent to elements
  l_nEl = l_tet.getNEl();
  int_el (*l_elFa)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];

  l_tet.getElFa( l_elFa );

  REQUIRE( l_tet.m_hexes[2524].fa[2]  == 15836 ); // 1st newly created face
  REQUIRE( l_tet.m_hexes[2524].fa[3]  == 15837 ); // 2nd
  REQUIRE( l_tet.m_hexes[2524].fa[6]  == 15838 ); // 3rd
  REQUIRE( l_tet.m_hexes[2524].fa[7]  == 15839 ); // 4th
  REQUIRE( l_tet.m_hexes[2524].fa[10] == 15840 ); // 5th
  REQUIRE( l_tet.m_hexes[2524].fa[11] == 15841 ); // 6th
  REQUIRE( l_tet.m_hexes[2524].fa[12] == 15842 ); // 7th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[13] == 15843 ); // 8th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[14] == 15844 ); // 9th (inner)
  REQUIRE( l_tet.m_hexes[2524].fa[15] == 15845 ); // 10th (inner)

  // left hex id  : 2524-1    = 2523 
  // front hex id:  2524-19   = 2505
  // bottom hex id: 2524-19*9 = 2353
  REQUIRE( l_tet.m_hexes[2524].fa[0]  == l_tet.m_hexes[2523].fa[2] );  // copied from left element
  REQUIRE( l_tet.m_hexes[2524].fa[1]  == l_tet.m_hexes[2523].fa[3] );  // copied from left element
  REQUIRE( l_tet.m_hexes[2524].fa[4]  == l_tet.m_hexes[2505].fa[6] );  // copied from front element
  REQUIRE( l_tet.m_hexes[2524].fa[5]  == l_tet.m_hexes[2505].fa[7] );  // copied from front element
  REQUIRE( l_tet.m_hexes[2524].fa[8]  == l_tet.m_hexes[2353].fa[10] ); // copied from bottom element
  REQUIRE( l_tet.m_hexes[2524].fa[9]  == l_tet.m_hexes[2353].fa[11] ); // copied from bottom element

  REQUIRE( l_elFa[6793][0] == l_tet.m_hexes[2353].fa[10] );
  REQUIRE( l_elFa[6793][1] == l_tet.m_hexes[2505].fa[ 6] );
  REQUIRE( l_elFa[6793][2] == l_tet.m_hexes[2523].fa[ 2] );
  REQUIRE( l_elFa[6793][3] == l_tet.m_hexes[2524].fa[12] );

  REQUIRE( l_elFa[6793+2][0] == l_tet.m_hexes[2524].fa[12] );
  REQUIRE( l_elFa[6793+2][1] == l_tet.m_hexes[2524].fa[13] );
  REQUIRE( l_elFa[6793+2][2] == l_tet.m_hexes[2524].fa[14] );
  REQUIRE( l_elFa[6793+2][3] == l_tet.m_hexes[2524].fa[15] );

  delete[] l_elFa;

  // now get the elements adjacent to elements through faces
  l_nEl = l_tet.getNEl();
  int_el (*l_elFaEl)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];

  l_tet.getElFaEl( l_elFaEl );

  // back to our favorite tet 6793
  REQUIRE( l_elFaEl[6793][0] == l_tet.m_hexes[2353].el[3][0] );
  REQUIRE( l_elFaEl[6793][1] == l_tet.m_hexes[2505].el[1][0] );
  REQUIRE( l_elFaEl[6793][2] == 6793-5                       );
  REQUIRE( l_elFaEl[6793][3] == 6793+2                       );

  delete[] l_elFaEl;

  // simple case with non-MPI periodic bnds
  l_mpiNe[0][0] = 0; l_mpiNe[0][1] = 0;
  l_mpiNe[1][0] = 0; l_mpiNe[1][1] = 0;
  l_mpiNe[2][0] = 0; l_mpiNe[2][1] = 0;

  l_nHex[0] = 16; l_nHex[1] = 32; l_nHex[2] = 6;

  l_tet.m_periodic = true;
  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  // get face adj
  l_nEl = l_tet.getNEl();
  l_elFa = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  l_tet.getElFa( l_elFa );

  // check that the periodic faces translate
  REQUIRE( l_elFa[1][1] == l_elFa[15*5+1   ][2] );
  REQUIRE( l_elFa[3][1] == l_elFa[15*5+3   ][2] );
  REQUIRE( l_elFa[0][1] == l_elFa[31*16*5+1][3] );
  REQUIRE( l_elFa[3][0] == l_elFa[31*16*5+4][2] );
  delete[] l_elFa;

  // another simple case with non-MPI periodic bnds
  l_nHex[0] = 8; l_nHex[1] = 8; l_nHex[2] = 8;

  l_tet.m_periodic = true;
  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_dummy, l_dummy );

  // get face adj
  l_nEl = l_tet.getNEl();
  l_elFa = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  l_tet.getElFa( l_elFa );

  l_nFa = l_tet.getNFa();
  l_faEl = ( int_el(*)[2] ) new int_el[2 * l_nFa];
  l_tet.getFaEl( l_faEl );

  // check that the periodic faces translate
  REQUIRE( l_elFa[1][1] == l_elFa[7*5+1  ][2] );
  REQUIRE( l_elFa[3][1] == l_elFa[7*5+3  ][2] );
  REQUIRE( l_elFa[0][1] == l_elFa[7*8*5+1][3] );
  REQUIRE( l_elFa[3][0] == l_elFa[7*8*5+4][2] );

  REQUIRE( l_elFa[0][1] == 4 );
  REQUIRE( l_faEl[4][0] == 0 );

  delete[] l_elFa;
  delete[] l_faEl;
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Vertex chars", "[regTet][veChars]" ) {
  edge::mesh::regular::Tet l_tet;

  // set up lower left coords and mesh width
  double l_crds[3] = { 3.5, 12.9, 952.3 };
  double l_dx[3]    = { 0.2,  0.5,   0.1 };

  unsigned int l_nHex[3] = {17,7,23};
  int l_mpiNe[3][2] = { {1,4}, {2,5}, {3,6} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_crds, l_dx );

  // check faces adjacent to elements
  int_el l_nVe = l_tet.getVeLayout().nEnts;

  t_vertexChars *l_veChars = new t_vertexChars[l_nVe];

  for( int_el l_ve = 0; l_ve < l_nVe; l_ve++ ) {
    l_veChars[l_ve].coords[0] = -999.9;
    l_veChars[l_ve].coords[1] = -999.9;
    l_veChars[l_ve].coords[2] = -999.9;
  }

  l_tet.getVeChars( l_veChars );

  for( int_el l_ve = 0; l_ve < l_nVe; l_ve++ ) {
    REQUIRE( l_veChars[l_ve].coords[0] != Approx(-999.9) );
    REQUIRE( l_veChars[l_ve].coords[1] != Approx(-999.9) );
    REQUIRE( l_veChars[l_ve].coords[2] != Approx(-999.9) );
  }


  // get vertices of elements
  int_el l_nEl = l_tet.getNEl();
  int_el (*l_elVe)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];

  l_tet.getElVe( &l_elVe[0] );

  // check vertices of tet0 (6793) in hex 15, 5, 13 (type B)
  REQUIRE( l_veChars[ l_elVe[6793][0] ].coords[0] == Approx(  3.5 + 0.2 * 15) );
  REQUIRE( l_veChars[ l_elVe[6793][0] ].coords[1] == Approx( 12.9 + 0.5 *  5) );
  REQUIRE( l_veChars[ l_elVe[6793][0] ].coords[2] == Approx(952.3 + 0.1 * 13) );

  REQUIRE( l_veChars[ l_elVe[6793][1] ].coords[0] == Approx(  3.5 + 0.2 * 16) );
  REQUIRE( l_veChars[ l_elVe[6793][1] ].coords[1] == Approx( 12.9 + 0.5 *  5) );
  REQUIRE( l_veChars[ l_elVe[6793][1] ].coords[2] == Approx(952.3 + 0.1 * 13) );

  REQUIRE( l_veChars[ l_elVe[6793][2] ].coords[0] == Approx(  3.5 + 0.2 * 15) );
  REQUIRE( l_veChars[ l_elVe[6793][2] ].coords[1] == Approx( 12.9 + 0.5 *  6) );
  REQUIRE( l_veChars[ l_elVe[6793][2] ].coords[2] == Approx(952.3 + 0.1 * 13) );

  REQUIRE( l_veChars[ l_elVe[6793][3] ].coords[0] == Approx(  3.5 + 0.2 * 15) );
  REQUIRE( l_veChars[ l_elVe[6793][3] ].coords[1] == Approx( 12.9 + 0.5 *  5) );
  REQUIRE( l_veChars[ l_elVe[6793][3] ].coords[2] == Approx(952.3 + 0.1 * 14) );

  delete[] l_veChars;
  delete[] l_elVe;
}

TEST_CASE_METHOD( RegTet,"Regular meshes with tets: Face chars", "[regTet][faChars]" ) {
  edge::mesh::regular::Tet l_tet;

  // set up lower left coords and mesh width
  double l_crds[3] = { 3.5, 12.9, 952.3 };
  double l_dx[3]    = { 0.2,  0.5,   0.1 };

  unsigned int l_nHex[3] = {17,7,23};
  int l_mpiNe[3][2] = { {1,4}, {2,5}, {3,6} };

  l_tet.m_periodic = false;
  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_crds, l_dx );

  // setup face chars
  int_el l_nFa = l_tet.getFaLayout().nEnts;
  REQUIRE( l_nFa == (17*7*23*10) + (17*7+17*23+7*23)*2 );

  t_faceChars *l_faChars = new t_faceChars[l_nFa];
  l_tet.getFaChars( l_faChars );

  // check newly created, coord aligned tet faces of hex 15, 5, 13 (type B)
  // we have 16,470 non-mpi faces by default, however the 'right' boundaries are not accounted for
  // removes (13*17 + 13*7 + 5)*2 = 634 faces for the boundaries
  //
  REQUIRE( l_faChars[15836].area == Approx(0.5 * 0.1 / 2.0 ) ); // 1st newly created face
  REQUIRE( l_faChars[15837].area == Approx(0.5 * 0.1 / 2.0 ) ); // 2nd
  REQUIRE( l_faChars[15838].area == Approx(0.2 * 0.1 / 2.0 ) ); // 3rd
  REQUIRE( l_faChars[15839].area == Approx(0.2 * 0.1 / 2.0 ) ); // 4th
  REQUIRE( l_faChars[15840].area == Approx(0.2 * 0.5 / 2.0 ) ); // 5th
  REQUIRE( l_faChars[15841].area == Approx(0.2 * 0.5 / 2.0 ) ); // 6th

  // check the out pointing normals
  l_dx[0] = 333; l_dx[1] = 333; l_dx[2] = 333;
  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_crds, l_dx );
  l_tet.getFaChars( l_faChars );

  REQUIRE( l_faChars[15836].outNormal[0] == Approx(1.0) );
  REQUIRE( l_faChars[15836].outNormal[1] == Approx(0.0) );
  REQUIRE( l_faChars[15836].outNormal[2] == Approx(0.0) );

  REQUIRE( l_faChars[15837].outNormal[0] == Approx(1.0) );
  REQUIRE( l_faChars[15837].outNormal[1] == Approx(0.0) );
  REQUIRE( l_faChars[15837].outNormal[2] == Approx(0.0) );

  REQUIRE( l_faChars[15838].outNormal[0] == Approx(0.0) );
  REQUIRE( l_faChars[15838].outNormal[1] == Approx(1.0) );
  REQUIRE( l_faChars[15838].outNormal[2] == Approx(0.0) );

  REQUIRE( l_faChars[15839].outNormal[0] == Approx(0.0) );
  REQUIRE( l_faChars[15839].outNormal[1] == Approx(1.0) );
  REQUIRE( l_faChars[15839].outNormal[2] == Approx(0.0) );

  REQUIRE( l_faChars[15840].outNormal[0] == Approx(0.0) );
  REQUIRE( l_faChars[15840].outNormal[1] == Approx(0.0) );
  REQUIRE( l_faChars[15840].outNormal[2] == Approx(1.0) );

  REQUIRE( l_faChars[15841].outNormal[0] == Approx(0.0) );
  REQUIRE( l_faChars[15841].outNormal[1] == Approx(0.0) );
  REQUIRE( l_faChars[15841].outNormal[2] == Approx(1.0) );

  // points in center el
  REQUIRE( l_faChars[15842].outNormal[0] == Approx(  1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15842].outNormal[1] == Approx(  1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15842].outNormal[2] == Approx(  1/std::sqrt(3) ) );

  // points in center el
  REQUIRE( l_faChars[15843].outNormal[0] == Approx( -1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15843].outNormal[1] == Approx( -1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15843].outNormal[2] == Approx(  1/std::sqrt(3) ) );

  // points out center el
  REQUIRE( l_faChars[15844].outNormal[0] == Approx(  1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15844].outNormal[1] == Approx( -1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15844].outNormal[2] == Approx(  1/std::sqrt(3) ) );

  // points out center el
  REQUIRE( l_faChars[15845].outNormal[0] == Approx( -1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15845].outNormal[1] == Approx(  1/std::sqrt(3) ) );
  REQUIRE( l_faChars[15845].outNormal[2] == Approx(  1/std::sqrt(3) ) );


  delete[] l_faChars;
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Element chars", "[regTet][elChars]" ) {
  edge::mesh::regular::Tet l_tet;

  // set up lower left coords and mesh width
  double l_crds[3] = { 3.5, 12.9, 952.3 };
  double l_dx[3]    = { 333,  333,   333 };

  unsigned int l_nHex[3] = {17,7,23};
  int l_mpiNe[3][2] = { {1,4}, {2,5}, {3,6} };

  l_tet.init( false,
              l_nHex,
              l_mpiNe,
              l_crds, l_dx );

  // get el chars
  int_el l_nEl = l_tet.getElLayout().nEnts;
  t_elementChars *l_elChars = new t_elementChars[l_nEl];

  l_tet.getElChars( l_elChars );

  REQUIRE( l_elChars[6793].volume == Approx(6154339.5)  );
  REQUIRE( l_elChars[6794].volume == Approx(6154339.5)  );
  REQUIRE( l_elChars[6794].volume == Approx(6154339.5)  );
  REQUIRE( l_elChars[6795].volume == Approx(12308679.0) );
  REQUIRE( l_elChars[6796].volume == Approx(6154339.5)  );
  REQUIRE( l_elChars[6797].volume == Approx(6154339.5)  );

  REQUIRE( l_elChars[6795].inDia  == Approx(96.1288198200644*2) );

  delete[] l_elChars;
}

TEST_CASE_METHOD( RegTet, "Regular meshes with tets: Init", "[regTet][init]" ) {
  // disable info logger
  el::Configurations l_conf;
  l_conf.parseFromText("* INFO:\n ENABLED = false");
  el::Loggers::reconfigureLogger("default", l_conf);

  // create tet-object
  edge::mesh::regular::Tet l_tet;

  // set up the mesh
  unsigned int l_nXHex[3] = {522, 632, 428};
  int l_rank = 721;
  int l_nRanks = 2048;
  double l_crds[3] = { -3.5,  4.25, 5.75 };
  double l_dx[3]    = { 0.5,  0.25, 0.75 };

  l_tet.init( l_nXHex, l_rank, l_nRanks, l_crds, l_dx, true );

  /*
   * 2048 -> 8 x 32 x 8 as partitioning
   *
   * 523 / 8 = 65.25
   * -> 0, 1 have 66 elements in x-direction
   * -> 2-7 have 65 elements
   *
   * 632 / 32 = 19.75
   * -> 0-23 have 20 elements in y-direction
   * -> 24-31 have 19 elements
   *
   * 428 / 8 = 53.5
   * -> 0-3 have 54 elements
   * -> 4-7 have 53 elements
   *
   */

  /*
   * our local rank is 721, this is partition 1,26,2
   * -> 66 x 19 x 54 is the local domain size.
   */
  REQUIRE( l_tet.m_nHex[0] == 66 );
  REQUIRE( l_tet.m_nHex[1] == 19 );
  REQUIRE( l_tet.m_nHex[2] == 54 );

  /*
   * the lowest left corner is -3.5, 4.25, 5.75
   * we have
   *   66 hexes in front of us in x-direction
   *   24*20 + 2*19=518 hexes in front of us in y-direction
   *   54*2=108 hex in front of us in z-direction
   */
  REQUIRE( l_tet.m_corner[0] == Approx( -3.5  +  66 * 0.5  ) );
  REQUIRE( l_tet.m_corner[1] == Approx(  4.25 + 518 * 0.25 ) );
  REQUIRE( l_tet.m_corner[2] == Approx(  5.75 + 108 * 0.75 ) );

  /*
   * Our MPI neighbors are
   *  720, 722 in x-direction
   *  713, 729 in y-direction
   *  465, 977 in z-direction
   */
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[0] == 720 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[1] == 722 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[2] == 713 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[3] == 729 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[4] == 465 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].neRanks[5] == 977 );

  // check the "messages"
  /**
   * we have 5*66*19*54=338,580 owned tets in the mesh (no duplication) and
   * 66*19*4 + 66*54*4 + 19*54*4 = 23,376 ghost tets
   *
   * for every "edge" of the domain, we share a tet, which gets duplicated:
   * 66*4+19*4+54*4= 556
   *
   * Since we have 19 hexes in y-direction, the bottom corner-types are mirrored at the top-bnd.
   * Therfore the 4 corners have each one tet which has 3 MPI-neighbors and therefore gets duplicate three times.
   */
  // inner elements
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.first == 0 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].inner.size  == 338580-23376+556-4 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].nEntsOwn    == 338580+556-4 );
  // send-elements
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].first == 338580-23376+556-4 );
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[0].size  == 19*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[1].first == l_tet.m_elLayout.timeGroups[0].send[0].first + 19*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[1].size  == 19*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[2].first == l_tet.m_elLayout.timeGroups[0].send[1].first + 19*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[2].size  == 66*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[3].first == l_tet.m_elLayout.timeGroups[0].send[2].first + 66*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[3].size  == 66*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[4].first == l_tet.m_elLayout.timeGroups[0].send[3].first + 66*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[4].size  == 66*19*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[5].first == l_tet.m_elLayout.timeGroups[0].send[4].first + 66*19*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].send[5].size  == 66*19*2 );
  // receive-elements
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].first == l_tet.m_elLayout.timeGroups[0].send[5].first + 66*19*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[0].size  == 19*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[1].first == l_tet.m_elLayout.timeGroups[0].receive[0].first + 19*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[1].size  == 19*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[2].first == l_tet.m_elLayout.timeGroups[0].receive[1].first + 19*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[2].size  == 66*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[3].first == l_tet.m_elLayout.timeGroups[0].receive[2].first + 66*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[3].size  == 66*54*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[4].first == l_tet.m_elLayout.timeGroups[0].receive[3].first + 66*54*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[4].size  == 66*19*2 );

  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[5].first == l_tet.m_elLayout.timeGroups[0].receive[4].first + 66*19*2);
  REQUIRE( l_tet.m_elLayout.timeGroups[0].receive[5].size  == 66*19*2 );

  // reset the logger
  l_conf.setToDefault();
  el::Loggers::reconfigureLogger("default", l_conf);
}

/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Tests the interface to MOAB.
 **/
#include <catch.hpp>
#define private public
#include "Moab.h"
#undef private

// directory containing the files of the unit tests
namespace edge {
  namespace test {
    extern std::string g_files;
  }
}

TEST_CASE( "Test adjacency of elements through vertices.", "[moab][elVeEl]" ) {
  // return silently if unit test files are not set
  if( edge::test::g_files == "" ) return;

  // gmsh mesh for consisting of tria3 elements
  std::string l_mshFile1 = edge::test::g_files + "/meshes/tria3.msh";

  // create non-periodic moab database
  edge::mesh::Moab l_moab1( 2 );

  // read mesh
  l_moab1.read( l_mshFile1, "" );

  int l_raw1[104*10];
  int *l_elVeEl1[104+1];
  l_elVeEl1[0] = l_raw1;

  // derive adjacency info
  l_moab1.getElVeEl( l_elVeEl1 );

  // check some sizes
  REQUIRE( (l_elVeEl1[  1] - l_elVeEl1[  0]) ==  8 );
  REQUIRE( (l_elVeEl1[  2] - l_elVeEl1[  1]) ==  9 );
  REQUIRE( (l_elVeEl1[  3] - l_elVeEl1[  2]) ==  9 );
  REQUIRE( (l_elVeEl1[  4] - l_elVeEl1[  3]) ==  8 );
  REQUIRE( (l_elVeEl1[  5] - l_elVeEl1[  4]) ==  5 );
  REQUIRE( (l_elVeEl1[  6] - l_elVeEl1[  5]) ==  5 );
  REQUIRE( (l_elVeEl1[  7] - l_elVeEl1[  6]) ==  7 );
  REQUIRE( (l_elVeEl1[  8] - l_elVeEl1[  7]) ==  9 );
  REQUIRE( (l_elVeEl1[  9] - l_elVeEl1[  8]) ==  8 );
  REQUIRE( (l_elVeEl1[ 10] - l_elVeEl1[  9]) ==  7 );

  REQUIRE( (l_elVeEl1[ 11] - l_elVeEl1[ 10]) ==  8 );
  REQUIRE( (l_elVeEl1[ 12] - l_elVeEl1[ 11]) ==  5 );
  REQUIRE( (l_elVeEl1[ 13] - l_elVeEl1[ 12]) ==  8 );
  REQUIRE( (l_elVeEl1[ 14] - l_elVeEl1[ 13]) ==  9 );
  REQUIRE( (l_elVeEl1[ 15] - l_elVeEl1[ 14]) ==  8 );
  REQUIRE( (l_elVeEl1[ 16] - l_elVeEl1[ 15]) ==  5 );
  REQUIRE( (l_elVeEl1[ 17] - l_elVeEl1[ 16]) ==  8 );
  REQUIRE( (l_elVeEl1[ 18] - l_elVeEl1[ 17]) ==  8 );
  REQUIRE( (l_elVeEl1[ 19] - l_elVeEl1[ 18]) == 12 );
  REQUIRE( (l_elVeEl1[ 20] - l_elVeEl1[ 19]) == 12 );

  REQUIRE( (l_elVeEl1[ 26] - l_elVeEl1[ 25]) == 11 );
  REQUIRE( (l_elVeEl1[ 27] - l_elVeEl1[ 26]) == 12 );
  REQUIRE( (l_elVeEl1[ 28] - l_elVeEl1[ 27]) ==  9 );
  REQUIRE( (l_elVeEl1[ 29] - l_elVeEl1[ 28]) == 11 );
  REQUIRE( (l_elVeEl1[ 30] - l_elVeEl1[ 29]) ==  7 );

  REQUIRE( (l_elVeEl1[ 31] - l_elVeEl1[ 30]) ==  9 );
  REQUIRE( (l_elVeEl1[ 32] - l_elVeEl1[ 31]) ==  7 );
  REQUIRE( (l_elVeEl1[ 33] - l_elVeEl1[ 32]) == 12 );
  REQUIRE( (l_elVeEl1[ 34] - l_elVeEl1[ 33]) ==  8 );
  REQUIRE( (l_elVeEl1[ 35] - l_elVeEl1[ 34]) == 12 );
  REQUIRE( (l_elVeEl1[ 36] - l_elVeEl1[ 35]) == 12 );
  REQUIRE( (l_elVeEl1[ 37] - l_elVeEl1[ 36]) ==  8 );
  REQUIRE( (l_elVeEl1[ 38] - l_elVeEl1[ 37]) == 13 );
  REQUIRE( (l_elVeEl1[ 39] - l_elVeEl1[ 38]) == 13 );

  REQUIRE( (l_elVeEl1[ 41] - l_elVeEl1[ 40]) == 13 );

  REQUIRE( (l_elVeEl1[ 66] - l_elVeEl1[ 65]) == 13 );
  REQUIRE( (l_elVeEl1[104] - l_elVeEl1[103]) ==  5 );

  // check some elements (see tria3.png, next to tria3.msh)
  REQUIRE( l_elVeEl1[0][0] ==  1 );
  REQUIRE( l_elVeEl1[0][1] ==  2 );
  REQUIRE( l_elVeEl1[0][2] ==  3 );
  REQUIRE( l_elVeEl1[0][3] ==  9 );
  REQUIRE( l_elVeEl1[0][4] == 10 );
  REQUIRE( l_elVeEl1[0][5] == 19 );
  REQUIRE( l_elVeEl1[0][6] == 24 );
  REQUIRE( l_elVeEl1[0][7] == 25 );

  REQUIRE( l_elVeEl1[65][ 0] == 50 );
  REQUIRE( l_elVeEl1[65][ 1] == 51 );
  REQUIRE( l_elVeEl1[65][ 2] == 54 );
  REQUIRE( l_elVeEl1[65][ 3] == 61 );
  REQUIRE( l_elVeEl1[65][ 4] == 62 );
  REQUIRE( l_elVeEl1[65][ 5] == 66 );
  REQUIRE( l_elVeEl1[65][ 6] == 67 );
  REQUIRE( l_elVeEl1[65][ 7] == 75 );
  REQUIRE( l_elVeEl1[65][ 8] == 76 );
  REQUIRE( l_elVeEl1[65][ 9] == 77 );
  REQUIRE( l_elVeEl1[65][10] == 81 );
  REQUIRE( l_elVeEl1[65][11] == 82 );
  REQUIRE( l_elVeEl1[65][12] == 83 );

  REQUIRE( l_elVeEl1[103][0] == 89  );
  REQUIRE( l_elVeEl1[103][1] == 90  );
  REQUIRE( l_elVeEl1[103][2] == 91  );
  REQUIRE( l_elVeEl1[103][3] == 93  );
  REQUIRE( l_elVeEl1[103][4] == 101 );

  // check #(adjacent elements) query
  int l_nElVeEl1 = 0;
  for( int l_el = 0; l_el < 104; l_el++ )
    l_nElVeEl1 += l_elVeEl1[l_el+1] - l_elVeEl1[l_el];
  REQUIRE( l_moab1.getNelVeEl() == l_nElVeEl1 );

  // return for MPI-builds, not supporting periodic boundaries
#ifdef PP_USE_MPI
  return;
#endif

  // create periodic moab database
  int l_bndVals1[1] = {106};
  edge::mesh::Moab l_moab2( 2, 1, l_bndVals1, 106 );

  // read mesh
  l_moab2.read( l_mshFile1, "" );

  int l_raw2[104*20];
  int *l_elVeEl2[104+1];
  l_elVeEl2[0] = l_raw2;

  l_moab2.getElVeEl( l_elVeEl2 );

  // check some sizes
  REQUIRE( (l_elVeEl2[  1] - l_elVeEl2[  0]) == 8+6 );
  REQUIRE( (l_elVeEl2[ 46] - l_elVeEl2[ 45]) == 14 );
  REQUIRE( (l_elVeEl2[ 47] - l_elVeEl2[ 46]) == 13 );
  REQUIRE( (l_elVeEl2[ 66] - l_elVeEl2[ 65]) == 13 );
  REQUIRE( (l_elVeEl2[ 93] - l_elVeEl2[ 92]) == 14 );
  REQUIRE( (l_elVeEl2[104] - l_elVeEl2[103]) == 14 );

  // check adj of some elements
  REQUIRE( l_elVeEl2[0][ 0] ==  1 );
  REQUIRE( l_elVeEl2[0][ 1] ==  2 );
  REQUIRE( l_elVeEl2[0][ 2] ==  3 );
  REQUIRE( l_elVeEl2[0][ 3] ==  9 );
  REQUIRE( l_elVeEl2[0][ 4] == 10 );
  REQUIRE( l_elVeEl2[0][ 5] == 19 );
  REQUIRE( l_elVeEl2[0][ 6] == 24 );
  REQUIRE( l_elVeEl2[0][ 7] == 25 );
  REQUIRE( l_elVeEl2[0][ 8] == 94 );
  REQUIRE( l_elVeEl2[0][ 9] == 95 );
  REQUIRE( l_elVeEl2[0][10] == 96 );
  REQUIRE( l_elVeEl2[0][11] == 97 );
  REQUIRE( l_elVeEl2[0][12] == 98 );
  REQUIRE( l_elVeEl2[0][13] == 100 );

  REQUIRE( l_elVeEl2[45][ 0] == 29 );
  REQUIRE( l_elVeEl2[45][ 1] == 31 );
  REQUIRE( l_elVeEl2[45][ 2] == 33 );
  REQUIRE( l_elVeEl2[45][ 3] == 35 );
  REQUIRE( l_elVeEl2[45][ 4] == 36 );
  REQUIRE( l_elVeEl2[45][ 5] == 40 );
  REQUIRE( l_elVeEl2[45][ 6] == 44 );
  REQUIRE( l_elVeEl2[45][ 7] == 46 );
  REQUIRE( l_elVeEl2[45][ 8] == 48 );
  REQUIRE( l_elVeEl2[45][ 9] == 49 );
  REQUIRE( l_elVeEl2[45][10] == 55 );
  REQUIRE( l_elVeEl2[45][11] == 56 );
  REQUIRE( l_elVeEl2[45][12] == 57 );
  REQUIRE( l_elVeEl2[45][13] == 58 );

  REQUIRE( l_elVeEl2[46][ 0] == 33 );
  REQUIRE( l_elVeEl2[46][ 1] == 34 );
  REQUIRE( l_elVeEl2[46][ 2] == 37 );
  REQUIRE( l_elVeEl2[46][ 3] == 44 );
  REQUIRE( l_elVeEl2[46][ 4] == 45 );
  REQUIRE( l_elVeEl2[46][ 5] == 47 );
  REQUIRE( l_elVeEl2[46][ 6] == 49 );
  REQUIRE( l_elVeEl2[46][ 7] == 52 );
  REQUIRE( l_elVeEl2[46][ 8] == 55 );
  REQUIRE( l_elVeEl2[46][ 9] == 56 );
  REQUIRE( l_elVeEl2[46][10] == 57 );
  REQUIRE( l_elVeEl2[46][11] == 58 );
  REQUIRE( l_elVeEl2[46][12] == 59 );

  REQUIRE( l_elVeEl2[65][ 0] == 50 );
  REQUIRE( l_elVeEl2[65][ 1] == 51 );
  REQUIRE( l_elVeEl2[65][ 2] == 54 );
  REQUIRE( l_elVeEl2[65][ 3] == 61 );
  REQUIRE( l_elVeEl2[65][ 4] == 62 );
  REQUIRE( l_elVeEl2[65][ 5] == 66 );
  REQUIRE( l_elVeEl2[65][ 6] == 67 );
  REQUIRE( l_elVeEl2[65][ 7] == 75 );
  REQUIRE( l_elVeEl2[65][ 8] == 76 );
  REQUIRE( l_elVeEl2[65][ 9] == 77 );
  REQUIRE( l_elVeEl2[65][10] == 81 );
  REQUIRE( l_elVeEl2[65][11] == 82 );
  REQUIRE( l_elVeEl2[65][12] == 83 );

  REQUIRE( l_elVeEl2[92][ 0] ==   4 );
  REQUIRE( l_elVeEl2[92][ 1] ==   5 );
  REQUIRE( l_elVeEl2[92][ 2] ==  11 );
  REQUIRE( l_elVeEl2[92][ 3] ==  15 );
  REQUIRE( l_elVeEl2[92][ 4] ==  71 );
  REQUIRE( l_elVeEl2[92][ 5] ==  72 );
  REQUIRE( l_elVeEl2[92][ 6] ==  84 );
  REQUIRE( l_elVeEl2[92][ 7] ==  85 );
  REQUIRE( l_elVeEl2[92][ 8] ==  86 );
  REQUIRE( l_elVeEl2[92][ 9] ==  87 );
  REQUIRE( l_elVeEl2[92][10] ==  91 );
  REQUIRE( l_elVeEl2[92][11] ==  93 );
  REQUIRE( l_elVeEl2[92][12] == 102 );
  REQUIRE( l_elVeEl2[92][13] == 103 );

  REQUIRE( l_elVeEl2[103][ 0] ==   4 );
  REQUIRE( l_elVeEl2[103][ 1] ==   5 );
  REQUIRE( l_elVeEl2[103][ 2] ==   6 );
  REQUIRE( l_elVeEl2[103][ 3] ==   7 );
  REQUIRE( l_elVeEl2[103][ 4] ==   8 );
  REQUIRE( l_elVeEl2[103][ 5] ==  11 );
  REQUIRE( l_elVeEl2[103][ 6] ==  15 );
  REQUIRE( l_elVeEl2[103][ 7] ==  89 );
  REQUIRE( l_elVeEl2[103][ 8] ==  90 );
  REQUIRE( l_elVeEl2[103][ 9] ==  91 );
  REQUIRE( l_elVeEl2[103][10] ==  92 );
  REQUIRE( l_elVeEl2[103][11] ==  93 );
  REQUIRE( l_elVeEl2[103][12] == 101 );
  REQUIRE( l_elVeEl2[103][13] == 102 );

  // check #(adjacent elements) query
  int l_nElVeEl2 = 0;
  for( int l_el = 0; l_el < 104; l_el++ )
    l_nElVeEl2 += l_elVeEl2[l_el+1] - l_elVeEl2[l_el];
  REQUIRE( l_moab2.getNelVeEl() == l_nElVeEl2 );
}

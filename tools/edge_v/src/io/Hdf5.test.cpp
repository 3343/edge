/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Tests the HDF5 interface.
 **/
#include <catch.hpp>
#define private public
#include "Hdf5.h"
#undef private

#include <cstdio>

TEST_CASE( "Tests set/get data HDF5 interface.", "[hdf5][setGetData]" ) {
  // dummy data
  std::string l_name0 = "data_0";
  edge_v::t_idx l_dataIn0[10] = { 0,  1,  2,  3,  4,
                                  10, 11, 12, 13, 14 };

  std::string l_name1 = "data_1";
  unsigned short l_dataIn1[7] = { 7, 7, 5, 5, 3,
                                  9, 2 };

  std::string l_name2 = "data_2";
  double l_dataIn2[3] = { 2.4, 3.8, 4.1 };

  // write to HDF5
  std::string l_path = std::tmpnam(nullptr);

  // create the file
  hid_t l_file = H5Fcreate( l_path.c_str(),
                            H5F_ACC_EXCL,
                            H5P_DEFAULT,
                            H5P_DEFAULT );
  REQUIRE( l_file >= 0 );
  herr_t l_err = H5Fclose( l_file );
  REQUIRE( l_err >= 0 );

  {
    edge_v::io::Hdf5 l_hdf( l_path,
                            false );

    l_hdf.set( l_name0,
               10,
               l_dataIn0 );

    l_hdf.set( l_name1,
               7,
               l_dataIn1 );

    l_hdf.set( l_name2,
               3,
               l_dataIn2 );
  }

  // read from HDF5 and check result
  edge_v::io::Hdf5 l_hdf( l_path );

  REQUIRE( l_hdf.exists( l_name0 ) );
  REQUIRE( l_hdf.exists( l_name1 ) );
  REQUIRE( l_hdf.exists( l_name2 ) );

  REQUIRE( l_hdf.nVas( l_name0 ) == 10 );
  REQUIRE( l_hdf.nVas( l_name1 ) ==  7 );
  REQUIRE( l_hdf.nVas( l_name2 ) ==  3 );

  edge_v::t_idx l_dataOut0[10] = {0};
  unsigned short l_dataOut1[7] = {0};
  double l_dataOut2[3] = {0};

  l_hdf.get( l_name0,
             l_dataOut0 );
  l_hdf.get( l_name1,
             l_dataOut1 );
  l_hdf.get( l_name2,
             l_dataOut2 );

  for( unsigned short l_va = 0; l_va < 10; l_va++ ) {
    REQUIRE( l_dataOut0[l_va] == l_dataIn0[l_va] );
  }
  for( unsigned short l_va = 0; l_va < 7; l_va++ ) {
    REQUIRE( l_dataOut1[l_va] == l_dataIn1[l_va] );
  }
  for( unsigned short l_va = 0; l_va < 3; l_va++ ) {
    REQUIRE( l_dataOut2[l_va] == l_dataIn2[l_va] );
  }
}
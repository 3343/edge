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
 * Unit tests for the initialization of kinematic sources.
 **/

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64

#include <catch.hpp>
#define private public
#include "KinematicsInit.hpp"
#undef private
#include <fstream>
#include <iostream>

TEST_CASE( "Kinematics: Derivation of local sources.", "[kinematics][local]" ) {
  /*
   * Our dummy quad4r-example:
   *
   * time group 0: 0-2
   * time group 1: 3-6
   *
   * 0, 3, 4 are inner elements
   * 1 is send-eleemnt
   * 2 is receive-element
   * 5/6 is a duplicated send element
   *
   * global ids are given by:
   *
   * 0: 3
   * 1: 2
   * 2: 4
   * 3: 1
   * 4: 0
   * 5: 5
   * 6: 5
   *
   *     0.0  0.5 1.0 1.5
   *      |   |   |   |
   *  2.0-0***3***6   |
   *      * 0 * 2 *   |
   *  1.0-1***4***7**10
   *      * 1 * 3 * 4 *
   *  0.0-2***5***8***9
   *      *5/6*
   * -1.0-11*12
   */

  // setup dummy element layout
  t_enLayout l_elLayout;
  l_elLayout.nEnts = 7;
  l_elLayout.timeGroups.resize(2);

  l_elLayout.timeGroups[0].nEntsOwn    = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[0].inner.first = 0;
  l_elLayout.timeGroups[0].inner.size  = 1;

  l_elLayout.timeGroups[0].send.resize(2);
  l_elLayout.timeGroups[0].send[0].first = 1;
  l_elLayout.timeGroups[0].send[0].size  = 1;
  l_elLayout.timeGroups[0].send[1].first = 2;
  l_elLayout.timeGroups[0].send[1].size  = 0;

  l_elLayout.timeGroups[0].receive.resize(2);
  l_elLayout.timeGroups[0].receive[0].first = 2;
  l_elLayout.timeGroups[0].receive[0].size  = 1;
  l_elLayout.timeGroups[0].receive[1].first = 3;
  l_elLayout.timeGroups[0].receive[1].size  = 0;

  l_elLayout.timeGroups[1].nEntsOwn    = 4;
  l_elLayout.timeGroups[1].nEntsNotOwn = 0;
  l_elLayout.timeGroups[1].inner.first = 3;
  l_elLayout.timeGroups[1].inner.size  = 2;

  l_elLayout.timeGroups[1].send.resize(2);
  l_elLayout.timeGroups[1].send[0].first = 5;
  l_elLayout.timeGroups[1].send[0].size  = 1;
  l_elLayout.timeGroups[1].send[1].first = 6;
  l_elLayout.timeGroups[1].send[1].size  = 1;

  l_elLayout.timeGroups[1].receive.resize(2);
  l_elLayout.timeGroups[1].receive[0].first = 7;
  l_elLayout.timeGroups[1].receive[0].size  = 0;
  l_elLayout.timeGroups[1].receive[1].first = 7;
  l_elLayout.timeGroups[1].receive[1].size  = 0;

  // assign vertex chars
  t_vertexChars l_veChars[13];
  l_veChars[ 0].coords[0] = 0.0; l_veChars[ 0].coords[1] =  2.0; l_veChars[ 0].coords[2] = 0.0;
  l_veChars[ 1].coords[0] = 0.0; l_veChars[ 1].coords[1] =  1.0; l_veChars[ 1].coords[2] = 0.0;
  l_veChars[ 2].coords[0] = 0.0; l_veChars[ 2].coords[1] =  0.0; l_veChars[ 2].coords[2] = 0.0;
  l_veChars[ 3].coords[0] = 0.5; l_veChars[ 3].coords[1] =  2.0; l_veChars[ 3].coords[2] = 0.0;
  l_veChars[ 4].coords[0] = 0.5; l_veChars[ 4].coords[1] =  1.0; l_veChars[ 4].coords[2] = 0.0;
  l_veChars[ 5].coords[0] = 0.5; l_veChars[ 5].coords[1] =  0.0; l_veChars[ 5].coords[2] = 0.0;
  l_veChars[ 6].coords[0] = 1.0; l_veChars[ 6].coords[1] =  2.0; l_veChars[ 6].coords[2] = 0.0;
  l_veChars[ 7].coords[0] = 1.0; l_veChars[ 7].coords[1] =  1.0; l_veChars[ 7].coords[2] = 0.0;
  l_veChars[ 8].coords[0] = 1.0; l_veChars[ 8].coords[1] =  0.0; l_veChars[ 8].coords[2] = 0.0;
  l_veChars[ 9].coords[0] = 1.5; l_veChars[ 9].coords[1] =  0.0; l_veChars[ 9].coords[2] = 0.0;
  l_veChars[10].coords[0] = 1.5; l_veChars[10].coords[1] =  1.0; l_veChars[10].coords[2] = 0.0;
  l_veChars[11].coords[0] = 0.0; l_veChars[11].coords[1] = -1.0; l_veChars[11].coords[2] = 0.0;
  l_veChars[12].coords[0] = 0.5; l_veChars[12].coords[1] = -1.0; l_veChars[12].coords[2] = 0.0;

  // assign el-ve adjacency
  int l_elVe[7][4];
  l_elVe[0][0] =  1; l_elVe[0][1] = 4;  l_elVe[0][2] =  3; l_elVe[0][3] = 0;
  l_elVe[1][0] =  2; l_elVe[1][1] = 5;  l_elVe[1][2] =  4; l_elVe[1][3] = 1;
  l_elVe[2][0] =  4; l_elVe[2][1] = 7;  l_elVe[2][2] =  6; l_elVe[2][3] = 3;
  l_elVe[3][0] =  5; l_elVe[3][1] = 8;  l_elVe[3][2] =  7; l_elVe[3][3] = 4;
  l_elVe[4][0] =  8; l_elVe[4][1] = 9;  l_elVe[4][2] = 10; l_elVe[4][3] = 7;
  l_elVe[5][0] = 11; l_elVe[5][1] = 12; l_elVe[5][2] =  5; l_elVe[5][3] = 2;
  l_elVe[6][0] = 11; l_elVe[6][1] = 12; l_elVe[6][2] =  5; l_elVe[6][3] = 2;

  // assign global ids
  int l_gIdsEl[7];
  l_gIdsEl[0] = 3; l_gIdsEl[1] = 2; l_gIdsEl[2] = 4; l_gIdsEl[3] = 1;
  l_gIdsEl[4] = 0; l_gIdsEl[5] = 5; l_gIdsEl[6] = 5;

  // source crds
  double   l_srcCrds[20][2];
  // pointer to source locations
  double (*l_srcPts[8])[2] = { l_srcCrds, l_srcCrds, l_srcCrds, l_srcCrds,
                               l_srcCrds, l_srcCrds, l_srcCrds, l_srcCrds };

  // elements holding source terms
  std::vector< int > l_srcEl[8];
  // ids of sources assigned to elements
  std::vector< int > l_srcId[8];

  // number of sources
  int l_nSrcs[8];
  l_nSrcs[0] = 0;

  // dry-run without source terms
  edge::elastic::setups::KinematicsInit< int, double, double, QUAD4R, 1, 1, 1 >::getLocal( l_nSrcs,
                                                                                           l_srcPts,
                                                                                           l_elLayout,
                                                                                           l_elVe,
                                                                                           l_veChars,
                                                                                           l_gIdsEl,
                                                                                           l_srcEl, l_srcId );
  for( unsigned short l_is = 0; l_is < 8; l_is++ ) {
    REQUIRE( l_srcEl[l_is].size() == 0 );
    REQUIRE( l_srcId[l_is].size() == 0 );
  }


  // single source
  l_srcCrds[0][0] = 0.25; l_srcCrds[0][1] = 0.25;
  l_nSrcs[0] = 1; l_nSrcs[1] = 0; l_nSrcs[2] = 0; l_nSrcs[3] = 0;
  l_nSrcs[4] = 0; l_nSrcs[5] = 0; l_nSrcs[6] = 0; l_nSrcs[7] = 0;
  edge::elastic::setups::KinematicsInit< int, double, double, QUAD4R, 1, 1, 1 >::getLocal( l_nSrcs,
                                                                                           l_srcPts,
                                                                                           l_elLayout,
                                                                                           l_elVe,
                                                                                           l_veChars,
                                                                                           l_gIdsEl,
                                                                                           l_srcEl, l_srcId );
  // check sources
  REQUIRE( l_srcEl[0][0] == 1 );
  REQUIRE( l_srcId[0][0] == 0 );
  for( unsigned short l_is = 1; l_is < 8; l_is++ ) {
    REQUIRE( l_srcEl[l_is].size() == 0 );
    REQUIRE( l_srcId[l_is].size() == 0 );
  }


  // single source one for each of four sims
  l_srcCrds[0][0] = 0.25; l_srcCrds[0][1] = 0.25;
  l_nSrcs[0] = 1; l_nSrcs[1] = 1; l_nSrcs[2] = 1; l_nSrcs[3] = 1;
  l_nSrcs[4] = 0; l_nSrcs[5] = 0; l_nSrcs[6] = 0; l_nSrcs[7] = 0;
  edge::elastic::setups::KinematicsInit< int, double, double, QUAD4R, 1, 4, 1 >::getLocal( l_nSrcs,
                                                                                           l_srcPts,
                                                                                           l_elLayout,
                                                                                           l_elVe,
                                                                                           l_veChars,
                                                                                           l_gIdsEl,
                                                                                           l_srcEl, l_srcId );
  // check sources
  for( unsigned short l_is = 0; l_is < 4; l_is++ ) {
    REQUIRE( l_srcEl[l_is].size() == 1 );
    REQUIRE( l_srcId[l_is].size() == 1 );

    REQUIRE( l_srcEl[0][0] == 1 );
    REQUIRE( l_srcId[0][0] == 0 );
  }


  // multiple souces
  // (0.5, 1.0)   -> gId 1, el 3
  // (0.25, -0.5) -> gId 5, el 5/6
  // (0, 0.25)    -> gId 2, el 1
  // (0.75, 0.5)  -> gId 1, el 3
  // (0.75, 1.5)  -> gId 4, el 2
  l_srcCrds[0][0] = 0.5;
  l_srcCrds[0][1] = 1.0;
  l_srcCrds[1][0] = 0.25;
  l_srcCrds[1][1] = -0.5;
  l_srcCrds[2][0] = 0;
  l_srcCrds[2][1] = 0.25;
  l_srcCrds[3][0] = 0.75;
  l_srcCrds[3][1] = 0.5;
  l_srcCrds[4][0] = 0.75;
  l_srcCrds[4][1] = 1.5;
  l_nSrcs[0] = 5;
  edge::elastic::setups::KinematicsInit< int, double, double, QUAD4R, 1, 1, 1 >::getLocal( l_nSrcs,
                                                                                           l_srcPts,
                                                                                           l_elLayout,
                                                                                           l_elVe,
                                                                                           l_veChars,
                                                                                           l_gIdsEl,
                                                                                           l_srcEl, l_srcId );
  REQUIRE( l_srcEl[0].size() == 5 );
  REQUIRE( l_srcId[0].size() == 5 );

  REQUIRE( l_srcEl[0][0] == 3 );
  REQUIRE( l_srcId[0][0] == 0 );
  REQUIRE( l_srcEl[0][1] == 5 );
  REQUIRE( l_srcId[0][1] == 1 );
  REQUIRE( l_srcEl[0][2] == 6 );
  REQUIRE( l_srcId[0][2] == 1 );
  REQUIRE( l_srcEl[0][3] == 1 );
  REQUIRE( l_srcId[0][3] == 2 );
  REQUIRE( l_srcEl[0][4] == 3 );
  REQUIRE( l_srcId[0][4] == 3 );
}

#if defined PP_HAS_NETCDF
#include "impl/elastic/io/Nrf.h"
#include "data/SparseEntities.hpp"

TEST_CASE( "Kinematics: Init NRF.", "[kinematics][initNrf]" ) {
  std::string l_ncFile;

  /*
   * Our dummy quad4r-example:
   *
   * time group 0: 0-2
   * time group 1: 3-6
   *
   * 0, 3, 4 are inner elements
   * 1 is send-eleemnt
   * 2 is receive-element
   * 5/6 is a duplicated send element
   *
   * global ids are given by:
   *
   * 0: 3
   * 1: 2
   * 2: 4
   * 3: 1
   * 4: 0
   * 5: 5
   * 6: 5
   *
   *     0.0  0.5 1.0 1.5
   *      |   |   |   |
   *  2.0-0***3***6   |
   *      * 0 * 2 *   |
   *  1.0-1***4***7**10
   *      * 1 * 3 * 4 *
   *  0.0-2***5***8***9
   *      *5/6*
   * -1.0-11*12
   */

  // setup dummy element layout
  t_enLayout l_elLayout;
  l_elLayout.nEnts = 7;
  l_elLayout.timeGroups.resize(2);

  l_elLayout.timeGroups[0].nEntsOwn    = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[0].inner.first = 0;
  l_elLayout.timeGroups[0].inner.size  = 1;

  l_elLayout.timeGroups[0].neRanks.resize(2);
  l_elLayout.timeGroups[0].neRanks[0] = 1;
  l_elLayout.timeGroups[0].neRanks[1] = 3;

  l_elLayout.timeGroups[0].send.resize(2);
  l_elLayout.timeGroups[0].send[0].first = 1;
  l_elLayout.timeGroups[0].send[0].size  = 1;
  l_elLayout.timeGroups[0].send[1].first = 2;
  l_elLayout.timeGroups[0].send[1].size  = 0;

  l_elLayout.timeGroups[0].receive.resize(2);
  l_elLayout.timeGroups[0].receive[0].first = 2;
  l_elLayout.timeGroups[0].receive[0].size  = 1;
  l_elLayout.timeGroups[0].receive[1].first = 3;
  l_elLayout.timeGroups[0].receive[1].size  = 0;

  l_elLayout.timeGroups[0].neTgs.resize(2);
  l_elLayout.timeGroups[0].neTgs[0] = 0;
  l_elLayout.timeGroups[0].neTgs[1] = 2;

  l_elLayout.timeGroups[1].nEntsOwn    = 4;
  l_elLayout.timeGroups[1].nEntsNotOwn = 0;
  l_elLayout.timeGroups[1].inner.first = 3;
  l_elLayout.timeGroups[1].inner.size  = 2;

  l_elLayout.timeGroups[1].neRanks.resize(2);
  l_elLayout.timeGroups[1].neRanks[0] = 0;
  l_elLayout.timeGroups[1].neRanks[1] = 2;

  l_elLayout.timeGroups[1].send.resize(2);
  l_elLayout.timeGroups[1].send[0].first = 5;
  l_elLayout.timeGroups[1].send[0].size  = 1;
  l_elLayout.timeGroups[1].send[1].first = 6;
  l_elLayout.timeGroups[1].send[1].size  = 1;

  l_elLayout.timeGroups[1].receive.resize(2);
  l_elLayout.timeGroups[1].receive[0].first = 7;
  l_elLayout.timeGroups[1].receive[0].size  = 0;
  l_elLayout.timeGroups[1].receive[1].first = 7;
  l_elLayout.timeGroups[1].receive[1].size  = 0;

  l_elLayout.timeGroups[1].neTgs.resize(2);
  l_elLayout.timeGroups[1].neTgs[0] = 1;
  l_elLayout.timeGroups[1].neTgs[1] = 0;

  // assign vertex chars
  t_vertexChars l_veChars[13];
  l_veChars[ 0].coords[0] = 0.0; l_veChars[ 0].coords[1] =  2.0; l_veChars[ 0].coords[2] = 0.0;
  l_veChars[ 1].coords[0] = 0.0; l_veChars[ 1].coords[1] =  1.0; l_veChars[ 1].coords[2] = 0.0;
  l_veChars[ 2].coords[0] = 0.0; l_veChars[ 2].coords[1] =  0.0; l_veChars[ 2].coords[2] = 0.0;
  l_veChars[ 3].coords[0] = 0.5; l_veChars[ 3].coords[1] =  2.0; l_veChars[ 3].coords[2] = 0.0;
  l_veChars[ 4].coords[0] = 0.5; l_veChars[ 4].coords[1] =  1.0; l_veChars[ 4].coords[2] = 0.0;
  l_veChars[ 5].coords[0] = 0.5; l_veChars[ 5].coords[1] =  0.0; l_veChars[ 5].coords[2] = 0.0;
  l_veChars[ 6].coords[0] = 1.0; l_veChars[ 6].coords[1] =  2.0; l_veChars[ 6].coords[2] = 0.0;
  l_veChars[ 7].coords[0] = 1.0; l_veChars[ 7].coords[1] =  1.0; l_veChars[ 7].coords[2] = 0.0;
  l_veChars[ 8].coords[0] = 1.0; l_veChars[ 8].coords[1] =  0.0; l_veChars[ 8].coords[2] = 0.0;
  l_veChars[ 9].coords[0] = 1.5; l_veChars[ 9].coords[1] =  0.0; l_veChars[ 9].coords[2] = 0.0;
  l_veChars[10].coords[0] = 1.5; l_veChars[10].coords[1] =  1.0; l_veChars[10].coords[2] = 0.0;
  l_veChars[11].coords[0] = 0.0; l_veChars[11].coords[1] = -1.0; l_veChars[11].coords[2] = 0.0;
  l_veChars[12].coords[0] = 0.5; l_veChars[12].coords[1] = -1.0; l_veChars[12].coords[2] = 0.0;

  // element chars
  t_elementChars l_elChars[7];

  // assign el-ve adjacency
  int l_elVe[7][4];

  l_elVe[0][0] =  1; l_elVe[0][1] = 4;  l_elVe[0][2] =  3; l_elVe[0][3] = 0;
  l_elVe[1][0] =  2; l_elVe[1][1] = 5;  l_elVe[1][2] =  4; l_elVe[1][3] = 1;
  l_elVe[2][0] =  4; l_elVe[2][1] = 7;  l_elVe[2][2] =  6; l_elVe[2][3] = 3;
  l_elVe[3][0] =  5; l_elVe[3][1] = 8;  l_elVe[3][2] =  7; l_elVe[3][3] = 4;
  l_elVe[4][0] =  8; l_elVe[4][1] = 9;  l_elVe[4][2] = 10; l_elVe[4][3] = 7;
  l_elVe[5][0] = 11; l_elVe[5][1] = 12; l_elVe[5][2] =  5; l_elVe[5][3] = 2;
  l_elVe[6][0] = 11; l_elVe[6][1] = 12; l_elVe[6][2] =  5; l_elVe[6][3] = 2;

  // assign global ids
  int l_gIdsEl[7];
  l_gIdsEl[0] = 3; l_gIdsEl[1] = 2; l_gIdsEl[2] = 4; l_gIdsEl[3] = 1;
  l_gIdsEl[4] = 0; l_gIdsEl[5] = 5; l_gIdsEl[6] = 5;

  // dynamic memory
  edge::data::Dynamic l_mem;

  // velocity model
  t_bgPars l_vel[7];
  for( unsigned short l_el = 0; l_el < 7; l_el++ ) {
    l_vel[l_el].mu = l_el+1;
    l_vel[l_el].lam = l_el+2;
  }

  /*
   * single point source
   */
  l_ncFile = "examples/ci/elastic/sources/unit_tests/kinematics_init_1.nc";

  // only run unit tests if input file exists
  if( std::ifstream(l_ncFile) ) {
    // NRF reader with 2-source buffer
    edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

    // init NRC-reader
    l_nrf1.init( l_ncFile );

    // reset element chars
    for( unsigned int l_el = 0; l_el < 7; l_el++ ) l_elChars[l_el].spType = 1;

    // inverse mass matrix
    double l_massI[4] = {0.1, 0.2, 0.3, 0.4};

    // first pt-sources of the sparse source elements
    int (*l_k1Sp)[1] = nullptr;
    // source solvers
    edge::elastic::solvers::t_Kinematics< 2,
                                          1,
                                          1,
                                          double,
                                          int > l_k1So[1];

    // get solvers
    edge::elastic::setups::KinematicsInit<
      int,
      double,
      double,
      QUAD4R,
      1,
      1,
      1 >::solvers( l_elLayout,
                    l_elVe,
                    (int_spType) 332,
                    l_veChars,
                    l_gIdsEl,
                    l_nrf1,
                    l_massI,
                    l_vel,
                    l_elChars,
                    l_mem,
                    l_k1Sp,
                    l_k1So );
    for( unsigned short l_el = 0; l_el < 7; l_el++ ) {
      if( l_el != 1 ) REQUIRE( l_elChars[l_el].spType == 1 );
      else REQUIRE( l_elChars[1].spType == 333 );
    }

    // add another kinematic source description
    l_nrf1.init( l_ncFile );

    // reset element chars
    for( unsigned int l_el = 0; l_el < 7; l_el++ ) l_elChars[l_el].spType = 1;

    // first pt-sources of the sparse source elements
    int (*l_ki2Sp)[2] = nullptr;
    // source solvers
    edge::elastic::solvers::t_Kinematics< 2,
                                          1,
                                          1,
                                          double,
                                          int > l_ki2So[2];

    // get solvers
    edge::elastic::setups::KinematicsInit<
      int,
      double,
      double,
      QUAD4R,
      1,
      2,
      1 >::solvers( l_elLayout,
                    l_elVe,
                    (int_spType) 332,
                    l_veChars,
                    l_gIdsEl,
                    l_nrf1,
                    l_massI,
                    l_vel,
                    l_elChars,
                    l_mem,
                    l_ki2Sp,
                    l_ki2So );

    for( unsigned short l_el = 0; l_el < 7; l_el++ ) {
      if( l_el != 1 ) REQUIRE( l_elChars[l_el].spType == 1 );
      else REQUIRE( l_elChars[1].spType == 333 );
    }
  }

  l_ncFile = "examples/ci/elastic/sources/unit_tests/kinematics_init_3.nc";
  if( std::ifstream(l_ncFile) ) {
    // NRF reader with 2-source buffer
    edge::elastic::io::Nrf< 2 > l_nrf2( 2 );

    // init NRC-reader
    l_nrf2.init( l_ncFile );
    l_nrf2.init( l_ncFile );
    l_nrf2.init( l_ncFile );

    // reset element chars
    for( unsigned int l_el = 0; l_el < 7; l_el++ ) l_elChars[l_el].spType = 1;

    // inverse mass matrix
    double l_massI[9] = {1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    // first pt-sources of the sparse source elements
    int (*l_ki3Sp)[3]= nullptr;
    // source solvers
    edge::elastic::solvers::t_Kinematics< 2,
                                          9,
                                          1,
                                          double,
                                          int > l_ki3So[3];

    // get solvers
    edge::elastic::setups::KinematicsInit<
      int,
      double,
      double,
      QUAD4R,
      3,
      3,
      1 >::solvers( l_elLayout,
                    l_elVe,
                    (int_spType) 332,
                    l_veChars,
                    l_gIdsEl,
                    l_nrf2,
                    l_massI,
                    l_vel,
                    l_elChars,
                    l_mem,
                    l_ki3Sp,
                    l_ki3So );
    REQUIRE( l_elChars[0].spType ==   1 );
    REQUIRE( l_elChars[1].spType == 333 );
    REQUIRE( l_elChars[2].spType ==   1 ); // receive element -> no src term
    REQUIRE( l_elChars[3].spType == 333 );
    REQUIRE( l_elChars[4].spType ==   1 );
    REQUIRE( l_elChars[5].spType == 333 ); // duplicated send element -> src in 5/6
    REQUIRE( l_elChars[6].spType == 333 ); // duplicated send element -> src in 5/6

    // check first entries and size
    for( unsigned short l_is = 0; l_is < 3; l_is++ ) {
      REQUIRE( l_ki3So[l_is].nSrcs == 5 );
      REQUIRE( l_ki3So[l_is].aSlip[0] == false );
      REQUIRE( l_ki3So[l_is].aSlip[1] == true  );

      // check sizes (Remark: the ordering follows the local ids)
      REQUIRE( l_ki3So[l_is].first[1][0] ==  0 ); // global src 3
      REQUIRE( l_ki3So[l_is].first[1][1] ==  4 ); // global src 0
      REQUIRE( l_ki3So[l_is].first[1][2] ==  6 ); // global src 6
      REQUIRE( l_ki3So[l_is].first[1][3] ==  8 ); // global src 1
      REQUIRE( l_ki3So[l_is].first[1][4] ==  9 ); // global src 1
      REQUIRE( l_ki3So[l_is].first[1][5] == 10 ); // ghost entry

      // check onset times
      REQUIRE( l_ki3So[l_is].onSet[0] == Approx( 3.0 ) );
      REQUIRE( l_ki3So[l_is].onSet[1] == Approx( 3.5 ) );
      REQUIRE( l_ki3So[l_is].onSet[2] == Approx( 2.5 ) );
      REQUIRE( l_ki3So[l_is].onSet[3] == Approx( 1.5 ) );
      REQUIRE( l_ki3So[l_is].onSet[4] == Approx( 1.5 ) );

      // check dt of samples
      REQUIRE( l_ki3So[l_is].dt[0] == Approx( 0.1 ) );
      REQUIRE( l_ki3So[l_is].dt[1] == Approx( 0.1 ) );
      REQUIRE( l_ki3So[l_is].dt[2] == Approx( 0.1 ) );
      REQUIRE( l_ki3So[l_is].dt[3] == Approx( 0.2 ) );
      REQUIRE( l_ki3So[l_is].dt[4] == Approx( 0.2 ) );

      // check slip-rates
      REQUIRE( l_ki3So[l_is].sr[1][0][0] == Approx( -0.4 ) );
      REQUIRE( l_ki3So[l_is].sr[1][1][0] == Approx(  0.5 ) );
      REQUIRE( l_ki3So[l_is].sr[1][2][0] == Approx( -0.6 ) );
      REQUIRE( l_ki3So[l_is].sr[1][3][0] == Approx(  0.7 ) );

      REQUIRE( l_ki3So[l_is].sr[1][4][0] == Approx(  0.1 ) );
      REQUIRE( l_ki3So[l_is].sr[1][5][0] == Approx( -0.2 ) );

      REQUIRE( l_ki3So[l_is].sr[1][6][0] == Approx( -0.8 ) );
      REQUIRE( l_ki3So[l_is].sr[1][7][0] == Approx(  0.9 ) );

      REQUIRE( l_ki3So[l_is].sr[1][8][0] == Approx(  0.3 ) );

      REQUIRE( l_ki3So[l_is].sr[1][9][0] == Approx(  0.3 ) );

      // check coefficiencts
      double l_coeff;
      l_coeff = 12000 * 2.0 * 3368000000 * 3.0 * 0.3;
      REQUIRE( l_ki3So[l_is].sSca[1][0][0][0] == Approx(l_coeff) ); // m_11
      l_coeff = 12000 * 2.0 * 3368000000 * 3.0 * 1.3;
      REQUIRE( l_ki3So[l_is].sSca[1][0][1][0] == Approx(l_coeff) ); // m_22
      l_coeff = 12000       * 3368000000 * (3.0 * 1.3 + 3.0 * 0.3 );
      REQUIRE( l_ki3So[l_is].sSca[1][0][2][0] == Approx(l_coeff) ); // m_12

      l_coeff = 10000 * 2.0 * 3168000000 * 1.0 * 0.1;
      REQUIRE( l_ki3So[l_is].sSca[1][1][0][0] == Approx(l_coeff) ); // m_11
      l_coeff = 10000 * 2.0 * 3168000000 * 5.0 * 1.5;
      REQUIRE( l_ki3So[l_is].sSca[1][1][1][0] == Approx(l_coeff) ); // m_22
      l_coeff = 10000       * 3168000000 * ( 1.0 * 1.5 + 5.0 * 0.1 );
      REQUIRE( l_ki3So[l_is].sSca[1][1][2][0] == Approx(l_coeff) ); // m_12

      l_coeff = 13000 * 2.0 * 3468000000 * 4.0 * 0.4;
      REQUIRE( l_ki3So[l_is].sSca[1][2][0][0] == Approx(l_coeff) ); // m_11
      l_coeff = 13000 * 2.0 * 3468000000 * 2.0 * 1.2;
      REQUIRE( l_ki3So[l_is].sSca[1][2][1][0] == Approx(l_coeff) ); // m_22
      l_coeff = 13000       * 3468000000 * ( 4.0 * 1.2 + 0.4 * 2.0 );
      REQUIRE( l_ki3So[l_is].sSca[1][2][2][0] == Approx(l_coeff) ); // m_12

      l_coeff = 11000 * 2.0 * 3268000000 * 2.0 * 0.2;
      REQUIRE( l_ki3So[l_is].sSca[1][3][0][0] == Approx(l_coeff) ); // m_11
      l_coeff = 11000 * 2.0 * 3268000000 * 4.0 * 1.4;
      REQUIRE( l_ki3So[l_is].sSca[1][3][1][0] == Approx(l_coeff) ); // m_22
      l_coeff = 11000       * 3268000000 * ( 2.0 * 1.4 + 0.2 * 4.0 );
      REQUIRE( l_ki3So[l_is].sSca[1][3][2][0] == Approx(l_coeff) ); // m_12

      l_coeff = 11000 * 2.0 * 3268000000 * 2.0 * 0.2;
      REQUIRE( l_ki3So[l_is].sSca[1][4][0][0] == Approx(l_coeff) ); // m_11
      l_coeff = 11000 * 2.0 * 3268000000 * 4.0 * 1.4;
      REQUIRE( l_ki3So[l_is].sSca[1][4][1][0] == Approx(l_coeff) ); // m_22
      l_coeff = 11000       * 3268000000 * ( 2.0 * 1.4 + 0.2 * 4.0 );
      REQUIRE( l_ki3So[l_is].sSca[1][4][2][0] == Approx(l_coeff) ); // m_12
    }

    // check dense ids of the sources
    for( unsigned short l_is = 0; l_is < 3; l_is++ ) {
      REQUIRE( l_ki3So[l_is].soElDe[0] == 1 );
      REQUIRE( l_ki3So[l_is].soElDe[1] == 3 );
      REQUIRE( l_ki3So[l_is].soElDe[2] == 3 );
      REQUIRE( l_ki3So[l_is].soElDe[3] == 5 );
      REQUIRE( l_ki3So[l_is].soElDe[4] == 6 );
    }

    // check the mapping from sparse source-elements to the point sources
    for( unsigned short l_is = 0; l_is < 3; l_is++ ) {
      REQUIRE( l_ki3Sp[0][l_is] == 0 ); // gsrc 3
      REQUIRE( l_ki3Sp[1][l_is] == 1 ); // gsrc 0/6
      REQUIRE( l_ki3Sp[2][l_is] == 3 ); // gsrc 5
      REQUIRE( l_ki3Sp[3][l_is] == 4 ); // gsrc 5
      REQUIRE( l_ki3Sp[4][l_is] == 5 ); // ghost
    }

    // derive a sparse entity layout of the elements with sources
    t_enLayout l_srcLayout;
    edge::data::SparseEntities::denseToSparse( (int_spType) 332,
                                               l_elChars,
                                               l_elLayout,
                                               l_srcLayout );

    // check the layout
    REQUIRE( l_srcLayout.nEnts                          == 4 );
    REQUIRE( l_srcLayout.timeGroups[0].nEntsOwn         == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].nEntsNotOwn      == 0 );
    REQUIRE( l_srcLayout.timeGroups[0].inner.first      == 0 );
    REQUIRE( l_srcLayout.timeGroups[0].inner.size       == 0 );

    REQUIRE( l_srcLayout.timeGroups[0].send[0].first    == 0 );
    REQUIRE( l_srcLayout.timeGroups[0].send[0].size     == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].send[1].first    == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].send[1].size     == 0 );

    REQUIRE( l_srcLayout.timeGroups[0].receive[0].first == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].receive[0].size  == 0 );
    REQUIRE( l_srcLayout.timeGroups[0].receive[1].first == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].receive[1].size  == 0 );

    REQUIRE( l_srcLayout.timeGroups[0].neRanks[0]       == 1 );
    REQUIRE( l_srcLayout.timeGroups[0].neRanks[1]       == 3 );

    REQUIRE( l_srcLayout.timeGroups[0].neTgs[0]         == 0 );
    REQUIRE( l_srcLayout.timeGroups[0].neTgs[1]         == 2 );


    REQUIRE( l_srcLayout.timeGroups[1].nEntsOwn         == 3 );
    REQUIRE( l_srcLayout.timeGroups[1].nEntsNotOwn      == 0 );
    REQUIRE( l_srcLayout.timeGroups[1].inner.first      == 1 );
    REQUIRE( l_srcLayout.timeGroups[1].inner.size       == 1 );

    REQUIRE( l_srcLayout.timeGroups[1].send[0].first    == 2 );
    REQUIRE( l_srcLayout.timeGroups[1].send[0].size     == 1 );
    REQUIRE( l_srcLayout.timeGroups[1].send[1].first    == 3 );
    REQUIRE( l_srcLayout.timeGroups[1].send[1].size     == 1 );

    REQUIRE( l_srcLayout.timeGroups[1].receive[0].first == 4 );
    REQUIRE( l_srcLayout.timeGroups[1].receive[0].size  == 0 );
    REQUIRE( l_srcLayout.timeGroups[1].receive[1].first == 4 );
    REQUIRE( l_srcLayout.timeGroups[1].receive[1].size  == 0 );

    REQUIRE( l_srcLayout.timeGroups[1].neRanks[0]       == 0 );
    REQUIRE( l_srcLayout.timeGroups[1].neRanks[1]       == 2 );

    REQUIRE( l_srcLayout.timeGroups[1].neTgs[0]         == 1 );
    REQUIRE( l_srcLayout.timeGroups[1].neTgs[1]         == 0 );
  }
}

#endif

#endif

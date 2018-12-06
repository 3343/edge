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
 * Unit tests of the dynamic load balancing.
 **/
#include <catch.hpp>
#define private public
#include "LoadBalancing.h"
#undef private

TEST_CASE( "Load balancing: balance function.", "[balance][loadBalancing]" ) {
  // init the load balancing with 7 workers
  edge::parallel::LoadBalancing l_lb1;
  l_lb1.init( 7 );

  // register four work regions
  l_lb1.regWrkRgn( 0,
                   200,
                   50211 );

  l_lb1.regWrkRgn( 0,
                   301,
                   0 );

  l_lb1.regWrkRgn( 2,
                   302,
                   12330 );

  l_lb1.regWrkRgn( 3,
                   111,
                   5 );

  /*
   * equal distribution (nothing elapsed)
   */
  l_lb1.balanceWrkRgn( 0 );

  // first work region is empty
  for( unsigned short l_wo = 0; l_wo < 7; l_wo++ ) {
    REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[l_wo].first == 301 );
    REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[l_wo].size  ==   0 );
  }

  // second work region is divisible by the workers
  l_lb1.balanceWrkRgn( 1 );

  REQUIRE( l_lb1.m_wrkRgns[1].wrkPkgs[0].first ==  200 );
  REQUIRE( l_lb1.m_wrkRgns[1].wrkPkgs[0].size  == 7173 );
  for( unsigned short l_wo = 1; l_wo < 7; l_wo++ ) {
    REQUIRE( l_lb1.m_wrkRgns[1].wrkPkgs[l_wo].first == l_lb1.m_wrkRgns[1].wrkPkgs[l_wo-1].first + 7173 );
    REQUIRE( l_lb1.m_wrkRgns[1].wrkPkgs[l_wo].size  == 7173 );
  }

  // third work region is not divisible by workers, first three have an additional entry
  l_lb1.balanceWrkRgn( 2 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[0].first ==   302 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[0].size  ==  1762 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[1].first ==  2064 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[1].size  ==  1762 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[2].first ==  3826 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[2].size  ==  1762 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[3].first ==  5588 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[3].size  ==  1761 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[4].first ==  7349 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[4].size  ==  1761 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[5].first ==  9110 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[5].size  ==  1761 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[6].first == 10871 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[6].size  ==  1761 );

  // fourth work region has less work than workers
  l_lb1.balanceWrkRgn( 3 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[0].first == 111 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[0].size  ==   1 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[1].first == 112 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[1].size  ==   1 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[2].first == 113 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[2].size  ==   1 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[3].first == 114 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[3].size  ==   1 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[4].first == 115 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[4].size  ==   1 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[5].first == 116 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[5].size  ==   0 );

  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[6].first == 116 );
  REQUIRE( l_lb1.m_wrkRgns[3].wrkPkgs[6].size  ==   0 );

  /*
   * distribution with imbalance
   *
   * total load summed load of 100s for all workers.
   */
  for( unsigned short l_rg = 0; l_rg < 4; l_rg++ ) {
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[0].timer.m_elapsed = 14.29;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[1].timer.m_elapsed = 15;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[2].timer.m_elapsed = 13;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[3].timer.m_elapsed = 11.71;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[4].timer.m_elapsed = 18;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[5].timer.m_elapsed = 5;
    l_lb1.m_wrkRgns[l_rg].wrkPkgs[6].timer.m_elapsed = 23;
  }

  l_lb1.balanceWrkRgn( 0 );

  for( unsigned short l_wo = 0; l_wo < 7; l_wo++ ) {
    REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[l_wo].first == 301 );
    REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[l_wo].size  ==   0 );
  }

  /*
   * third work region: 12330 entities
   *
   * initial throughput balanced
   * 1762    123.31     1443     +1
   * 1762    117.47     1375     +1
   * 1762    135.54     1586     +1
   * 1761    150.38     1760
   * 1761     97.83     1145
   * 1761    352.20     4122
   * 1761     76.57      896
   *        1053.3     12327
   */
  l_lb1.balanceWrkRgn( 2 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[0].first ==   302 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[0].size  ==  1444 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[1].first ==  1746 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[1].size  ==  1376 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[2].first ==  3122 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[2].size  ==  1587 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[3].first ==  4709 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[3].size  ==  1760 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[4].first ==  6469 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[4].size  ==  1145 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[5].first ==  7614 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[5].size  ==  4122 );

  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[6].first == 11736 );
  REQUIRE( l_lb1.m_wrkRgns[2].wrkPkgs[6].size  ==   896 );
}

TEST_CASE( "Load balancing: sparse entities.", "[spEn][loadBalancing]" ) {
  // our sparse types
  struct {
    int spType;
  } l_chars[100];

  for( unsigned short l_en = 0; l_en < 100; l_en++ ) l_chars[l_en].spType = 0;

  l_chars[ 4].spType = 7; // 0
  l_chars[ 9].spType = 7; // 1
  l_chars[13].spType = 7; // 2
  l_chars[24].spType = 7; // 3
  l_chars[29].spType = 7; // 4
  l_chars[38].spType = 7; // 5
  l_chars[41].spType = 7; // 6
  l_chars[61].spType = 7; // 7
  l_chars[77].spType = 7; // 8
  l_chars[93].spType = 7; // 9

  l_chars[ 8].spType = 8; // 0
  l_chars[21].spType = 8; // 1
  l_chars[59].spType = 8; // 2
  l_chars[60].spType = 8; // 3
  l_chars[62].spType = 8; // 4
  l_chars[81].spType = 8; // 5
  l_chars[94].spType = 8; // 6

  // init the load balancing with 3 workers
  edge::parallel::LoadBalancing l_lb1;
  l_lb1.init( 3 );

  int l_spTypes[2] = {7,8};

  // register first work region
  l_lb1.regWrkRgn( 0,
                   15,
                   75,
                   2,
                   l_spTypes,
                   l_chars );

  REQUIRE( l_lb1.m_wrkRgns[0].firstSp[0] == 3 );
  REQUIRE( l_lb1.m_wrkRgns[0].firstSp[1] == 1 );

  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[0].firstSp[0] == 3 );
  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[1].firstSp[0] == 6 );
  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[2].firstSp[0] == 8 );

  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[0].firstSp[1] == 1 );
  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[1].firstSp[1] == 2 );
  REQUIRE( l_lb1.m_wrkRgns[0].wrkPkgs[2].firstSp[1] == 5 );
}
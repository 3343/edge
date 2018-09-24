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
 * Unit tests for the point sources.
 **/

#include <catch.hpp>
#define private public
#include "PointSources.hpp"
#undef private

TEST_CASE( "Point sources, applied to quad4r elements.", "[pointSources][quad4r]" ) {
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
  struct{
    double coords[3];
  } l_veChars[13];
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

  struct{
    int spType;
  } l_elChars[7] = { {0},{0},{0},{0},{0},{0},{0} };

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

  std::string l_h5File = "cont/unit_tests/elastic/sources/point_sources_2d_0.h5";

  // DOFs
  float l_dofs[7][5][4][3];
  for( unsigned short l_el = 0; l_el < 7; l_el++ )
    for( unsigned short l_qt = 0; l_qt < 5; l_qt++ )
      for( unsigned short l_md = 0; l_md < 4; l_md++ )
        for( unsigned short l_cr = 0; l_cr < 3; l_cr++ )
          l_dofs[l_el][l_qt][l_md][l_cr] = 0;

  float l_massI[4] = {1.0, 1.0, 1.0, 1.0};

  // return silently if file does not exist
  if( std::ifstream(l_h5File) ) {
    std::string l_h5Files[3];
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ )
      l_h5Files[l_cr] = l_h5File;

    edge::data::Dynamic l_dynMem;

    edge::seismic::setups::PointSources<
     float,
     int,
     QUAD4R,
     2,
     3 > l_pss;
    l_pss.init( l_h5Files,
                l_gIdsEl,
                99,
                l_elLayout,
                l_elVe,
                l_veChars,
                l_massI,
                l_elChars,
                l_dynMem );

    /*
     * check the first example, which defines two point forces.
     */
    REQUIRE( l_elChars[0].spType == 99 );
    REQUIRE( l_elChars[1].spType ==  0 );
    REQUIRE( l_elChars[2].spType ==  0 );
    REQUIRE( l_elChars[3].spType ==  0 );
    REQUIRE( l_elChars[4].spType == 99 );
    REQUIRE( l_elChars[5].spType ==  0 );
    REQUIRE( l_elChars[6].spType ==  0 );

    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_pss.m_elSpPs[0][l_cr] == 0 );
      REQUIRE( l_pss.m_elSpPs[1][l_cr] == 1 );
      REQUIRE( l_pss.m_elSpPs[2][l_cr] == 2 ); // ghost

      REQUIRE( l_pss.m_psFused[l_cr].nPts == 2 );

      REQUIRE( l_pss.m_psFused[l_cr].el[0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].el[1] == 4 );

      REQUIRE( l_pss.m_psFused[l_cr].times[0] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].times[1] == Approx(2.3) );

      REQUIRE( l_pss.m_psFused[l_cr].dts[0] == Approx(0.001) );
      REQUIRE( l_pss.m_psFused[l_cr].dts[1] == Approx(0.002) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[0][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][3] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][4] == Approx(2.3) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[1][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][3] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][4] == Approx(3.3) );

      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[1] == 5 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[2] == 8 );

      REQUIRE( l_pss.m_psFused[l_cr].tss[0] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[1] == Approx(0.1) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[2] == Approx(0.2) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[3] == Approx(0.1) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[4] == Approx(0.15) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[5] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[6] == Approx(0.3) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[7] == Approx(0.2) );
    }

    l_pss.apply( 0,
                 2,
                 0.2,
                 1.0,
                 l_dofs );
  }

  // reset element characteristics
  for( unsigned short l_el = 0; l_el < 6; l_el++ ) l_elChars[l_el].spType = 0;
  l_h5File = "cont/unit_tests/elastic/sources/point_sources_2d_1.h5";
  // return silently if file does not exist
  if( std::ifstream(l_h5File) ) {
    std::string l_h5Files[3];
    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ )
      l_h5Files[l_cr] = l_h5File;

    edge::data::Dynamic l_dynMem;

    edge::seismic::setups::PointSources<
     float,
     int,
     QUAD4R,
     2,
     3 > l_pss;
    l_pss.init( l_h5Files,
                l_gIdsEl,
                99,
                l_elLayout,
                l_elVe,
                l_veChars,
                l_massI,
                l_elChars,
                l_dynMem );

    /*
     * check the first example, which defines two point forces.
     */
    REQUIRE( l_elChars[0].spType ==  0 );
    REQUIRE( l_elChars[1].spType == 99 );
    REQUIRE( l_elChars[2].spType ==  0 );
    REQUIRE( l_elChars[3].spType ==  0 );
    REQUIRE( l_elChars[4].spType == 99 );
    REQUIRE( l_elChars[5].spType == 99 );
    REQUIRE( l_elChars[6].spType == 99 );

    for( unsigned short l_cr = 0; l_cr < 3; l_cr++ ) {
      REQUIRE( l_pss.m_elSpPs[0][l_cr] == 0 ); // src 0 + 3
      REQUIRE( l_pss.m_elSpPs[1][l_cr] == 2 ); // src 2
      REQUIRE( l_pss.m_elSpPs[2][l_cr] == 3 ); // src 1
      REQUIRE( l_pss.m_elSpPs[3][l_cr] == 4 ); // src 1 (dupl)
      REQUIRE( l_pss.m_elSpPs[4][l_cr] == 5 ); // ghost

     /*
      * global source 0 is in element 1
      * global source 1 is in element 5 and 6
      * global source 2 is in element 4
      * global source 3 is in element 1
      */
      REQUIRE( l_pss.m_psFused[l_cr].nPts == 5 );

      REQUIRE( l_pss.m_psFused[l_cr].el[0] == 1 );
      REQUIRE( l_pss.m_psFused[l_cr].el[1] == 1 );
      REQUIRE( l_pss.m_psFused[l_cr].el[2] == 4 );
      REQUIRE( l_pss.m_psFused[l_cr].el[3] == 5 );
      REQUIRE( l_pss.m_psFused[l_cr].el[4] == 6 );

      REQUIRE( l_pss.m_psFused[l_cr].times[0] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].times[1] == Approx(0.4) );
      REQUIRE( l_pss.m_psFused[l_cr].times[2] == Approx(0.2) );
      REQUIRE( l_pss.m_psFused[l_cr].times[3] == Approx(2.3) );
      REQUIRE( l_pss.m_psFused[l_cr].times[4] == Approx(2.3) );

      REQUIRE( l_pss.m_psFused[l_cr].dts[0] == Approx(0.001) );
      REQUIRE( l_pss.m_psFused[l_cr].dts[1] == Approx(0.001) );
      REQUIRE( l_pss.m_psFused[l_cr].dts[2] == Approx(0.001) );
      REQUIRE( l_pss.m_psFused[l_cr].dts[3] == Approx(0.002) );
      REQUIRE( l_pss.m_psFused[l_cr].dts[4] == Approx(0.002) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[0][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][3] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[0][4] == Approx(2.3) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[1][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][1] == Approx(2.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][3] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[1][4] == Approx(3.3) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[2][0] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[2][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[2][2] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[2][3] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[2][4] == Approx(2.3) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[3][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[3][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[3][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[3][3] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[3][4] == Approx(3.3) );

      REQUIRE( l_pss.m_psFused[l_cr].scas[4][0] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[4][1] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[4][2] == 0 );
      REQUIRE( l_pss.m_psFused[l_cr].scas[4][3] == Approx(1.0) );
      REQUIRE( l_pss.m_psFused[l_cr].scas[4][4] == Approx(3.3) );
 
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[0] ==  0 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[1] ==  5 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[2] ==  8 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[3] == 12 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[4] == 15 );
      REQUIRE( l_pss.m_psFused[l_cr].tsPtrs[5] == 18 );
 
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 0] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 1] == Approx(0.1) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 2] == Approx(0.2) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 3] == Approx(0.1) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 4] == Approx(0.15) );

      REQUIRE( l_pss.m_psFused[l_cr].tss[ 5] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 6] == Approx(0.2) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 7] == Approx(0.5) );
  
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 8] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[ 9] == Approx(0.3) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[10] == Approx(0.4) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[11] == Approx(0.5) );
  
      REQUIRE( l_pss.m_psFused[l_cr].tss[12] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[13] == Approx(0.3) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[14] == Approx(0.2) );

      REQUIRE( l_pss.m_psFused[l_cr].tss[15] == Approx(0) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[16] == Approx(0.3) );
      REQUIRE( l_pss.m_psFused[l_cr].tss[17] == Approx(0.2) );
    }

    l_pss.apply( 0,
                 2,
                 0.2,
                 1.0,
                 l_dofs );
  }
}
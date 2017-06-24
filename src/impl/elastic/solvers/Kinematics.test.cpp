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
 * Unit tests for kinematic sources.
 **/

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64

#include <catch.hpp>
#include "Kinematics.hpp"
#include "impl/elastic/io/Nrf.h"

//: TODO: Remove this
#define private public
#include "impl/elastic/setups/KinematicsInit.hpp"
#undef private

TEST_CASE( "Kinematics: Apply dirac delta distribution.", "[kinematics][applyDirac]" ) {
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

  // velocity model
  t_bgPars l_vel[7];
  for( unsigned short l_el = 0; l_el < 7; l_el++ ) {
    l_vel[l_el].mu = l_el+1;
    l_vel[l_el].lam = l_el+2;
  }

  std::string l_ncFile = "cont/elastic/sources/unit_tests/kinematics_init_3.nc";

  // dynamic memory
  edge::data::Dynamic l_mem;

  if( std::ifstream(l_ncFile) ) {
    // NRF reader with 2-source buffer
    edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

    for( unsigned short l_kId = 0; l_kId < 4; l_kId++ ) l_nrf1.init( l_ncFile );

    // inverse mass matrix
    double l_massI[9] = {1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    // first pt-sources of the sparse source elements
    int (*l_kiSp)[4] = nullptr;;
    // source solvers
    edge::elastic::solvers::t_Kinematics< 2,
                                          9,
                                          1,
                                          double,
                                          int > l_kiSo[4];

    // get solvers
    edge::elastic::setups::KinematicsInit<
      int,
      double,
      double,
      QUAD4R,
      3,
      4,
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
                  l_kiSp,
                  l_kiSo );

    double l_dofs[7][5][9][4];
    for( unsigned short l_el = 0; l_el < 7; l_el++ )
      for( unsigned short l_qt = 0; l_qt < 5; l_qt++ )
        for( unsigned short l_md = 0; l_md < 9; l_md++ )
          for( unsigned short l_ru = 0; l_ru < 4; l_ru++ )
            l_dofs[l_el][l_qt][l_md][l_ru] = 0;

    // call the solver
    edge::elastic::solvers::Kinematics< QUAD4R,
                                        5,
                                        3,
                                        4,
                                        1 >::applyDirac( 1,
                                                         2,
                                                         1.0,
                                                         3.7,
                                                         l_kiSp,
                                                         l_kiSo, 
                                                         l_dofs );
  }
}

#endif

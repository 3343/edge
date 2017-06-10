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
 * Unit tests for output of receivers at quadrature points.
 **/

// TODO: unit tests only valid for double-precision arithmetic since dg::Basis is not templatized.
#if PP_PRECISION == 64

#include <catch.hpp>
#define private public
#define protected public
#include "ReceiversQuad.hpp"
#undef protected
#undef private
#include "dg/QuadratureEval.hpp"

TEST_CASE( "Receivers at quad points: Initialization, points on faces, tria3.", "[RecvsQuad][initFaTria3]" ) {
  /*
   * Sketch:
   *                                                         5  -- 2.5
   *                                                       **
   *                                                     * *
   *                                                   *  *
   *                 4                               *   *      --  2  
   *               *    *                          7    *
   *             0  r3     *                     a     2
   *           a              f                f      a
   *         f       el2        a            *       f
   *       *                       6       *   el1  *
   *     *                            *  *         * 
   *   2****r0********fa1********r1***3           *             --  0
   *    5                            *  *        *
   *      a        el0      fa4 *        f      *
   *        f              *      el3     a    *            r2
   *          *       *                    8  *
   *            0*************fa3***********1                   -- -1
   *
   *   |        |     |                |     |                |
   *  -2       -1.4  -1                0    0.5               3
   *
   *
   * #receivers: 4
   * order: 3
   * #vertices: 6
   * #faces: 9
   * #elements: 4
   *
   * sp faces: fa1, fa8
   *
   *   el | fa0 | fa1 | fa3 || fa | el0 | el1 || ve |  x  |  y  || el | ve0 | ve1 | ve2 || re |  x  |  y  |
   *    0 |   4 |   1 |   5 ||  0 |   2 | xxx ||  0 |-1.4 |-1.0 ||  0 |   0 |   3 |   2 ||  0 |-1.6 | 0.0 |
   *    1 |   8 |   2 |   7 ||  1 |   0 |   2 ||  1 | 0.5 |-1.0 ||  1 |   3 |   1 |   5 ||  1 |-0.3 | 0.0 |
   *    2 |   0 |   1 |   6 ||  2 | xxx |   1 ||  2 |-2.0 | 0.0 ||  2 |   4 |   2 |   3 ||  2 | 2.9 |-0.9 |
   *    3 |   4 |   3 |   8 ||  3 |   3 | xxx ||  3 | 0.0 | 0.0 ||  3 |   3 |   0 |   1 ||  3 |-1.1 | 1.3 |
   *                            4 |   3 |   0 ||  4 |-1.0 | 2.0 ||    |     |     |
   *                            5 | xxx |   0 ||  5 | 3.0 | 2.5 ||xx=4|   6 |   7 |   8
   *                            6 |   2 | xxx ||    |
   *                            7 |   1 | xxx ||  6 |  13 |  14
   *                            8 |   3 |   1 ||  7 |  15 |  16
   *                                              8 |  17 |  18
   *
   *  el | tg | own?
   *   0 |  0 | yes
   *   1 |  0 | yes
   *   2 |  0 |  no
   *   3 |  1 | yes
   *
   * Remark: Since the internal storage of the receivers follows the layout of the sparse
   *         entities, we store the receivers in the order r0, r1, r3, r2. Reason is
   *         fa1 < fa8.
   *
   */


  unsigned int l_nRecvs = 4;
  std::string l_outDir = "/tmp/recvs_fa_tria3";
  std::string l_names[4] = {  "recv_one",
                              "recv_2",
                              "third_recv",
                              "last_and_four" };

  double l_recvCrds[4][2] = { {-1.6,  0.0},
                              {-0.3,  0.0},
                              { 2.9, -0.9},
                              {-1.1, 1.3} };
  double l_freq = 0.3;

  t_enLayout l_elLayout;
  l_elLayout.timeGroups.resize( 2 );
  l_elLayout.timeGroups[0].nEntsOwn = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[1].nEntsOwn = 1;
  l_elLayout.timeGroups[1].nEntsNotOwn = 0;

  t_enLayout l_faLayout;
  l_faLayout.nEnts = 9;
  l_faLayout.timeGroups.resize( 2 );
  // all faces are adjacent to time group 0
  l_faLayout.timeGroups[0].nEntsOwn = 8;
  l_faLayout.timeGroups[0].nEntsNotOwn = 1; // fa8 is assumed to be recv-face
  l_faLayout.timeGroups[1].nEntsOwn = 0;
  l_faLayout.timeGroups[1].nEntsNotOwn = 0;

  int l_elFa[4][3] = { {4,1,5},
                       {8,2,7},
                       {0,1,6},
                       {4,3,8} };

  int l_faEl[9][2] = { {2,4},
                       {0,2},
                       {4,1},
                       {3,4},
                       {3,0},
                       {4,0},
                       {2,4},
                       {1,4},
                       {3,1} };

  int l_elVe[5][3] = { {0,3,2},
                       {3,1,5},
                       {4,2,3},
                       {3,0,1},
                       {6,7,8} };

  struct { double coords[2]; } l_veChars[9];
  l_veChars[0].coords[0] = -1.4; l_veChars[0].coords[1] = -1.0;
  l_veChars[1].coords[0] =  0.5; l_veChars[1].coords[1] = -1.0;
  l_veChars[2].coords[0] = -2.0; l_veChars[2].coords[1] =  0.0;
  l_veChars[3].coords[0] =  0.0; l_veChars[3].coords[1] =  0.0;
  l_veChars[4].coords[0] = -1.0; l_veChars[4].coords[1] =  2.0;
  l_veChars[5].coords[0] =  3.0; l_veChars[5].coords[1] =  2.5;
  l_veChars[6].coords[0] = 13.0; l_veChars[6].coords[1] = 14.0;
  l_veChars[7].coords[0] = 14.0; l_veChars[7].coords[1] = 16.0;
  l_veChars[8].coords[0] = 15.0; l_veChars[8].coords[1] = 18.0;

  struct { int spType; } l_faChars[9];
  l_faChars[0].spType = 1;
  l_faChars[1].spType = 1024;
  l_faChars[2].spType = 2;
  l_faChars[3].spType = 3;
  l_faChars[4].spType = 4;
  l_faChars[5].spType = 5;
  l_faChars[6].spType = 6;
  l_faChars[7].spType = 7;
  l_faChars[8].spType = 1024;


  // get the quad points
  double l_pts[6][3][2];
  double l_weights[3];
  double l_basisEval[6][3][6];
  edge::dg::QuadratureEval<TRIA3, 3, 1>::faces( l_pts, l_weights, l_basisEval );

  edge::io::ReceiversQuad< double, TRIA3, 3, 8 > l_recvs;
  l_recvs.init( l_nRecvs,
                1024,
                3,
                l_outDir,
                l_names,
                l_recvCrds,
                l_freq,
                l_pts,
                l_faLayout,
                l_elLayout,
                l_faEl,
                l_elVe,
                l_elFa,
                l_veChars,
                l_faChars,
                13 );

  // check the update face type
  REQUIRE( (l_faChars[1].spType & 1024) == 1024 );
  REQUIRE( (l_faChars[8].spType & 1024) == 1024 );

  REQUIRE( l_recvs.m_recvs.size() == 4 );
  REQUIRE( l_recvs.m_recvsQuad.size() == 4 );

  REQUIRE( l_recvs.m_buffSize == 13 );

  for( unsigned short l_re = 0; l_re < 4; l_re++ ) {
    REQUIRE( l_recvs.m_recvs[l_re].nBuff == 0 );
    REQUIRE( l_recvs.m_recvs[l_re].buffer.size() == 13*3*8 ); 
    REQUIRE( l_recvs.m_recvs[l_re].time == Approx(0.0) );
    REQUIRE( l_recvs.m_recvs[l_re].tg == 0 );
  }

  REQUIRE( l_recvs.m_recvs[0].tg == 0 );
  REQUIRE( l_recvs.m_recvs[1].tg == 0 );
  REQUIRE( l_recvs.m_recvs[2].tg == 0 );
  REQUIRE( l_recvs.m_recvs[3].tg == 0 );

  REQUIRE( l_recvs.m_recvs[0].en == 0 );
  REQUIRE( l_recvs.m_recvs[1].en == 0 );
  REQUIRE( l_recvs.m_recvs[2].en == 0 );
  REQUIRE( l_recvs.m_recvs[3].en == 1 );

  REQUIRE( l_recvs.m_recvs[0].enTg == 0 );
  REQUIRE( l_recvs.m_recvs[1].enTg == 0 );
  REQUIRE( l_recvs.m_recvs[2].enTg == 0 );
  REQUIRE( l_recvs.m_recvs[3].enTg == 1 );

  REQUIRE( l_recvs.m_recvs[0].path == "/tmp/recvs_fa_tria3/recv_one.csv" );
  REQUIRE( l_recvs.m_recvs[1].path == "/tmp/recvs_fa_tria3/recv_2.csv" );
  REQUIRE( l_recvs.m_recvs[2].path == "/tmp/recvs_fa_tria3/last_and_four.csv" );
  REQUIRE( l_recvs.m_recvs[3].path == "/tmp/recvs_fa_tria3/third_recv.csv" );

  double l_lineQps[3];
  l_lineQps[0] = ( 1.0 - std::sqrt( 3.0 / 5.0 ) ) / 2.0;
  l_lineQps[1] = 0.5;
  l_lineQps[2] = ( 1.0 + std::sqrt( 3.0 / 5.0 ) ) / 2.0;

  REQUIRE( l_recvs.m_recvsQuad[0].qp == 2 );
  REQUIRE( l_recvs.m_recvsQuad[1].qp == 0 );
  REQUIRE( l_recvs.m_recvsQuad[2].qp == 1 );
  REQUIRE( l_recvs.m_recvsQuad[3].qp == 0 );

  REQUIRE( l_recvs.m_recvsQuad[0].crds[0] == Approx(-2.0 + 2.0*l_lineQps[0]) ); 
  REQUIRE( l_recvs.m_recvsQuad[0].crds[1] == Approx( 0.0 ) );

  REQUIRE( l_recvs.m_recvsQuad[1].crds[0] == Approx( 0.0 - 2.0*l_lineQps[0]) );
  REQUIRE( l_recvs.m_recvsQuad[1].crds[1] == Approx( 0.0 ) );

  REQUIRE( l_recvs.m_recvsQuad[2].crds[0] == Approx(-2.0 + 2.0*l_lineQps[1]) ); 
  REQUIRE( l_recvs.m_recvsQuad[2].crds[1] == Approx( 0.0 ) );

  REQUIRE( l_recvs.m_recvsQuad[3].crds[0] == Approx( 0.5 * l_lineQps[2]) );
  REQUIRE( l_recvs.m_recvsQuad[3].crds[1] == Approx(-1.0 * l_lineQps[2]) );

  REQUIRE( l_recvs.m_spEnToRecv[0] == 0 );
  REQUIRE( l_recvs.m_spEnToRecv[1] == 3 );

  // write receivers
  double l_data[2][3][3][8];

  // init receiver data
  for( unsigned short l_qt = 0; l_qt < 3; l_qt++ ) {
    for( unsigned short l_qp = 0; l_qp < 3; l_qp++ ) {
      for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
        l_data[0][l_qt][l_qp][l_ru] = (l_qt+1) * 10000 + (l_qp+1) * 100 + l_ru+1;
        l_data[1][l_qt][l_qp][l_ru] = (l_qt+1) * 20000 + (l_qp+1) * 200 + l_ru+2;
      }
    }
  }

  // write time 0
  l_recvs.writeRecvAll( 0.0, 0.2, 0, l_data[0] );
  l_recvs.writeRecvAll( 0.0, 0.2, 1, l_data[1] );

  // check get time
  REQUIRE( l_recvs.getRecvTimeRel( 0, 0.0, 0.2 ) < 0 );

  // write time
  l_recvs.writeRecvAll( 0.9, 0.1, 0, l_data[1] );
  l_recvs.writeRecvAll( 0.9, 0.1, 1, l_data[0] );

  // check get time
  REQUIRE( l_recvs.getRecvTimeRel( 0, 1.0, 0.6 ) == Approx( 0.2 ) );
  REQUIRE( l_recvs.getRecvTimeRel( 1, 1.0, 0.6 ) == Approx( 0.2 ) );

  // write again
  l_recvs.writeRecvAll( 1.25, 0.15, 0, l_data[0] );
  l_recvs.writeRecvAll( 2.35, 0.15, 1, l_data[1] );

  // flush output
  l_recvs.flushIf( 0 );
}

#endif

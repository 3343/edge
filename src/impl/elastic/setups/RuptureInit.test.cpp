/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Tests the initialization of rupture physics.
 **/

#include <catch.hpp>
#include "RuptureInit.hpp"
#include "linalg/Domain.hpp"
#include "linalg/HalfSpace.hpp"

TEST_CASE( "Init of rupture physics.", "[RuptureInit][LSW2D]" ) {
  /*
   * 2D setup:
   *
   * Remark: This is a "coarse" illustration only,
   *         the setup is mostly random and inconsistent.
   *
   *
   *            |                  [...]                   |  [...]
   *                 |          nucleation 2         | 30,000
   *                     |      nucleation 1     | 20,000
   *                         |  nucleation 0  | 10,000
   *    *****************************+*****************************
   *    *             f0             +     f1                     *
   *    *            v1              +        v2 el3              *
   *    *        *  *            plus side   el2 *                *
   *    *     el3  x                 +              x             *
   *    *         * el1              +            *    *          *
   *    *        v0                  +                    v3      *
   *    *      ------------------------------------------ "fault" * 
   *    *                   el0  *   - (0,0)                 v7   *
   *    *                 f2 v4****v5-                       x    *
   *    *                   el2  minus side             el0 *   * *
   *    *                            -                     * el3  *
   *    *                            -                   v6       *
   *    *                            -                  f3        *
   *    x****************************-****************************x
   *    (-100,000, -25,000)                      (100,000, -25,000)
   *
   *
   *      normal | in nucleation patch | left | right | minus | plus | vertices (don't match normals)       |
   * f0: (-1, 1) | > 2                 | el1  | el3   | el1   | el3  | (-42,000,  10,000) (-31,000, 30,000) |
   * f1  (-1,-1) | > 1                 | el3  | el2   | el2   | el3  | ( 22,000,   5,000) ( 23,000,  9,000) |
   * f2: ( 0, 1) | not a rupture face  | el2  | el0   | ---   | ---  | (  -3000,   -2000) (  1,000, -2,000) |
   * f3: ( 1,-1) | > 3                 | el0  | el3   | el3   | el0  | ( 43,000, -33,000) ( 47,000,   -100) |
   *
   * friction parameters:
   *   mus:   run + 0.1
   *   mud:   run * 2 + 0.01
   *   dcInv: run * 3 + 0.3
   *
   * nucleation patches:
   *   sn0: run * 1.0E6
   *   ss0: run * -2.0E6
   *
   * remainder of the fault:
   *   sn0: -3.0E6
   *   ss0:  4.0E6
   *
   * element-wise lame parameters:
   *   lam: 20 + el
   *   mu:  30 + el
   *   rho: 40 + el
   *
   *   therefore:
   *     element | csDmu
   *     0       | 0.02886751346
   *     1       | 0.02804963567
   *     2       | 0.02727723628
   *     3       | 0.02654659366
   */
  // set friction parameters
  double l_lswPars[8][3];
  for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
    l_lswPars[l_ru][0] = l_ru     + 0.1;
    l_lswPars[l_ru][1] = l_ru * 2 + 0.01;
    l_lswPars[l_ru][2] = l_ru * 3 + 0.3;
  }

  // create domains
  std::vector< edge::linalg::Domain< double, 2, edge::linalg::HalfSpace > > l_doms[8];

  for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
    l_doms[l_ru].resize(2);

    double l_origin[2], l_normal[2];

    // add halfspace to dom 0
    l_origin[0] = -10000 - 10000*l_ru;
    l_origin[1] =     0;
    l_normal[0] =     1;
    l_normal[1] =     0;
    edge::linalg::HalfSpace<double, 2> l_hs( l_origin, l_normal );
    l_doms[l_ru][0].add( l_hs );

    // add second half space to dom 0
    l_origin[0] =  10000 + 10000*l_ru;
    l_origin[1] =     0;
    l_normal[0] =    -1;
    l_normal[1] =     0;
    l_hs.setOrigin( l_origin );
    l_hs.setNormal( l_normal );
    l_doms[l_ru][0].add( l_hs );

    // add half space to dom 1
    l_origin[0] =  100000;
    l_origin[1] =       0;
    l_normal[0] =      -1;
    l_normal[1] =       0;
    l_hs.setOrigin( l_origin );
    l_hs.setNormal( l_normal );
    l_doms[l_ru][1].add( l_hs );
  }

  // set up stresses
  std::vector< std::array< double, 2 > > l_stress[8];

  for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
    l_stress[l_ru].resize(2);

    // "nucleation"
    l_stress[l_ru][0][0] = l_ru *  1.0E6;
    l_stress[l_ru][0][1] = l_ru * -2.0E6;

    // "fault"
    l_stress[l_ru][1][0] = l_ru * -3.0E6;
    l_stress[l_ru][1][1] = l_ru *  4.0E6;
  }

  // set up element background paramters
  t_bgPars l_bgPars[4];
  for( unsigned short l_el = 0; l_el < 4; l_el++ ) {
    l_bgPars[l_el].lam = 20 + l_el;
    l_bgPars[l_el].mu  = 30 + l_el;
    l_bgPars[l_el].rho = 40 + l_el;
  }

  int l_nFaces = 4;
  int l_faVe[4][2] = { {0, 1},
                       {2, 3},
                       {4, 5},
                       {6, 7} };
  int l_faEl[4][2] = { {1, 3},
                       {3, 2},
                       {2, 0},
                       {0, 3} };

  // set up coordinates of vertices
  t_vertexChars l_veChars[8];
  l_veChars[0].coords[0] = -42000;
  l_veChars[0].coords[1] =  10000;
  l_veChars[1].coords[0] = -31000;
  l_veChars[1].coords[1] =  30000;
  l_veChars[2].coords[0] =  22000;
  l_veChars[2].coords[1] =   5000;
  l_veChars[3].coords[0] =  23000;
  l_veChars[3].coords[1] =   9000;
  l_veChars[4].coords[0] =  -3000;
  l_veChars[4].coords[1] =  -2000;
  l_veChars[5].coords[0] =   1000;
  l_veChars[5].coords[1] =  -2000;
  l_veChars[6].coords[0] =  43000;
  l_veChars[6].coords[1] = -33000;
  l_veChars[7].coords[0] =  47000;
  l_veChars[7].coords[1] =   -100;

  // set up outer pointing normals of the faces
  t_faceChars l_faChars[4];
  l_faChars[0].spType = 1234;
  l_faChars[0].outNormal[0] = -1;
  l_faChars[0].outNormal[1] =  1;

  l_faChars[1].spType = 1234;
  l_faChars[1].outNormal[0] = -1;
  l_faChars[1].outNormal[1] = -1;

  l_faChars[2].spType = 4321;
  l_faChars[2].outNormal[0] =  0;
  l_faChars[2].outNormal[1] =  1;

  l_faChars[3].spType = 1234;
  l_faChars[3].outNormal[0] =  1;
  l_faChars[3].outNormal[1] = -1;

  double l_faultCrds[2][2] = { {0, 1}, {1, 0} };

  // output
  edge::elastic::solvers::t_LinSlipWeakGlobal<double, 8>     l_lswGlobal;
  edge::elastic::solvers::t_LinSlipWeakFace<double>          l_lswFace[3];
  edge::elastic::solvers::t_LinSlipWeakSubFace<double, 2, 8> l_lswSf[3][5];
  edge::elastic::solvers::t_LinSlipWeak< double, TRIA3, 3, 8 > l_lsw;
  l_lsw.gl = l_lswGlobal;
  l_lsw.fa = l_lswFace;
  l_lsw.sf = l_lswSf;

  // init the rupture physics
  edge::elastic::setups::RuptureInit< TRIA3, 3, 8>::linSlipWeak( l_nFaces,
                                                                 1234,
                                                                 l_faVe,
                                                                 l_faEl,
                                                                 l_veChars,
                                                                 l_faChars,
                                                                 l_bgPars,
                                                                 l_faultCrds,
                                                                 l_lswPars,
                                                                 l_doms,
                                                                 l_stress,
                                                                 l_lsw );

  // check the resulting friction parameters
  for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
    REQUIRE( l_lswGlobal.mus[l_ru]   == Approx( l_ru + 0.1 ) );
    REQUIRE( l_lswGlobal.mud[l_ru]   == Approx( l_ru * 2 + 0.01 ) );
    REQUIRE( l_lswGlobal.dcInv[l_ru] == Approx( 1.0 / (l_ru * 3 + 0.3) ) );
  }

  // check the stress setup
  for( unsigned short l_ru = 0; l_ru < 8; l_ru++ ) {
    for( unsigned short l_sp = 0; l_sp < 3; l_sp++ ) {
      for( unsigned short l_sf = 0; l_sf < 5; l_sf++ ) {
        if(    (l_sp == 0 && l_ru > 2)
            || (l_sp == 1 && l_ru > 1)
            || (l_sp == 2 && l_ru > 3) ) {
          REQUIRE( l_lswSf[l_sp][l_sf].sn0[l_ru]    == Approx(l_ru *  1.0E6) );
          REQUIRE( l_lswSf[l_sp][l_sf].ss0[0][l_ru] == Approx(l_ru * -2.0E6) );
          REQUIRE( l_lswSf[l_sp][l_sf].muf[l_ru]    == l_lswGlobal.mus[l_ru] );
        }
        else {
          REQUIRE( l_lswSf[l_sp][l_sf].sn0[l_ru]    == Approx(l_ru * -3.0E6) );
          REQUIRE( l_lswSf[l_sp][l_sf].ss0[0][l_ru] == Approx(l_ru *  4.0E6) );
          REQUIRE( l_lswSf[l_sp][l_sf].muf[l_ru]    == l_lswGlobal.mus[l_ru] );
        }
      }
    }
  }

  // check plus/minus and csDmu
  REQUIRE( l_lswFace[0].lEqM == false );
  REQUIRE( l_lswFace[1].lEqM == true  );
  REQUIRE( l_lswFace[2].lEqM == true  );

  REQUIRE( l_lswFace[0].csDmuM == Approx( 0.02654659366 ) );
  REQUIRE( l_lswFace[0].csDmuP == Approx( 0.02804963567 ) );

  REQUIRE( l_lswFace[1].csDmuM == Approx( 0.02654659366 ) );
  REQUIRE( l_lswFace[1].csDmuP == Approx( 0.02727723628 ) );

  REQUIRE( l_lswFace[2].csDmuM == Approx( 0.02886751346 ) );
  REQUIRE( l_lswFace[2].csDmuP == Approx( 0.02654659366 ) );
}

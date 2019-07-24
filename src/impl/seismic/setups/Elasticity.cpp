/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2019, Regents of the University of California
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
 * Elasticity.
 **/
#include "Elasticity.h"

void edge::seismic::setups::Elasticity::star( double       i_rho,
                                              double       i_lam,
                                              double       i_mu, 
                                              double const i_jacInv[2][2],
                                              double       o_starE[2][5 * 5] ) {
  /*
  * Jacobians in q_t = A(x,y) * q_x + B(x,y) * q_y
  *
  * A:
  *    _____0__1_______2_______________3____4
  *  0|     0, 0,      0, -lambda - 2*mu,   0|0
  *  1|     0, 0,      0,        -lambda,   0|1
  *  2|     0, 0,      0,              0, -mu|2
  *  3|-1/rho, 0,      0,              0,   0|3
  *  4|     0, 0, -1/rho,              0,   0|4
  *    -----0--1-------2---------------3----4
  *
  * B:
  *    0_______1_______2____3_______________4
  *  0|0,      0,      0,   0,        -lambda|0
  *  1|0,      0,      0,   0, -lambda - 2*mu|1
  *  2|0,      0,      0, -mu,              0|2
  *  3|0,      0, -1/rho,   0,              0|3
  *  4|0, -1/rho,      0,   0,              0|4
  *    0-------1-------2----3---------------4
  */
  // init to zero
  for( unsigned short l_di = 0; l_di < 2; l_di++ )
    for( unsigned short l_m0 = 0; l_m0 < 5; l_m0++ )
      for( unsigned short l_m1 = 0; l_m1 < 5; l_m1++ )
        o_starE[l_di][l_m0*5 + l_m1] = 0;

  // iterate over reference dimension
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
    // set non-zeros of first jacobian
    o_starE[l_di][0*5 + 3] += (-i_lam - 2*i_mu) * i_jacInv[0][l_di];
    o_starE[l_di][1*5 + 3] += (-i_lam         ) * i_jacInv[0][l_di];
    o_starE[l_di][2*5 + 4] += (-i_mu          ) * i_jacInv[0][l_di];
    o_starE[l_di][3*5 + 0] += (-1/i_rho       ) * i_jacInv[0][l_di];
    o_starE[l_di][4*5 + 2] += (-1/i_rho       ) * i_jacInv[0][l_di];

    // set non-zeros of second jacobian
    o_starE[l_di][0*5 + 4] += (-i_lam         ) * i_jacInv[1][l_di];
    o_starE[l_di][1*5 + 4] += (-i_lam - 2*i_mu) * i_jacInv[1][l_di];
    o_starE[l_di][2*5 + 3] += (-i_mu          ) * i_jacInv[1][l_di];
    o_starE[l_di][3*5 + 2] += (-1/i_rho       ) * i_jacInv[1][l_di];
    o_starE[l_di][4*5 + 1] += (-1/i_rho       ) * i_jacInv[1][l_di];
  }
}

void edge::seismic::setups::Elasticity::star( double       i_rho,
                                              double       i_lam,
                                              double       i_mu, 
                                              double const i_jacInv[2][2],
                                              double       o_starE[2][10] ) {
  double l_starE[2][5*5];
  star( i_rho,
        i_lam,
        i_mu,
        i_jacInv,
        l_starE );

  // iterate over reference dimension
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
    o_starE[l_di][0] = l_starE[l_di][0*5 + 3];
    o_starE[l_di][1] = l_starE[l_di][0*5 + 4];

    o_starE[l_di][2] = l_starE[l_di][1*5 + 3];
    o_starE[l_di][3] = l_starE[l_di][1*5 + 4];

    o_starE[l_di][4] = l_starE[l_di][2*5 + 3];
    o_starE[l_di][5] = l_starE[l_di][2*5 + 4];

    o_starE[l_di][6] = l_starE[l_di][3*5 + 0];
    o_starE[l_di][7] = l_starE[l_di][3*5 + 2];

    o_starE[l_di][8] = l_starE[l_di][4*5 + 1];
    o_starE[l_di][9] = l_starE[l_di][4*5 + 2];
  }
}

void edge::seismic::setups::Elasticity::star( double       i_rho,
                                              double       i_lam,
                                              double       i_mu, 
                                              double const i_jacInv[3][3],
                                              double       o_starE[3][9 * 9] ) {
  /*
    * Jacobians in q_t = A(x,y,z) * q_x + B(x,y,z) * q_y + C(x,y,z) * q_z
    *
    * A:
    *   _____0__1__2_______3__4_______5_______________6____7____8
    * 0|     0, 0, 0,      0  0,      0, -lambda - 2*mu,   0,   0|0
    * 1|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|1
    * 2|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|2
    * 3|     0, 0, 0,      0, 0,      0,              0, -mu,   0|3
    * 4|     0, 0, 0,      0, 0,      0,              0,   0,   0|4
    * 5|     0, 0, 0,      0, 0,      0,              0,   0, -mu|5
    * 6|-1/rho, 0, 0,      0, 0,      0,              0,   0,   0|6
    * 7|     0, 0, 0, -1/rho, 0,      0,              0,   0,   0|7
    * 8|     0, 0, 0,      0, 0, -1/rho,              0,   0,   0|8
    *   -----0--1--2-------3--4-------5---------------6----7----8
    *
    * B:
    *   0_______1__2_______3_______4__5____6_________________7______8
    * 0|0,      0, 0,      0       0, 0,   0,          -lambda,     0|0
    * 1|0,      0, 0,      0,      0, 0,   0,   -lambda - 2*mu,     0|1
    * 2|0,      0, 0,      0,      0, 0,   0,          -lambda,     0|2
    * 3|0,      0, 0,      0,      0, 0, -mu,                0,     0|3
    * 4|0,      0, 0,      0,      0, 0,   0,                0,   -mu|4
    * 5|0,      0, 0,      0,      0, 0,   0,                0,     0|5
    * 6|0,      0, 0, -1/rho,      0, 0,   0,                0,     0|6
    * 7|0, -1/rho, 0,      0,      0, 0,   0,                0,     0|7
    * 8|0,      0, 0,      0, -1/rho, 0,   0,                0,     0|8
    *   0-------1--2-------3-------4--5---------------6------7------8
    *
    * C:
    *   0__1_______2__3_______4_______5____6____7_______________8
    * 0|0, 0,      0, 0       0,      0,   0,   0,        -lambda|0
    * 1|0, 0,      0, 0,      0,      0,   0,   0,        -lambda|1
    * 2|0, 0,      0, 0,      0,      0,   0,   0, -lambda - 2*mu|2
    * 3|0, 0,      0, 0,      0,      0,   0,   0,              0|3
    * 4|0, 0,      0, 0,      0,      0,   0, -mu,              0|4
    * 5|0, 0,      0, 0,      0,      0, -mu,   0,              0|5
    * 6|0, 0,      0, 0,      0, -1/rho,   0,   0,              0|6
    * 7|0, 0,      0, 0, -1/rho,      0,   0,   0,              0|7
    * 8|0, 0, -1/rho, 0,      0,      0,   0,   0,              0|8
    *   0--1-------2--3-------4-------5----6----7---------------8
    */
  // init to zero
  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    for( unsigned short l_m0 = 0; l_m0 < 9; l_m0++ )
      for( unsigned short l_m1 = 0; l_m1 < 9; l_m1++ )
        o_starE[l_di][l_m0*9 + l_m1] = 0;

  // iterate over reference dimension
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    // set non-zeros of first jacobian
    o_starE[l_di][0*9 + 6] += (-i_lam - 2*i_mu) * i_jacInv[0][l_di];
    o_starE[l_di][1*9 + 6] += (-i_lam         ) * i_jacInv[0][l_di];
    o_starE[l_di][2*9 + 6] += (-i_lam         ) * i_jacInv[0][l_di];
    o_starE[l_di][3*9 + 7] += (-i_mu          ) * i_jacInv[0][l_di];
    o_starE[l_di][5*9 + 8] += (-i_mu          ) * i_jacInv[0][l_di];
    o_starE[l_di][6*9 + 0] += (-1 / i_rho     ) * i_jacInv[0][l_di];
    o_starE[l_di][7*9 + 3] += (-1 / i_rho     ) * i_jacInv[0][l_di];
    o_starE[l_di][8*9 + 5] += (-1 / i_rho     ) * i_jacInv[0][l_di];

    // set non-zeros of second jacobian
    o_starE[l_di][0*9 + 7] += (-i_lam         ) * i_jacInv[1][l_di];
    o_starE[l_di][1*9 + 7] += (-i_lam - 2*i_mu) * i_jacInv[1][l_di];
    o_starE[l_di][2*9 + 7] += (-i_lam         ) * i_jacInv[1][l_di];
    o_starE[l_di][3*9 + 6] += (-i_mu          ) * i_jacInv[1][l_di];
    o_starE[l_di][4*9 + 8] += (-i_mu          ) * i_jacInv[1][l_di];
    o_starE[l_di][6*9 + 3] += (-1 / i_rho     ) * i_jacInv[1][l_di];
    o_starE[l_di][7*9 + 1] += (-1 / i_rho     ) * i_jacInv[1][l_di];
    o_starE[l_di][8*9 + 4] += (-1 / i_rho     ) * i_jacInv[1][l_di];

    // set non-zeros of third jacobian
    o_starE[l_di][0*9 + 8] += (-i_lam         ) * i_jacInv[2][l_di];
    o_starE[l_di][1*9 + 8] += (-i_lam         ) * i_jacInv[2][l_di];
    o_starE[l_di][2*9 + 8] += (-i_lam - 2*i_mu) * i_jacInv[2][l_di];
    o_starE[l_di][4*9 + 7] += (-i_mu          ) * i_jacInv[2][l_di];
    o_starE[l_di][5*9 + 6] += (-i_mu          ) * i_jacInv[2][l_di];
    o_starE[l_di][6*9 + 5] += (-1 / i_rho     ) * i_jacInv[2][l_di];
    o_starE[l_di][7*9 + 4] += (-1 / i_rho     ) * i_jacInv[2][l_di];
    o_starE[l_di][8*9 + 2] += (-1 / i_rho     ) * i_jacInv[2][l_di];
  }
}

void edge::seismic::setups::Elasticity::star( double       i_rho,
                                              double       i_lam,
                                              double       i_mu, 
                                              double const i_jacInv[3][3],
                                              double       o_starE[3][24] ) {
  double l_starE[3][9*9];
  star( i_rho,
        i_lam,
        i_mu,
        i_jacInv,
        l_starE );

  // iterate over reference dimension
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    o_starE[l_di][ 0] = l_starE[l_di][0*9 + 6];
    o_starE[l_di][ 1] = l_starE[l_di][0*9 + 7];
    o_starE[l_di][ 2] = l_starE[l_di][0*9 + 8];

    o_starE[l_di][ 3] = l_starE[l_di][1*9 + 6];
    o_starE[l_di][ 4] = l_starE[l_di][1*9 + 7];
    o_starE[l_di][ 5] = l_starE[l_di][1*9 + 8];

    o_starE[l_di][ 6] = l_starE[l_di][2*9 + 6];
    o_starE[l_di][ 7] = l_starE[l_di][2*9 + 7];
    o_starE[l_di][ 8] = l_starE[l_di][2*9 + 8];

    o_starE[l_di][ 9] = l_starE[l_di][3*9 + 6];
    o_starE[l_di][10] = l_starE[l_di][3*9 + 7];

    o_starE[l_di][11] = l_starE[l_di][4*9 + 7];
    o_starE[l_di][12] = l_starE[l_di][4*9 + 8];

    o_starE[l_di][13] = l_starE[l_di][5*9 + 6];
    o_starE[l_di][14] = l_starE[l_di][5*9 + 8];

    o_starE[l_di][15] = l_starE[l_di][6*9 + 0];
    o_starE[l_di][16] = l_starE[l_di][6*9 + 3];
    o_starE[l_di][17] = l_starE[l_di][6*9 + 5];

    o_starE[l_di][18] = l_starE[l_di][7*9 + 1];
    o_starE[l_di][19] = l_starE[l_di][7*9 + 3];
    o_starE[l_di][20] = l_starE[l_di][7*9 + 4];

    o_starE[l_di][21] = l_starE[l_di][8*9 + 2];
    o_starE[l_di][22] = l_starE[l_di][8*9 + 4];
    o_starE[l_di][23] = l_starE[l_di][8*9 + 5];
  }
}
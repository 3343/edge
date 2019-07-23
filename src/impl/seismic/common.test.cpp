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
 * Unit tests for common function for the elastic wave equations.
 **/

#include <catch.hpp>
#include "common.hpp"

TEST_CASE( "Trafo3D", "[elasticCommon][trafo3D]" ) {
  // remote coordinate system
  double l_nx, l_ny, l_nz;
  double l_sx, l_sy, l_sz;
  double l_tx, l_ty, l_tz;

  double l_tI[9][9];

  l_nx = std::sqrt(2.0); l_ny = -std::sqrt(2.0); l_nz = 0.0;
  l_sx = std::sqrt(2.0); l_sy = std::sqrt(2.0); l_sz = 0.0;
  l_tx = 0; l_ty = 0; l_tz = 1;

  edge::seismic::common::setupTrafo3d( l_nx, l_ny, l_nz,
                                       l_sx, l_sy, l_sz,
                                       l_tx, l_ty, l_tz,
                                       l_tI );

  REQUIRE( l_tI[0][0] == Approx(  2.0 ) );
  REQUIRE( l_tI[0][1] == Approx(  2.0 ) );
  REQUIRE( l_tI[0][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[0][3] == Approx(  4.0 ) );
  REQUIRE( l_tI[0][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[0][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[0][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[0][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[0][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[1][0] == Approx(  2.0 ) );
  REQUIRE( l_tI[1][1] == Approx(  2.0 ) );
  REQUIRE( l_tI[1][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[1][3] == Approx( -4.0 ) );
  REQUIRE( l_tI[1][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[1][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[1][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[1][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[1][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[2][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][2] == Approx(  1.0 ) );
  REQUIRE( l_tI[2][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[2][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[3][0] == Approx( -2.0 ) );
  REQUIRE( l_tI[3][1] == Approx(  2.0 ) );
  REQUIRE( l_tI[3][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[3][8] == Approx(  0.0 ) );

  double l_s = std::sqrt(2.0);

  REQUIRE( l_tI[4][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][4] == Approx(  l_s ) );
  REQUIRE( l_tI[4][5] == Approx( -l_s ) );
  REQUIRE( l_tI[4][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[4][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[5][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][4] == Approx(  l_s ) );
  REQUIRE( l_tI[5][5] == Approx(  l_s ) );
  REQUIRE( l_tI[5][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[5][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[6][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[6][6] == Approx(  l_s ) );
  REQUIRE( l_tI[6][7] == Approx(  l_s ) );
  REQUIRE( l_tI[6][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[7][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[7][6] == Approx( -l_s ) );
  REQUIRE( l_tI[7][7] == Approx(  l_s ) );
  REQUIRE( l_tI[7][8] == Approx(  0.0 ) );

  REQUIRE( l_tI[8][0] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][1] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][2] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][3] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][4] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][5] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][6] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][7] == Approx(  0.0 ) );
  REQUIRE( l_tI[8][8] == Approx(  1.0 ) );
}

TEST_CASE( "TrafoInv3D", "[elasticCommon][trafoInv3D]" ) {
  // remote coordinate system
  double l_nx, l_ny, l_nz;
  double l_sx, l_sy, l_sz;
  double l_tx, l_ty, l_tz;

  double l_tI[9][9];

  l_nx = std::sqrt(2.0); l_ny = -std::sqrt(2.0); l_nz = 0.0;
  l_sx = std::sqrt(2.0); l_sy = std::sqrt(2.0); l_sz = 0.0;
  l_tx = 0; l_ty = 0; l_tz = 1;

  edge::seismic::common::setupTrafoInv3d( l_nx, l_ny, l_nz,
                                          l_sx, l_sy, l_sz,
                                          l_tx, l_ty, l_tz,
                                          l_tI );

  REQUIRE( l_tI[0][0] == Approx(  0.125 ) );
  REQUIRE( l_tI[0][1] == Approx(  0.125 ) );
  REQUIRE( l_tI[0][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[0][3] == Approx( -0.25  ) );
  REQUIRE( l_tI[0][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[0][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[0][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[0][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[0][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[1][0] == Approx(  0.125 ) );
  REQUIRE( l_tI[1][1] == Approx(  0.125 ) );
  REQUIRE( l_tI[1][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[1][3] == Approx(  0.25  ) );
  REQUIRE( l_tI[1][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[1][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[1][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[1][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[1][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[2][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][2] == Approx(  1.0   ) );
  REQUIRE( l_tI[2][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[2][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[3][0] == Approx(  0.125 ) );
  REQUIRE( l_tI[3][1] == Approx( -0.125 ) );
  REQUIRE( l_tI[3][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[3][8] == Approx(  0.0   ) );

  double l_s = 1.0 / ( 2.0 * std::sqrt(2.0) );

  REQUIRE( l_tI[4][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][4] == Approx(  l_s   ) );
  REQUIRE( l_tI[4][5] == Approx(  l_s   ) );
  REQUIRE( l_tI[4][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[4][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[5][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][4] == Approx( -l_s   ) );
  REQUIRE( l_tI[5][5] == Approx(  l_s   ) );
  REQUIRE( l_tI[5][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[5][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[6][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[6][6] == Approx(  l_s   ) );
  REQUIRE( l_tI[6][7] == Approx( -l_s   ) );
  REQUIRE( l_tI[6][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[7][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[7][6] == Approx(  l_s   ) );
  REQUIRE( l_tI[7][7] == Approx(  l_s   ) );
  REQUIRE( l_tI[7][8] == Approx(  0.0   ) );

  REQUIRE( l_tI[8][0] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][1] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][2] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][3] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][4] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][5] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][6] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][7] == Approx(  0.0   ) );
  REQUIRE( l_tI[8][8] == Approx(  1.0   ) );
}

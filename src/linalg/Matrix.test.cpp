/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONsTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Unit tests of the matrix module.
 **/
#include <catch.hpp>

#include "Matrix.h"

TEST_CASE( "Matrix: Tests the dense to coordinate format conversion", "[matrix][denseToCrd]" ) {
  double l_mat[5][9];

  for( unsigned short l_ro = 0; l_ro < 5; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
      l_mat[l_ro][l_co] = 0;
    }
  }

  t_matCrd l_res;

  edge::linalg::Matrix::denseToCrd( 5, 9, l_mat[0], l_res );
  REQUIRE( l_res.nz.size() == 0 );
  REQUIRE( l_res.ro.size() == 0 );
  REQUIRE( l_res.co.size() == 0 );

  l_mat[0][3] = 15.0;
  l_mat[0][7] = -3.0;
  l_mat[2][4] = -9.3;
  l_mat[3][0] = 13.2;
  edge::linalg::Matrix::denseToCrd( 5, 9, l_mat[0], l_res );
  REQUIRE( l_res.nz.size() == 4 );
  REQUIRE( l_res.ro.size() == 4 );
  REQUIRE( l_res.co.size() == 4 );

  REQUIRE( l_res.nz[0] == (real_base) 15.0 );
  REQUIRE( l_res.nz[1] == (real_base) -3.0 );
  REQUIRE( l_res.nz[2] == (real_base) -9.3 );
  REQUIRE( l_res.nz[3] == (real_base) 13.2 );

  REQUIRE( l_res.ro[0] == 0 );
  REQUIRE( l_res.ro[1] == 0 );
  REQUIRE( l_res.ro[2] == 2 );
  REQUIRE( l_res.ro[3] == 3 );

  REQUIRE( l_res.co[0] == 3 );
  REQUIRE( l_res.co[1] == 7 );
  REQUIRE( l_res.co[2] == 4 );
  REQUIRE( l_res.co[3] == 0 );

}

TEST_CASE( "Matrix: Tests the dense to compressed sparse row format conversion", "[matrix][denseToCsr]" ) {
  double l_mat[5][9];

  for( unsigned short l_ro = 0; l_ro < 5; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
      l_mat[l_ro][l_co] = 0;
    }
  }

  t_matCsr l_res;

  edge::linalg::Matrix::denseToCsr( 5, 9, l_mat[0], l_res );
  REQUIRE( l_res.val.size()    == 0 );
  REQUIRE( l_res.colIdx.size() == 0 );
  REQUIRE( l_res.rowPtr.size() == 6 );


  l_mat[0][3] = 15.0;
  l_mat[0][7] = -3.0;
  l_mat[2][4] = -9.3;
  l_mat[3][0] = 13.2;
  edge::linalg::Matrix::denseToCsr( 5, 9, l_mat[0], l_res );
  REQUIRE( l_res.val.size()    == 4 );
  REQUIRE( l_res.colIdx.size() == 4 );
  REQUIRE( l_res.rowPtr.size() == 6 );

  REQUIRE( l_res.val[0] == (real_base) 15.0 );
  REQUIRE( l_res.val[1] == (real_base) -3.0 );
  REQUIRE( l_res.val[2] == (real_base) -9.3 );
  REQUIRE( l_res.val[3] == (real_base) 13.2 );

  REQUIRE( l_res.rowPtr[0] == 0 );
  REQUIRE( l_res.rowPtr[1] == 2 );
  REQUIRE( l_res.rowPtr[2] == 2 );
  REQUIRE( l_res.rowPtr[3] == 3 );
  REQUIRE( l_res.rowPtr[4] == 4 );
  REQUIRE( l_res.rowPtr[5] == 4 );

  REQUIRE( l_res.colIdx[0] == 3 );
  REQUIRE( l_res.colIdx[1] == 7 );
  REQUIRE( l_res.colIdx[2] == 4 );
  REQUIRE( l_res.colIdx[3] == 0 );
}

TEST_CASE( "Matrix: Tests the dense to compressed sparse column format conversion", "[matrix][denseToCsc]") {
  double l_mat[3][3];

  for( unsigned short l_ro = 0; l_ro < 3; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 3; l_co++ ) {
      l_mat[l_ro][l_co] = 0;
    }
  }

  t_matCsc l_res;

  edge::linalg::Matrix::denseToCsc( 3, 3, l_mat[0], l_res );
  REQUIRE( l_res.val.size()    == 0 );
  REQUIRE( l_res.rowIdx.size() == 0 );
  REQUIRE( l_res.colPtr.size() == 4 );


  l_mat[0][0] = 1.0;
  l_mat[0][2] = 2.0;
  l_mat[1][2] = 3.0;
  l_mat[2][0] = 4.0;
  l_mat[2][1] = 5.0;
  l_mat[2][2] = 6.0;
  edge::linalg::Matrix::denseToCsc( 3, 3, l_mat[0], l_res );
  REQUIRE( l_res.val.size()    == 6 );
  REQUIRE( l_res.rowIdx.size() == 6 );
  REQUIRE( l_res.colPtr.size() == 4 );

  REQUIRE( l_res.val[0] == (real_base) 1.0 );
  REQUIRE( l_res.val[1] == (real_base) 4.0 );
  REQUIRE( l_res.val[2] == (real_base) 5.0 );
  REQUIRE( l_res.val[3] == (real_base) 2.0 );
  REQUIRE( l_res.val[4] == (real_base) 3.0 );
  REQUIRE( l_res.val[5] == (real_base) 6.0 );

  REQUIRE( l_res.colPtr[0] == 0 );
  REQUIRE( l_res.colPtr[1] == 2 );
  REQUIRE( l_res.colPtr[2] == 3 );
  REQUIRE( l_res.colPtr[3] == 6 );

  REQUIRE( l_res.rowIdx[0] == 0 );
  REQUIRE( l_res.rowIdx[1] == 2 );
  REQUIRE( l_res.rowIdx[2] == 2 );
  REQUIRE( l_res.rowIdx[3] == 0 );
  REQUIRE( l_res.rowIdx[4] == 1 );
  REQUIRE( l_res.rowIdx[5] == 2 );
}

TEST_CASE( "Matrix: Derivation of non-zero blocks in matrices", "[matrix][getBlockNz]" ) {
  // dense test matrix
  real_base l_matDense[9][9];
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
      l_matDense[l_ro][l_co] = 0;
    }
  }

  /* set non-zero values
   *    __                          __
   *  0|                              |
   *  1|    x                 x       |
   *  2| x                 x          |
   *  3|    x     x                   |
   *  4|                x             |
   *  5|                         x    |
   *  6|                              |
   *  7|                              |
   *  8|__                          __|
   *     0  1  2  3  4  5  6  7  8  9
   */

  l_matDense[1][1] = 1;
  l_matDense[1][7] = 1;
  l_matDense[2][0] = 1;
  l_matDense[2][6] = 1;
  l_matDense[3][1] = 1;
  l_matDense[3][3] = 1;
  l_matDense[4][5] = 1;
  l_matDense[5][8] = 1;

  // convert to coordinate format
  t_matCrd l_matCrd;
  edge::linalg::Matrix::denseToCrd( 9, 9,
                                    l_matDense[0],
                                    l_matCrd );

  unsigned int l_nzBlockIn[2][2];
  unsigned int l_nzBlockOut[2][2];
  for( unsigned short l_ro = 0; l_ro < 2; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 2; l_co++ ) {
      l_nzBlockIn[l_ro][l_co]  = std::numeric_limits<unsigned int>::max();
      l_nzBlockOut[l_ro][l_co] = std::numeric_limits<unsigned int>::max()-1;
    }
  }

  // determine blocks size of entire matrix (no bounds)
  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == 1 );
  REQUIRE( l_nzBlockOut[0][1] == 5 );
  REQUIRE( l_nzBlockOut[1][0] == 0 );
  REQUIRE( l_nzBlockOut[1][1] == 8 );

  // try bounds given by matrix-size
  l_nzBlockIn[0][0] = l_nzBlockIn[1][0] = 0;
  l_nzBlockIn[0][1] = l_nzBlockIn[1][1] = 8;

  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == 1 );
  REQUIRE( l_nzBlockOut[0][1] == 5 );
  REQUIRE( l_nzBlockOut[1][0] == 0 );
  REQUIRE( l_nzBlockOut[1][1] == 8 );

  // true bounds given by non-zero subblock
  l_nzBlockIn[0][0] = 1;
  l_nzBlockIn[0][1] = 5;
  l_nzBlockIn[1][0] = 0;
  l_nzBlockIn[1][1] = 8;

  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == 1 );
  REQUIRE( l_nzBlockOut[0][1] == 5 );
  REQUIRE( l_nzBlockOut[1][0] == 0 );
  REQUIRE( l_nzBlockOut[1][1] == 8 );

  // try submatrix
  l_nzBlockIn[0][0] = 2;
  l_nzBlockIn[0][1] = 4;
  l_nzBlockIn[1][0] = 1;
  l_nzBlockIn[1][1] = 9;

  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == 2 );
  REQUIRE( l_nzBlockOut[0][1] == 4 );
  REQUIRE( l_nzBlockOut[1][0] == 1 );
  REQUIRE( l_nzBlockOut[1][1] == 6 );

  // try submatrix with single nz-entry
  l_nzBlockIn[0][0] = 5;
  l_nzBlockIn[0][1] = std::numeric_limits<unsigned int>::max();
  l_nzBlockIn[1][0] = std::numeric_limits<unsigned int>::max();
  l_nzBlockIn[1][1] = std::numeric_limits<unsigned int>::max();

  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == 5 );
  REQUIRE( l_nzBlockOut[0][1] == 5 );
  REQUIRE( l_nzBlockOut[1][0] == 8 );
  REQUIRE( l_nzBlockOut[1][1] == 8 );

  // try submatrix without any nz-entries
  l_nzBlockIn[0][0] = 5;
  l_nzBlockIn[0][1] = std::numeric_limits<unsigned int>::max();
  l_nzBlockIn[1][0] = std::numeric_limits<unsigned int>::max();
  l_nzBlockIn[1][1] = 7;

  edge::linalg::Matrix::getBlockNz( l_matCrd,
                                    l_nzBlockIn,
                                    l_nzBlockOut );

  REQUIRE( l_nzBlockOut[0][0] == std::numeric_limits<unsigned int>::max() );
  REQUIRE( l_nzBlockOut[0][1] == std::numeric_limits<unsigned int>::max() );
  REQUIRE( l_nzBlockOut[1][0] == std::numeric_limits<unsigned int>::max() );
  REQUIRE( l_nzBlockOut[1][1] == std::numeric_limits<unsigned int>::max() );
}

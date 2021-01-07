/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Tests the generation of fused matrix kernels.
 **/
#include <catch.hpp>
#define private public
#include "MmXsmmFused.hpp"
#undef private


TEST_CASE( "Tests the fused matrix multiplications in libxsmm (dense tensor)x(sparse matrix)", "[data][MmXsmmFusedDenseSparse]" ) {
  // set up matrix structures
#include "MmSparse.test.inc"

  float l_input1[9][10][1] = {0};
  float l_result1[9][10][1] = {0};
  float l_input2[9][10][2] = {0};
  float l_result2[9][10][2] = {0};
  float l_input4[9][10][4] = {0};
  float l_result4[9][10][4] = {0};
  float l_input8[9][10][8] = {0};
  float l_result8[9][10][8] = {0};
  float l_input16[9][10][16] = {0};
  float l_result16[9][10][16] = {0};
  float l_input32[9][10][32] = {0};
  float l_result32[9][10][32] = {0};

  // init input
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 32; l_cr++ ) {
        if( l_cr < 1 ) {
          l_input1[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 2 ) {
          l_input2[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 4 ) {
          l_input4[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 8 ) {
          l_input8[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
        }


        if( l_cr < 16 ) {
          l_input16[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
        }

        l_input32[l_ro][l_co][l_cr] = l_denseSparseInput[l_ro][l_co] * (l_cr+1);
      }
    }
  }

  // init libxsmm-kernels
  edge::data::MmXsmmFused< float > m_mm;

  for( unsigned short l_nc = 1; l_nc <= 32; l_nc=l_nc*2 ) {
    m_mm.add( 0, // group
              l_nc, // nCrs
              false, // csc
              l_denseSparseCscColPtrs, // column pointer
              l_denseSparseCscRowIds, // row indices
              l_denseSparseCscNzs, // values
              9, // m
              10, // n
              10, // k
              10, // ldA
              0, // ldB
              10, // ldC
              1.0f, // alpha
              1.0f, // beta
              LIBXSMM_GEMM_PREFETCH_NONE ); 
  }

  // run 1-kernel
  m_mm.m_kernels[0][0]( l_input1[0][0],
                        l_denseSparseCscNzs,
                        l_result1[0][0] );

  // check 1-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      REQUIRE( l_result1[l_ro][l_co][0] == Approx( l_denseSparseRefResult[l_ro][l_co] * 1 ) );
    }
  }

  // run 2-kernel
  m_mm.m_kernels[0][1]( l_input2[0][0],
                        l_denseSparseCscNzs,
                        l_result2[0][0] );

  // check 2-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        REQUIRE( l_result2[l_ro][l_co][l_cr] == Approx( l_denseSparseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 4-kernel
  m_mm.m_kernels[0][2]( l_input4[0][0],
                        l_denseSparseCscNzs,
                        l_result4[0][0] );

  // check 4-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 4; l_cr++ ) {
        REQUIRE( l_result4[l_ro][l_co][l_cr] == Approx( l_denseSparseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 8-kernel
  m_mm.m_kernels[0][3]( l_input8[0][0],
                        l_denseSparseCscNzs,
                        l_result8[0][0] );

  // check 8-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 8; l_cr++ ) {
        REQUIRE( l_result8[l_ro][l_co][l_cr] == Approx( l_denseSparseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 16-kernel
  m_mm.m_kernels[0][4]( l_input16[0][0],
                        l_denseSparseCscNzs,
                        l_result16[0][0] );

  // check 16-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_result16[l_ro][l_co][l_cr] == Approx( l_denseSparseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 32-kernel
  m_mm.m_kernels[0][5]( l_input32[0][0],
                        l_denseSparseCscNzs,
                        l_result32[0][0] );

  // check 32-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 10; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 32; l_cr++ ) {
        REQUIRE( l_result32[l_ro][l_co][l_cr] == Approx( l_denseSparseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }
}

TEST_CASE( "Tests the fused matrix multiplications in libxsmm (sparse matrix)x(dense tensor)", "[data][MmXsmmFusedSparseDense]" ) {
  // set up matrix structures
#include "MmSparse.test.inc"

  float l_input1[9][20][1] = {0};
  float l_result1[9][20][1] = {0};
  float l_input2[9][20][2] = {0};
  float l_result2[9][20][2] = {0};
  float l_input4[9][20][4] = {0};
  float l_result4[9][20][4] = {0};
  float l_input8[9][20][8] = {0};
  float l_result8[9][20][8] = {0};
  float l_input16[9][20][16] = {0};
  float l_result16[9][20][16] = {0};
  float l_input32[9][20][32] = {0};
  float l_result32[9][20][32] = {0};

  // init input
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 32; l_cr++ ) {
        if( l_cr < 1 ) {
          l_input1[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 2 ) {
          l_input2[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 4 ) {
          l_input4[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
        }

        if( l_cr < 8 ) {
          l_input8[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
        }


        if( l_cr < 16 ) {
          l_input16[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
        }

        l_input32[l_ro][l_co][l_cr] = l_sparseDenseInput[l_ro][l_co] * (l_cr+1);
      }
    }
  }

  // init libxsmm-kernels
  edge::data::MmXsmmFused< float > m_mm;

  for( unsigned short l_nc = 1; l_nc <= 32; l_nc=l_nc*2 ) {
    m_mm.add( 0, // group
              l_nc, // nCrs
              true, // csr
              l_sparseDenseCsrRowPtrs, // row pointers
              l_sparseDenseCsrColIds, // row indices
              l_sparseDenseCsrNzs, // values
              9, // m
              20, // n
              9, // k
              0, // ldA
              20, // ldB
              20, // ldC
              1.0f, // alpha
              0.0f, // beta
              LIBXSMM_GEMM_PREFETCH_NONE );
  }

  // run 1-kernel
  m_mm.m_kernels[0][0]( l_sparseDenseCsrNzs,
                        l_input1[0][0],
                        l_result1[0][0] );

  // check 1-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      REQUIRE( l_result1[l_ro][l_co][0] == Approx( l_sparseDenseRefResult[l_ro][l_co] * 1 ) );
    }
  }

  // run 2-kernel
  m_mm.m_kernels[0][1]( l_sparseDenseCsrNzs,
                        l_input2[0][0],
                        l_result2[0][0] );


  // check 2-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ ) {
        REQUIRE( l_result2[l_ro][l_co][l_cr] == Approx( l_sparseDenseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 4-kernel
  m_mm.m_kernels[0][2]( l_sparseDenseCsrNzs,
                        l_input4[0][0],
                        l_result4[0][0] );


  // check 4-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 4; l_cr++ ) {
        REQUIRE( l_result4[l_ro][l_co][l_cr] == Approx( l_sparseDenseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 8-kernel
  m_mm.m_kernels[0][3]( l_sparseDenseCsrNzs,
                        l_input8[0][0],
                        l_result8[0][0] );


  // check 8-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 8; l_cr++ ) {
        REQUIRE( l_result8[l_ro][l_co][l_cr] == Approx( l_sparseDenseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 16-kernel
  m_mm.m_kernels[0][4]( l_sparseDenseCsrNzs,
                        l_input16[0][0],
                        l_result16[0][0] );


  // check 16-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_result16[l_ro][l_co][l_cr] == Approx( l_sparseDenseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }

  // run 32-kernel
  m_mm.m_kernels[0][5]( l_sparseDenseCsrNzs,
                        l_input32[0][0],
                        l_result32[0][0] );


  // check 32-results
  for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
    for( unsigned short l_co = 0; l_co < 20; l_co++ ) {
      for( unsigned short l_cr = 0; l_cr < 32; l_cr++ ) {
        REQUIRE( l_result32[l_ro][l_co][l_cr] == Approx( l_sparseDenseRefResult[l_ro][l_co] * (l_cr+1) ) );
      }
    }
  }
}
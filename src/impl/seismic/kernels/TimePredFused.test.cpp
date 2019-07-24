/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Tests the optimized time prediction for fused forward simulation.
 **/
#include <catch.hpp>
#define private public
#include "TimePredFused.hpp"
#undef private


TEST_CASE( "Optimized elastic ADER time prediction for fused forward simulations.", "[elastic][TimePredFused]" ) {
  // set up matrix structures
#include "TimePred.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredFused< float,
                                         0,
                                         TET4,
                                         4,
                                         4,
                                         16 > l_pred( nullptr,
                                                      l_dynMem );

  float l_scratch[9][20][16];
  float l_ders[4][9][20][16];
  float l_dofsFused[9][20][16];
  float l_tDofs[9][20][16];

  // duplicate DOFs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_dofsFused[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];
      }
    }
  }

  // compute time prediction
  l_pred.ck( 0.017,
             l_starSpE,
             nullptr,
             nullptr,
             l_dofsFused,
             nullptr,
             l_scratch,
             l_ders,
             nullptr,
             l_tDofs,
             nullptr );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_tDofs[l_qt][l_md][l_cr] == Approx( l_refEtDofs[l_qt][l_md] ) );
      }
    }
  }
}

TEST_CASE( "Optimized viscoelastic ADER time prediction for fused forward simulations.", "[visco][TimePredFused]" ) {
  // set up matrix structures
#include "TimePred.test.inc"

  float l_scratch[9][20][16];
  float l_dersE[4][9][20][16];
  float l_dersA[2][4][6][20][16];
  float l_dofsFusedE[9][20][16];
  float l_dofsFusedA[2][6][20][16];
  float l_tDofsE[9][20][16];
  float l_tDofsA[2][6][20][16];

  // duplicate DOFs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_dofsFusedE[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];

        if( l_qt < 6) {
          for( unsigned short l_rm = 0; l_rm < 2; l_rm++ ) {
            l_dofsFusedA[l_rm][l_qt][l_md][l_cr] = l_dofsA[l_rm][l_qt][l_md];
          }
        }
      }
    }
  }

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredFused< float,
                                         2,
                                         TET4,
                                         4,
                                         4,
                                         16 > l_pred( l_rfs,
                                                      l_dynMem );

  // compute time prediction
  l_pred.ck( 0.017,
             l_starSpE,
             l_starSpA,
             l_srcSpA,
             l_dofsFusedE,
             l_dofsFusedA,
             l_scratch,
             l_dersE,
             l_dersA,
             l_tDofsE,
             l_tDofsA );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_tDofsE[l_qt][l_md][l_cr] == Approx( l_refVtDofsE[l_qt][l_md] ) );
      }
    }
  }

  for( unsigned short l_rm = 0; l_rm < 2; l_rm++ ) {
    for( unsigned short l_qt = 0; l_qt < 6; l_qt++ ) {
      for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
        for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
          REQUIRE( l_tDofsA[l_rm][l_qt][l_md][l_cr] == Approx( l_refVtDofsA[l_rm][l_qt][l_md] ) );
        }
      }
    }
  }
}
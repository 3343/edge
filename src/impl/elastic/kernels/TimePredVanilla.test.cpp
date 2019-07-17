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
 * Tests the vanilla time prediction.
 **/
#include <catch.hpp>
#define private public
#include "TimePredVanilla.hpp"
#undef private


TEST_CASE( "Elastic ADER time prediction using vanilla kernels.", "[elastic][TimePredVanilla]" ) {
  // set up matrix structures
#include "TimePred.test.inc"


  float l_scratch[9][20][1];
  float l_dersE[4][9][20][1];
  float l_tDofsE[9][20][1];

  // set up kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredVanilla< float,
                                           0,
                                           TET4,
                                           4,
                                           4,
                                           1 > l_predElastic( nullptr,
                                                              l_dynMem );

  // compute time prediction
  l_predElastic.ck( 0.017,
                    l_starE,
                    nullptr,
                    nullptr,
                    (float (*)[20][1]) l_dofsE,
                    nullptr,
                    l_scratch,
                    l_dersE,
                    nullptr,
                    l_tDofsE,
                    nullptr );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      REQUIRE( l_tDofsE[l_qt][l_md][0] == Approx( l_refEtDofs[l_qt][l_md] ) );
    }
  }
}

TEST_CASE( "Viscoelastic ADER time prediction using vanilla kernels.", "[visco][TimePredVanilla]" ) {
  // setup matrix structures
#include "TimePred.test.inc"

  float l_scratch[9][20][1];
  float l_dersE[4][9][20][1];
  float l_tDofsE[9][20][1];
  float l_tDofsA[2][6][20][1];
  float l_dersA[2][4][6][20][1];

  // set up kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredVanilla< float,
                                           2,
                                           TET4,
                                           4,
                                           4,
                                           1 > l_pred( l_rfs,
                                                       l_dynMem );

  // compute time prediction
  l_pred.ck( 0.017,
             l_starE,
             (float (*)[18]) l_starA,
             (float (*)[36]) l_srcA,
             (float (*)[20][1]) l_dofsE,
             (float (*)[6][20][1]) l_dofsA,
             l_scratch,
             l_dersE,
             l_dersA,
             l_tDofsE,
             l_tDofsA );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      REQUIRE( l_tDofsE[l_qt][l_md][0] == Approx( l_refVtDofsE[l_qt][l_md] ) );
    }
  }

  for( unsigned short l_rm = 0; l_rm < 2; l_rm++ ) {
    for( unsigned short l_qt = 0; l_qt < 6; l_qt++ ) {
      for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
        REQUIRE( l_tDofsA[l_rm][l_qt][l_md][0] == Approx( l_refVtDofsA[l_rm][l_qt][l_md] ) );
      }
    }
  }
}
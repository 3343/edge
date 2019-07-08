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


TEST_CASE( "ADER time prediction using vanilla kernels.", "[TimePredVanilla][seismic]" ) {
  edge::data::Dynamic l_dynMem;

  edge::elastic::kernels::TimePredVanilla< float,
                                           TET4,
                                           4,
                                           4,
                                           1 > l_pred( l_dynMem );

  // setup matrix structures
  #include "TimePred.test.inc"

  float l_scratch[9][20][1];
  float l_ders[4][9][20][1];
  float l_tDofs[9][20][1];

  // compute time prediction
  l_pred.ck( 0.017, l_starMats, (float (*)[20][1])  l_dofs, l_scratch, l_ders, l_tDofs );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      REQUIRE( l_tDofs[l_qt][l_md][0] == Approx( l_tDofsRef[l_qt][l_md] ) );
    }
  }
}
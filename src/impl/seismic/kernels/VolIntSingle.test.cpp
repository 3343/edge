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
 * Tests the optimized volume integration for single forward simulations.
 **/
#include <catch.hpp>
#include "VolIntSingle.hpp"

TEST_CASE( "Optimized elastic volume integration for single forward simulations.", "[elastic][VolIntSingle]" ) {
  // set up matrix structures
#include "VolInt.test.inc"

  // volume kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::VolIntSingle< float,
                                        0,
                                        TET4,
                                        4 > l_vol( nullptr, l_dynMem );

  float l_scratch[9][20][1];

  // compute volume integration
  l_vol.apply( (float (*)[81])    l_starE,
                                  nullptr,
                                  nullptr,
               (float (*)[20][1]) l_tDofsE,
                                  nullptr,
               (float (*)[20][1]) l_dofsE,
                                  nullptr,
                                  l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      REQUIRE( l_dofsE[l_qt][l_md] == Approx( l_refEdofs[l_qt][l_md] ) );
    }
  }
}

TEST_CASE( "Optimized viscoelastic volume integration for single forward simulations.", "[visco][VolIntSingle]" ) {
  // set up matrix structures
#include "VolInt.test.inc"

  float l_scratch[9][20][1];

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::VolIntSingle< float,
                                        2,
                                        TET4,
                                        4 > l_vol( l_rfs,
                                                   l_dynMem );

  // compute solution
  l_vol.apply( (float (*)[81])       l_starE,
               (float (*)[6*3])      l_starA,
               (float (*)[6*6])      l_srcA,
               (float (*)[20][1])    l_tDofsE,
               (float (*)[6][20][1]) l_tDofsA,
               (float (*)[20][1])    l_dofsE,
               (float (*)[6][20][1]) l_dofsA,
                                     l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
      REQUIRE( l_dofsE[l_qt][l_md] == Approx( l_refVdofsE[l_qt][l_md] ) );
    }
  }

  for( unsigned short l_rm = 0; l_rm < 2; l_rm++ ) {
    for( unsigned short l_qt = 0; l_qt < 6; l_qt++ ) {
      for( unsigned short l_md = 0; l_md < 20; l_md++ ) {
        REQUIRE( l_dofsA[l_rm][l_qt][l_md] == Approx( l_refVdofsA[l_rm][l_qt][l_md] ) );
      }
    }
  }
}
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
 * Tests the fused surface integration.
 **/
#include <catch.hpp>
#define private public
#include "SurfIntFused.hpp"
#undef private


TEST_CASE( "Local surface integration for fused simulations.", "[SurfIntLocalFused][seismic]" ) {
  edge::data::Dynamic l_dynMem;

  edge::elastic::kernels::SurfIntFused< float,
                                        TET4,
                                        3,
                                        16 > l_surf( l_dynMem );

  // setup matrix structures
#include "SurfInt.test.inc"

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofs[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofs[l_qt][l_md];
      }
    }
  }

  // compute local surface integration
  l_surf.local( l_fSolv,
                l_tDofsFused,
                l_dofsFused,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_dofsLocalRef[l_qt][l_md] ) );
      }
    }
  }
}


TEST_CASE( "Neighboring surface integration for fused simulations.", "[SurfIntNeighFused][seismic]" ) {
  edge::data::Dynamic l_dynMem;

  edge::elastic::kernels::SurfIntFused< float,
                                        TET4,
                                        3,
                                        16 > l_surf( l_dynMem );

  // setup matrix structures
#include "SurfInt.test.inc"

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofs[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofs[l_qt][l_md];
      }
    }
  }


  // compute neighboring surface integration
  l_surf.neigh( 3,
                1,
                2,
                l_fSolv[0],
                l_tDofsFused,
                l_dofsFused,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_dofsNeighRef[l_qt][l_md] ) );
      }
    }
  }
}


TEST_CASE( "Neighboring surface integration in the presence of a free surface for used simulations.", "[SurfIntNeighFsFused][seismic]" ) {
  edge::data::Dynamic l_dynMem;

  edge::elastic::kernels::SurfIntFused< float,
                                        TET4,
                                        3,
                                        16 > l_surf( l_dynMem );

  // setup matrix structures
#include "SurfInt.test.inc"

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofs[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofs[l_qt][l_md];
      }
    }
  }

  // compute neighboring surface integration with at a free-surface
  l_surf.neigh( 2,
                std::numeric_limits< unsigned short >::max(),
                std::numeric_limits< unsigned short >::max(),
                l_fSolv[0],
                l_tDofsFused,
                l_dofsFused,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_dofsNeighFfRef[l_qt][l_md] ) );
      }
    }
  }
}
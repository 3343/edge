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


TEST_CASE( "Local elastic surface integration for fused simulations.", "[elastic][SurfIntLocalFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        0,
                                        TET4,
                                        3,
                                        16 > l_surf( nullptr,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];
      }
    }
  }

  // compute local surface integration
  l_surf.local( (float (*)[9*9]) l_fSolvE,
                                 nullptr,
                                 l_tDofsFused,
                                 l_dofsFused,
                                 nullptr,
                                 l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_refEdofs[l_qt][l_md] ) );
      }
    }
  }
}

TEST_CASE( "Neighboring elastic surface integration for fused simulations.", "[elastic][SurfIntNeighFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        0,
                                        TET4,
                                        3,
                                        16 > l_surf( nullptr,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];
      }
    }
  }


  // compute neighboring surface integration
  l_surf.neigh( 3,
                1,
                2,
                l_fSolvE[0][0],
                nullptr,
                l_tDofsFused,
                l_dofsFused,
                nullptr,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_refEneighDofs[l_qt][l_md] ) );
      }
    }
  }
}


TEST_CASE( "Neighboring elastic surface integration in the presence of a free surface for used simulations.", "[elastic][SurfIntNeighFsFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        0,
                                        TET4,
                                        3,
                                        16 > l_surf( nullptr,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFused[9][10][16];
  float l_tDofsFused[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFused[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFused[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];
      }
    }
  }

  // compute neighboring surface integration with at a free-surface
  l_surf.neigh( 2,
                std::numeric_limits< unsigned short >::max(),
                std::numeric_limits< unsigned short >::max(),
                l_fSolvE[0][0],
                nullptr,
                l_tDofsFused,
                l_dofsFused,
                nullptr,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFused[l_qt][l_md][l_cr] == Approx( l_refEneighFdofs[l_qt][l_md] ) );
      }
    }
  }
}

TEST_CASE( "Local viscoelastic surface integration for fused simulations.", "[visco][SurfIntLocalFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        3,
                                        TET4,
                                        3,
                                        16 > l_surf( l_rfs,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFusedE[9][10][16];
  float l_dofsFusedA[3][6][10][16];
  float l_tDofsFusedE[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFusedE[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFusedE[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];

        if( l_qt < 6 ) {
          for( unsigned short l_rm = 0; l_rm < 3; l_rm++ ) {
            l_dofsFusedA[l_rm][l_qt][l_md][l_cr] = l_dofsA[l_rm][l_qt][l_md];
          }
        }
      }
    }
  }

  // compute local surface integration
  l_surf.local( (float (*) [9*9]) l_fSolvE,
                (float (*) [6*9]) l_fSolvA,
                                  l_tDofsFusedE,
                                  l_dofsFusedE,
                                  l_dofsFusedA,
                                  l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFusedE[l_qt][l_md][l_cr] == Approx( l_refEdofs[l_qt][l_md] ) );

        if( l_qt < 6) {
          for( unsigned short l_rm = 0; l_rm < 3; l_rm++ ) {
            REQUIRE( l_dofsFusedA[l_rm][l_qt][l_md][l_cr] == Approx( l_refVdofsA[l_rm][l_qt][l_md] ) );
          }
        }
      }
    }
  }
}

TEST_CASE( "Neighboring viscoelastic surface integration for fused simulations.", "[visco][SurfIntNeighFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        3,
                                        TET4,
                                        3,
                                        16 > l_surf( l_rfs,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFusedE[9][10][16];
  float l_upsFusedA[6][10][16];
  float l_tDofsFusedE[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFusedE[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFusedE[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];

        if( l_qt < 6 ) {
          l_upsFusedA[l_qt][l_md][l_cr] = 0;
        }
      }
    }
  }

  // compute neighboring surface integration
  l_surf.neigh( 3,
                1,
                2,
                l_fSolvE[0][0],
                l_fSolvA[0][0],
                l_tDofsFusedE,
                l_dofsFusedE,
                l_upsFusedA,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFusedE[l_qt][l_md][l_cr] == Approx( l_refEneighDofs[l_qt][l_md] ) );
      }
    }
  }

  // check the results
  for( unsigned short l_qt = 0; l_qt < 6; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_upsFusedA[l_qt][l_md][l_cr] == Approx( l_refVneighDofsA[l_qt][l_md] ) );
      }
    }
  }
}

TEST_CASE( "Neighboring viscoelastic surface integration in the presence of a free surface for fused simulations.", "[visco][SurfIntNeighFsFused]" ) {
  // set up matrix structures
#include "SurfInt.test.inc"

  // kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::SurfIntFused< float,
                                        3,
                                        TET4,
                                        3,
                                        16 > l_surf( l_rfs,
                                                     l_dynMem );

  float l_scratch[2][9][6][16];
  float l_dofsFusedE[9][10][16];
  float l_upsFusedA[6][10][16];
  float l_tDofsFusedE[9][10][16];

  // duplicate DOFs and tDofs for fused config
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        l_tDofsFusedE[l_qt][l_md][l_cr] = l_tDofsE[l_qt][l_md];
        l_dofsFusedE[l_qt][l_md][l_cr] = l_dofsE[l_qt][l_md];

        if( l_qt < 6 ) {
          l_upsFusedA[l_qt][l_md][l_cr] = 0;
        }
      }
    }
  }

  // compute neighboring surface integration
  l_surf.neigh( 2,
                std::numeric_limits< unsigned short >::max(),
                std::numeric_limits< unsigned short >::max(),
                l_fSolvE[0][0],
                l_fSolvA[0][0],
                l_tDofsFusedE,
                l_dofsFusedE,
                l_upsFusedA,
                l_scratch );

  // check the results
  for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_dofsFusedE[l_qt][l_md][l_cr] == Approx( l_refEneighFdofs[l_qt][l_md] ) );
      }
    }
  }

  // check the results
  for( unsigned short l_qt = 0; l_qt < 6; l_qt++ ) {
    for( unsigned short l_md = 0; l_md < 10; l_md++ ) {
      for( unsigned short l_cr = 0; l_cr < 16; l_cr++ ) {
        REQUIRE( l_upsFusedA[l_qt][l_md][l_cr] == Approx( l_refVneighFdofsA[l_qt][l_md] ) );
      }
    }
  }
}
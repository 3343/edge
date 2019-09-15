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
 * Time prediction unit tests.
 **/
#include <catch.hpp>
#define private public
#include "TimePred.hpp"
#undef private


TEST_CASE( "Time prediction unit tests for the advection equation.", "[advection][TimePred]" ) {
  // set up unit test input
#include "TimePred.test.inc"

  // set up kernel
  edge::data::Dynamic l_dynMem;
  edge::advection::kernels::TimePred< float,
                                      TET4,
                                      3,
                                      3,
                                      2 > l_pred( l_dynMem );

  // set up data structures for time prediction
  float l_ders[3][10][2];
  float l_tInt[10][2];
  float l_dofs[10][2];
  for( unsigned short l_md = 0; l_md < 10; l_md++ )
    for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
      l_dofs[l_md][l_cr] = l_utDofs[l_md];

  // call kernel
  l_pred.ck( l_utDt,
             l_utStar,
             l_dofs,
             l_ders,
             l_tInt );

  // check results
  for( unsigned short l_de = 0; l_de < 3; l_de++ )
    for( unsigned short l_md = 0; l_md < 10; l_md++ )
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
        REQUIRE( l_ders[l_de][l_md][l_cr] == Approx( l_utDers[l_de][l_md] ) );

  for( unsigned short l_md = 0; l_md < 10; l_md++ )
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
        REQUIRE( l_tInt[l_md][l_cr] == Approx( l_utTint[l_md] ) );

  // perform rate-2 sub-sampling in time
  float l_subInt2[10][2];

  l_pred.integrate( l_utDt*0.5,
                    l_ders,
                    l_subInt2 );

  // check the results
  for( unsigned short l_md = 0; l_md < 10; l_md++ )
      for( unsigned short l_cr = 0; l_cr < 2; l_cr++ )
        REQUIRE( l_subInt2[l_md][l_cr] == Approx( l_utTint12[l_md] ) );
}
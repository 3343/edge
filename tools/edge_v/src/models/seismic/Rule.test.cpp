/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Tests the velocity rules
 **/

#include <catch.hpp>
#define private public
#include "Rule.h"
#undef private

TEST_CASE( "Tests the application of velocity rules for the 2018 HighF setting.", "[rule][highf2018]" ) {
  std::string l_expr = R"V0G0N(
  // Decription, following https://scec.usc.edu/scecpedia/HighF_2018:
  //
  //  1. Set Min Vs=500 m/s
  //  2. If Vs was lower than 500 m/s and adjusted, then adjust Vp with original Vp/Vs ratio (so that we don’t have the automatic Vs/Vp of 3). We may want to set a minimum value of Vp (Rob to check)
  //  3. Then set Max Vp/Vs= 3, if lower Vp to maintain the max of 3 ratio

  // limit vs to 500 m/s
  if( vs < 500 ) {
    //  derive vp/vs ratio
    var sca := vp / vs;
    // adjust velocities
    vs := 500;
    vp := 500*sca;
  };

  // limit vp/vs ratio to 3
  var sca := vp / vs;
  if( sca > 3 ) {
    vp := 3*vs;
  };
  )V0G0N";

  edge_v::models::seismic::Rule l_rule( l_expr );

  float l_vp = 700;
  float l_vs = 600;
  float l_rho = 2000;
  l_rule.apply( l_vp, l_vs, l_rho );


  REQUIRE( l_vs == Approx(600) );
  REQUIRE( l_vp == Approx(700) );
  REQUIRE( l_rho == Approx(2000) );

  l_vp = 700;
  l_vs = 200;
  l_rho = 2000;
  l_rule.apply( l_vp, l_vs, l_rho );

  REQUIRE( l_vs == Approx(500) );
  REQUIRE( l_vp == Approx( (l_vp/l_vs)*500 ) );
  REQUIRE( l_rho == Approx(2000) );

  l_vp = 900;
  l_vs = 200;
  l_rho = 2000;
  l_rule.apply( l_vp, l_vs, l_rho );

  REQUIRE( l_vs == Approx(500) );
  REQUIRE( l_vp == Approx(1500) );
  REQUIRE( l_rho == Approx(2000) );
}
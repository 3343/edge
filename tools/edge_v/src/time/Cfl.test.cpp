/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Tests the time step derivation, based on stability constraints.
 **/
#include <catch.hpp>
#include "models/Constant.h"
#define private public
#include "Cfl.h"
#undef private

TEST_CASE( "Tests the step derivation.", "[time][cfl]" ) {
  // set up dummy data
  edge_v::t_idx l_nVes = 4;
  edge_v::t_idx l_nEls = 4;
  edge_v::t_idx l_elVe[4 * 3] = {};
  double l_veCrds[4][3] = {};
  double l_inDia[4] = { 2.0, 3.0, 1.5, 2.5 };

  edge_v::models::Constant l_mod( 2.0 );

  // set up time steps
  edge_v::time::Cfl l_cfl( edge_v::TRIA3,
                           l_nVes,
                           l_nEls,
                           l_elVe,
                           l_veCrds,
                           l_inDia,
                           l_mod );

  // check the results
  REQUIRE( l_cfl.m_tsAbsMin == Approx(0.75) );
  REQUIRE( l_cfl.getTimeSteps()[0] == Approx( 1.0  / 0.75 ) );
  REQUIRE( l_cfl.getTimeSteps()[1] == Approx( 1.5  / 0.75 ) );
  REQUIRE( l_cfl.getTimeSteps()[2] == Approx( 0.75 / 0.75 ) );
  REQUIRE( l_cfl.getTimeSteps()[3] == Approx( 1.25 / 0.75 ) );
}
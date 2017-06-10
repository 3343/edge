/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Unit test of half space class.
 **/

#include <catch.hpp>
#include "HalfSpace.hpp"

#include <iostream>

TEST_CASE( "2D Half space: Inside/Outside.", "[inside2D][HalfSpace]" ) {
  double l_origin[2];
  double l_normal[2];
  double l_pt[2];

  // construct half-space
  l_origin[0] =  13;
  l_origin[1] = -12;
  l_normal[0] =  20;
  l_normal[1] =   0;

  edge::linalg::HalfSpace<double, 2> l_hs1( l_origin, l_normal );

  // check properties of half-space
  l_pt[0] =  13;
  l_pt[1] = -12;
  REQUIRE( l_hs1.inside( l_pt )        == true  );
  REQUIRE( l_hs1.inside( l_pt, -1E-8 ) == true  );
  REQUIRE( l_hs1.inside( l_pt,  1E-6 ) == false );

  l_pt[0] =  13;
  l_pt[1] = -99999;
  REQUIRE( l_hs1.inside( l_pt ) == true  );

  l_pt[0] =  15;
  l_pt[1] = -591;
  REQUIRE( l_hs1.inside( l_pt ) == true  );

  l_pt[0] =  11;
  l_pt[1] = 591;
  REQUIRE( l_hs1.inside( l_pt ) == false  );

  // construct new half-space
  l_origin[0] = -3;
  l_origin[1] =  2;
  l_normal[0] =  1;
  l_normal[1] = -1;
  edge::linalg::HalfSpace<double, 2> l_hs2( l_origin, l_normal );

  // check properties of half-space
  l_pt[0] = 0;
  l_pt[1] = 0;
  REQUIRE( l_hs2.inside( l_pt ) == true  );

  l_pt[0] = -4;
  l_pt[1] =  3;
  REQUIRE( l_hs2.inside( l_pt ) == false  );
}

TEST_CASE( "3D Half space: Inside/Outside.", "[inside3D][HalfSpace]" ) {
  double l_origin[3];
  double l_normal[3];
  double l_pt[3];

  // construct half-space
  l_origin[0] =  0;
  l_origin[1] =  0;
  l_origin[2] =  0;
  l_normal[0] =  1;
  l_normal[1] = -1;
  l_normal[2] =  3;

  edge::linalg::HalfSpace<double, 2> l_hs1( l_origin, l_normal );

  // check properties of half-space
  l_pt[0] = 0;
  l_pt[1] = 0;
  l_pt[2] = 0;
  REQUIRE( l_hs1.inside( l_pt ) == true );

  l_pt[0] = 2;
  l_pt[1] = -4;
  l_pt[2] = 5;
  REQUIRE( l_hs1.inside( l_pt ) == true );

  l_pt[0] = -2;
  l_pt[1] = 4;
  l_pt[2] = -5;
  REQUIRE( l_hs1.inside( l_pt ) == false );
}

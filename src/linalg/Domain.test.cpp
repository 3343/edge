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
 * Unit test of domain class.
 **/

#include <catch.hpp>
#include "Domain.hpp"
#include "HalfSpace.hpp"

TEST_CASE( "2D Domain: Inside/Outside.", "[inside2D][Domain]" ) {
  double l_pt[2];
  double l_origin[2];
  double l_normal[2];

  // create domain 1
  edge::linalg::Domain< double, 2, edge::linalg::HalfSpace > l_dom1;

  // no objects -> point is outside
  l_pt[0] = 1;
  l_pt[1] = 2;
  REQUIRE( l_dom1.inside ( l_pt ) == true );

  // create domain 1
  edge::linalg::Domain< double, 2, edge::linalg::HalfSpace > l_dom2;

  // add halfspace
  l_origin[0] = -41;
  l_origin[1] = -24;
  l_normal[0] =  90;
  l_normal[1] =   0;

  edge::linalg::HalfSpace<double, 2> l_hs1( l_origin, l_normal );
  l_dom2.add( l_hs1 );

  REQUIRE( l_dom2.inside( l_pt ) == true );

  l_origin[0] = -10;
  l_origin[1] = -994;
  l_normal[0] = -10;
  l_normal[1] =  0;
  l_hs1.setOrigin( l_origin );
  l_hs1.setNormal( l_normal );

  l_dom2.add( l_hs1 );
  REQUIRE( l_dom2.inside( l_pt ) == false );
}

/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Unit tests for receiver output.
 **/
#include <catch.hpp>
#include "Receivers.h"

TEST_CASE( "Receivers: Initialization", "[receivers][init]" ) {
  edge::io::Receivers l_recv;

  t_enLayout l_elLayout;

  l_elLayout.timeGroups.resize( 2 );
  l_elLayout.timeGroups[0].nEntsOwn    = 2;
  l_elLayout.timeGroups[0].nEntsNotOwn = 1;
  l_elLayout.timeGroups[1].nEntsOwn    = 3;
  l_elLayout.timeGroups[1].nEntsNotOwn = 4;
#ifdef PP_T_ELEMENTS_TET4
  // setup receiver coordinates
  real_mesh l_recvCrds[6][3] = { { -1.00, -1.00, -2.00 }, // not part of the "domain"
                                 {  0.15,  0.15,  0.15 }, // inside a tet
                                 {  0.00,  0.05,  0.00 }, // hits an edge
                                 {  0.00,  0.00,  0.00 }, // vertex
                                 {  2.00,  0.00,  0.00 }, // outside
                                 {  0.20,  0.10,  0.05 }  // inside
                               };

  t_vertexChars l_veChars[8] = { {{0.0, 0.0, 0.0}, 0}, // reference tet as dummy tet
                                 {{1.0, 0.0, 0.0}, 0},
                                 {{0.0, 1.0, 0.0}, 0},
                                 {{0.0, 0.0, 1.0}, 0}, 
                                 {{0.0, 0.0, 8.0}, 0}, // shifted ref-tet
                                 {{1.0, 0.0, 8.0}, 0},
                                 {{0.0, 1.0, 8.0}, 0},
                                 {{0.0, 1.0, 9.0}, 0} };

  // setup tets
  int_el l_enVe[10][4] = { {4,5,6,7}, // tg 1
                           {4,5,6,7},
                           {0,0,0,0},
                           {4,5,6,7}, // tg 2
                           {0,1,2,3},
                           {0,1,2,3},
                           {0,1,2,3}, // ghost tg 2
                           {4,5,6,7},
                           {4,5,6,7},
                           {0,1,2,3} };

  // For plotting in Mathematica, you can use:
  //  Graphics3D[ {{Opacity[.3], Tetrahedron[]}, {PointSize[0.1], Point[{0.15, 0.15, 0.15}] }}, Axes -> True ]

  std::string l_recvNames[6] = {"1", "2", "3", "4", "5", "6"};
  std::string l_oDir = "/tmp";

  l_recv.init( TET4, 6, l_oDir, l_recvNames, l_recvCrds, 0.01, l_elLayout, l_enVe[0], l_veChars, 10 );

  std::vector< int_el > l_enRecv;
  l_recv.getEnRecv( l_enRecv );

  REQUIRE( l_enRecv.size() == 1 );
  REQUIRE( l_enRecv[0]     == 4 );
#endif
}

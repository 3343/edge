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
 * Unit tests of the dynamic memory allocations.
 **/
#include <catch.hpp>
#include "Dynamic.h"

TEST_CASE( "Dynamic: Flex data structure.", "[Dynamic][flex]" ) {
  edge::data::Dynamic l_dyn;

  /*
   *   i_nEn:                        5
   *   i_nSp:                        3
   *   i_bSize                       2
   *   i_spTypes:                   (0, 1, 0), (1, 1, 0), (0, 0, 1)
   *   i_spSizes:                    5,         3,         8
   *   .spType of i_enChars:        (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0), (1, 1, 1)
   *
   *                [-------base---------|------(0,1,0)------|------(1,1,0)-------|-------(0,0,1)------]
   *     addresses: [ b0, b1, b2, b3, b4 |  sp00, sp01, sp02 |   sp10,     sp11   |    sp20,    sp21   ]
   *     sizes:     [      5*2=10        |       3*5=15      |       2*3=6        |        2*8=16      ]
   *     pointers:         p00: b0             p10: nullptr        p20: nullptr          p30: sp20 
   *                       p01: b1             p11: sp00           p21: nullptr          p31: nullptr
   *                       p02: b2             p12: sp01           p22: sp10             p32: nullptr
   *                       p03: b3             p13: nullptr        p32: nullptr          p33: nullptr
   *                       p04: b4             p14: sp02           p33: sp11             p34: sp21
   *     total size: 10+15+6+16=47
   *
   *   resulting sizes:              2+8,      2+5,         2+5+3,      2,        2+5+3+8
   *   offset of returned pointers:    0,       10,            17,     27,             29,       47
   */
  double (**l_flex1)[20][8];

  unsigned short l_spTypes1[3] = {2, 6, 1};
  std::size_t    l_spSizes1[3]  = {5, 3, 8};
  struct{ unsigned short spType; } l_chars1[5] = {
    {1}, {2}, {6}, {4}, {7}
  };

  l_flex1 = l_dyn.flex< double [20][8] >( 5,
                                          3,
                                          2,
                                          l_spTypes1,
                                          l_spSizes1,
                                          l_chars1 );

  // check the base pointers
  REQUIRE( (l_flex1[ 1] - l_flex1[ 0]) == 2       );
  REQUIRE( (l_flex1[ 2] - l_flex1[ 1]) == 2       );
  REQUIRE( (l_flex1[ 3] - l_flex1[ 2]) == 2       );
  REQUIRE( (l_flex1[ 4] - l_flex1[ 3]) == 2       );
  REQUIRE( (l_flex1[ 6] - l_flex1[ 4]) == 2       );

  // check the first sparse type
  REQUIRE(  l_flex1[ 5]                == nullptr );
  REQUIRE( (l_flex1[ 7] - l_flex1[ 6]) == 5       );
  REQUIRE( (l_flex1[ 9] - l_flex1[ 7]) == 5       );
  REQUIRE(  l_flex1[ 8]                == nullptr );
  REQUIRE( (l_flex1[12] - l_flex1[ 9]) == 5       );

  // check the second sparse type
  REQUIRE(  l_flex1[10]                == nullptr );
  REQUIRE(  l_flex1[11]                == nullptr );
  REQUIRE( (l_flex1[14] - l_flex1[12]) == 3       );
  REQUIRE(  l_flex1[13]                == nullptr );
  REQUIRE( (l_flex1[15] - l_flex1[14]) == 3       );

  // check third sparse type
  REQUIRE( (l_flex1[19] - l_flex1[15]) == 8       );
  REQUIRE(  l_flex1[16]                == nullptr );
  REQUIRE(  l_flex1[17]                == nullptr );
  REQUIRE(  l_flex1[18]                == nullptr );
  REQUIRE( (l_flex1[20] - l_flex1[19]) == 8       );

  // size per entry: rows, columns, 8 byte per double
  std::size_t l_sizeEntry = 20 * 8 * 8;

  // check the base pointers
  REQUIRE( ((char*) l_flex1[ 1] - (char*) l_flex1[ 0]) == 2 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[ 2] - (char*) l_flex1[ 1]) == 2 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[ 3] - (char*) l_flex1[ 2]) == 2 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[ 4] - (char*) l_flex1[ 3]) == 2 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[ 6] - (char*) l_flex1[ 4]) == 2 * l_sizeEntry );

  // check the first sparse type
  REQUIRE( ((char*) l_flex1[ 7] - (char*) l_flex1[ 6]) == 5 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[ 9] - (char*) l_flex1[ 7]) == 5 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[12] - (char*) l_flex1[ 9]) == 5 * l_sizeEntry );

  // check the second sparse type
  REQUIRE( ((char*) l_flex1[14] - (char*) l_flex1[12]) == 3 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[15] - (char*) l_flex1[14]) == 3 * l_sizeEntry );

  // check third sparse type
  REQUIRE( ((char*) l_flex1[19] - (char*) l_flex1[15]) == 8 * l_sizeEntry );
  REQUIRE( ((char*) l_flex1[20] - (char*) l_flex1[19]) == 8 * l_sizeEntry );

  // assign data to trigger memory debuggers
  for( double* l_ch = (double*) l_flex1[0]; l_ch != (double*) l_flex1[20]; l_ch++ ) {
    *l_ch = 33.43;
  }
}

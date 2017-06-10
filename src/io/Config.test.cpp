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
 * Unit tests for EDGE's config.
 **/

#include <catch.hpp>

#include "Config.h"

TEST_CASE( "Config: Parsing of domains, given as XML-tree.", "[Config][domsParse]" ) {
  char l_xmlRup[] = "                          \
    <domain>                                   \
     <half_space>                              \
       <origin>                                \
         <x>16000</x>                          \
         <y>0</y>                              \
       </origin>                               \
       <normal>                                \
         <x>1</x>                              \
         <y>0</y>                              \
       </normal>                               \
     </half_space>                             \
                                               \
     <half_space>                              \
       <origin>                                \
         <x>19000</x>                          \
         <y>0</y>                              \
       </origin>                               \
       <normal>                                \
         <x>-1</x>                             \
         <y>0</y>                              \
       </normal>                               \
     </half_space>                             \
                                               \
     <values>                                  \
       <sn0>-120.0E6</sn0>                     \
       <ss0>                                   \
         <normal_strike>81.6E6</normal_strike> \
       </ss0>                                  \
     </values>                                 \
    </domain>                                  \
  ";

  // XML doc
  pugi::xml_document l_docRup;
  // result of pugixml's parsing
  pugi::xml_parse_result l_resRup = l_docRup.load_string( l_xmlRup );
  REQUIRE( l_resRup );

  // dom node
  pugi::xml_node l_domXml = l_docRup.root().child("domain");

  // name of data which are read
  std::vector< std::vector< std::string > > l_datN;

  // hierarchy of values
  l_datN.resize( 2 );
  l_datN[0].push_back( "values" ); l_datN[0].push_back( "sn0" );
  l_datN[1].push_back( "values" ); l_datN[1].push_back( "ss0" ); l_datN[1].push_back( "normal_strike" );

  // geometric description of the domain
  edge::linalg::Domain<
    double,
    2,
    edge::linalg::HalfSpace
   > l_domRes;

  // data associated to the domain which is read
  std::vector< double > l_dat;

  // parse XML-config
  edge::io::ConfigDoms< 2 >::parse( l_domXml, l_datN, edge::io::ConfigDoms< 2 >::F64, l_dat, l_domRes );

  // check data
  REQUIRE( l_dat.size() == 2 );
  REQUIRE( l_dat[0] == Approx( -120.0E6 ) );
  REQUIRE( l_dat[1] == Approx(   81.6E6 ) );

  // check domain
  double l_pt[2];

  l_pt[0] = 15999;
  l_pt[1] = -1000;
  REQUIRE( !l_domRes.inside(l_pt) );

  l_pt[0] = 16001;
  l_pt[1] = 99999;
  REQUIRE(  l_domRes.inside(l_pt) );

  l_pt[0] = 18999;
  l_pt[1] = -9999;
  REQUIRE(  l_domRes.inside(l_pt) );

  l_pt[0] = 19001;
  l_pt[1] =  9999;
  REQUIRE( !l_domRes.inside(l_pt) );
}
